#!/usr/bin/env python
'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

LDSC is a command line tool for estimating
    1. LD Score
    2. heritability / partitioned heritability
    3. genetic covariance / correlation

'''
from __future__ import division
import ldscore.ldscore as ld
import ldscore.parse as ps
import prep as pr
import numpy as np
import pandas as pd
from subprocess import call
from itertools import product
import time, sys, traceback, argparse, glob

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.drop_duplicates(subset='A')
except TypeError:
    raise ImportError('LDSC requires pandas version > 0.15.2')

__version__ = '1.0.0'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* LD Score Regression (LDSC)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane\n"
MASTHEAD += "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f


def _remove_dtype(x):
    '''Removes dtype: float64 and dtype: int64 from pandas printouts'''
    x = str(x)
    x = x.replace('\ndtype: int64', '')
    x = x.replace('\ndtype: float64', '')
    return x


class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module

    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'wb')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.

        '''
        print msg

def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = ps.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        print f(c, len(x.IDList))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            print f(c, len_merged_list)
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list

def loj_bim(filter_df, array):
    r = filter_df.columns[1]
    l = array.IDList.columns[0]
    merge_df = filter_df.iloc[:,[1]]
    merge_df.loc[:,'keep'] = True
    z = pd.merge(array.IDList, merge_df, how='left', left_on=l, right_on=r, sort=False)
    ii = z['keep'] == True
    return np.nonzero(ii)[0]

def __filter_bim__(filter_df, array):
    merged_list = loj_bim(filter_df, array)
    len_merged_list = len(merged_list)
    if len_merged_list > 0:
        c = 'After merging, {0} SNPs remain'
        print c.format(len_merged_list)
    else:
        error_msg = 'No SNPs retained for analysis'
        raise ValueError(error_msg)

    return merged_list

def annot_sort_key(s):
    '''For use with --cts-bin. Fixes weird pandas crosstab column order.'''
    if type(s) == tuple:
        s = [x.split('_')[0] for x in s]
        s = map(lambda x: float(x) if x != 'min' else -float('inf'), s)
    else:  # type(s) = str:
        s = s.split('_')[0]
        if s == 'min':
            s = float('-inf')
        else:
            s = float(s)

    return s

def subset_annot_file(a_df, GWAS_df, kept_cols):
    GWAS_df.loc[:,'idx'] = pd.Series(range(len(GWAS_df.SNP.values)))
    a_df = pd.merge(a_df, GWAS_df, how="right", on=['SNP'])
    a_df = a_df.sort_values(['idx'], ascending=True)
    a_df.drop('idx', axis=1, inplace=True)
    a_df.rename(columns={'CHR_x':'CHR', 'BP_x':'BP', 'CM_x':'CM'}, inplace=True)
    a_df = a_df.iloc[:,0:kept_cols]
    return a_df

def remove_brackets(x):
    return x.replace('[', '').replace(']', '').strip()

def _ldscore(args, log):
    '''
    Wrapper function for estimating l1, l1^2, l2 and l4 (+ optionally standard errors) from
    reference panel genotypes.

    Annot format is
    chr snp bp cm <annotations>

    '''

    #args = set_defaults(args)
    if args.n_blocks <= 1:
        raise ValueError('--n-blocks must be an integer > 1.')
    if args.bfile:
        snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
        ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
        array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile
    GWASsnps_df = args.GWASsnps
    # read bim/snp
    array_snps = snp_obj(snp_file)
    # snp list
    m = len(array_snps.IDList)
    log.log('Read list of {m} SNPs from {f}'.format(m=m, f=snp_file))
    if args.annot is not None:  # read --annot
        try:
            annots = [pd.read_csv(f, delim_whitespace=True) for f in glob.glob(args.annot.replace('@', '*'))]
            annot = ps.AnnotFile(pd.concat(annots, ignore_index=True))
            n_annot, ma = len(annot.df.columns) - 4, len(annot.df)
            log.log("Read {A} annotations for {M} SNPs from {f}".format(f=args.annot,
                A=n_annot, M=ma))
            annot_matrix = np.array(annot.df.iloc[:,4:])
            annot_colnames = annot.df.columns[4:]
            keep_snps = None
            #take only annot SNPs in intersect
            kept_cols = len(annot.df.columns)
            annot.df = subset_annot_file(annot.df, GWASsnps_df, kept_cols)
            if np.any(annot.df.SNP.values != GWASsnps_df.SNP.values):
                raise ValueError('The .annot file must contain all SNPs in the study intersect in the same'+\
                    ' order as the .bim file.')
        except Exception:
            log.log('Error parsing .annot file')
            raise
    else:
        annot_matrix, annot_colnames, keep_snps = None, None, None,
        n_annot = 1

    keep_snps = __filter_bim__(GWASsnps_df, array_snps)


    # read fam
    array_indivs = ind_obj(ind_file)
    n = len(array_indivs.IDList)
    log.log('Read list of {n} individuals from {f}'.format(n=n, f=ind_file))
    # read keep_indivs
    if args.keep:
        keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
    else:
        keep_indivs = None

    # read genotype array
    log.log('Reading genotypes from {fname}'.format(fname=array_file))
    geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=args.maf)

    #determine block widths

    max_dist = args.ld_wind_cm
    coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

    block_left = ld.getBlockLefts(coords, max_dist)

    scale_suffix = ''

    lN = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_matrix)
    col_prefix = "L2"; file_suffix = "l2"

    if n_annot == 1:
        ldscore_colnames = [col_prefix+scale_suffix]
    else:
        ldscore_colnames =  [y+col_prefix+scale_suffix for y in annot_colnames]

    # print .ldscore. Output columns: CHR, BP, RS, [LD Scores]
    out_fname = args.out + '.' + file_suffix + '.ldscore'
    new_colnames = geno_array.colnames + ldscore_colnames
    df = pd.DataFrame.from_records(np.c_[geno_array.df, lN])
    df.columns = new_colnames
    l2_suffix = '.gz'
    log.log("Writing LD Scores for {N} SNPs to {f}.gz".format(f=out_fname, N=len(df)))
    df.drop(['CM','MAF'], axis=1)

    # print annot matrix
    if (args.cts_bin is not None) and not args.no_print_annot:
        out_fname_annot = args.out + '.annot'
        new_colnames = geno_array.colnames + ldscore_colnames
        annot_df = pd.DataFrame(np.c_[geno_array.df, annot_matrix])
        annot_df.columns = new_colnames
        del annot_df['MAF']
        log.log("Writing annot matrix produced by --cts-bin to {F}".format(F=out_fname+'.gz'))
        #annot.df

    # print LD Score summary
    pd.set_option('display.max_rows', 200)
    log.log('\nSummary of LD Scores in {F}'.format(F=out_fname+l2_suffix))
    t = df.ix[:,4:].describe()
    log.log( t.ix[1:,:] )

    np.seterr(divide='ignore', invalid='ignore')  # print NaN instead of weird errors
    # print correlation matrix including all LD Scores and sample MAF
    log.log('')
    log.log('MAF/LD Score Correlation Matrix')
    log.log( df.ix[:,4:].corr() )

    # print condition number
    if n_annot > 1: # condition number of a column vector w/ nonzero var is trivially one
        log.log('\nLD Score Matrix Condition Number')
        cond_num = np.linalg.cond(df.ix[:,5:])
        log.log(remove_brackets(str(np.matrix(cond_num))))
        if cond_num > 10000:
            log.log('WARNING: ill-conditioned LD Score Matrix!')

    # summarize annot matrix if there is one
    if annot_matrix is not None:
        # covariance matrix
        x = pd.DataFrame(annot_matrix, columns=annot_colnames)
        log.log('\nAnnotation Correlation Matrix')
        log.log( x.corr() )

        # column sums
        log.log('\nAnnotation Matrix Column Sums')
        log.log(_remove_dtype(x.sum(axis=0)))

        # row sums
        log.log('\nSummary of Annotation Matrix Row Sums')
        row_sums = x.sum(axis=1).describe()
        log.log(_remove_dtype(row_sums))

    np.seterr(divide='raise', invalid='raise')
    return df

def ldscore(args, log):
    df = None
    if '@' in args.bfile:
        b_flag = args.bfile
        annot_flag = args.annot
        all_dfs = []
        for i in range(1, 23):
            print '=== COMPUTING LD SCORES FOR CHROMOSOME {} ==='.format(i)
            args.bfile = b_flag.replace('@', str(i))
            args.annot = annot_flag.replace('@', str(i))
            all_dfs.append(_ldscore(args, log))
        args.bfile = bim_flag
        args.annot= annot_flag
        df = pd.concat(all_dfs)
    else:
        df = _ldscore(args, log)

    numeric = df._get_numeric_data()
    numeric[numeric < 0] = 0
    return df

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='ldsc', type=str,
    help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
    'default output filename prefix.')
# Basic LD Score Estimation Flags'
parser.add_argument('--bfile', default=None, type=str,
    help='Prefix for Plink .bed/.bim/.fam file')
# Filtering / Data Management for LD Score
parser.add_argument('--extract', default=None, type=str,
    help='File with SNPs to include in LD Score estimation. '
    'The file should contain one SNP ID per row.')
parser.add_argument('--keep', default=None, type=str,
    help='File with individuals to include in LD Score estimation. '
    'The file should contain one individual ID per row.')
parser.add_argument('--ld-wind-cm', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'centiMorgans (cM). You can only specify one --ld-wind-* option.')
# Fancy LD Score Estimation Flags
parser.add_argument('--annot', default=None, type=str,
    help='Filename prefix for annotation file for partitioned LD Score estimation. '
    'LDSC will automatically append .annot or .annot.gz to the filename prefix. '
    'See docs/file_formats_ld for a definition of the .annot format.')
parser.add_argument('--cts-bin', default=None, type=str,
    help='This flag tells LDSC to compute partitioned LD Scores, where the partition '
    'is defined by cutting one or several continuous variable[s] into bins. '
    'The argument to this flag should be the name of a single file or a comma-separated '
    'list of files. The file format is two columns, with SNP IDs in the first column '
    'and the continuous variable in the second column. ')
parser.add_argument('--no-print-annot', default=False, action='store_true',
    help='By defualt, seting --cts-bin or --cts-bin-add causes LDSC to print '
    'the resulting annot matrix. Setting --no-print-annot tells LDSC not '
    'to print the annot matrix. ')
parser.add_argument('--maf', default=None, type=float,
    help='Minor allele frequency lower bound. Default is MAF > 0.')

# Flags for both LD Score estimation and h2/gencor estimation
# Flags you should almost never use
parser.add_argument('--chunk-size', default=50, type=int,
    help='Chunk size for LD Score calculation. Use the default.')
parser.add_argument('--n-blocks', default=200, type=int,
    help='Number of block jackknife blocks.')

if __name__ == '__main__':
    args = parser.parse_args()
    args.bfile = "/net/zhao/ql68/GeneticCorrelation/1000G/EUR/eur_chr1_SNPmaf5"
    args.annot = "/net/zhao/ql68/GeneticCorrelation/realdata/CD_UC/All_SNPmaf5/tmpFiles/SNPmaf5.1.annot.gz"
    args.sumstats1 = "/net/zhao/ql68/GeneticCorrelation/sumstats/CD.sumstats.gz"
    args.sumstats2 = "/net/zhao/ql68/GeneticCorrelation/sumstats/UC.sumstats.gz"
    args.bimfile = "/net/zhao/ql68/GeneticCorrelation/1000G/EUR/eur_chr1_SNPmaf5.bim"
    args.munge1 = False
    args.munge2 = False
    args.N1 = 27726
    args.N2 = 28738
    args.out = "temp"

    args.ld_wind_cm = 1
    log = Logger(args.out + '.log')
    try:
        args.GWASsnps = pr.pre_function(args)
        LD_matrix = ldscore(args, log)
        print(LD_matrix)

    # bad flags
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
