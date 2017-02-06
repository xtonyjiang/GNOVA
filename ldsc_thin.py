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
    if args.bfile is not None:
        if args.l2 is None:
            raise ValueError('Must specify --l2 with --bfile.')
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
        bim_flag = args.bfile
        annot_flag = args.annot
        all_dfs = []
        for i in range(1, 23):
            print '=== COMPUTING LD SCORES FOR CHR {} ==='.format(i)
            args.bfile = bim_flag.replace('@', str(i))
            args.annot = annot_flag.replace('@', str(i))
            all_dfs.append(_ldscore(args, log))
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
parser.add_argument('--l2', default=False, action='store_true',
    help='Estimate l2. Compatible with both jackknife and non-jackknife.')
# Filtering / Data Management for LD Score
parser.add_argument('--extract', default=None, type=str,
    help='File with SNPs to include in LD Score estimation. '
    'The file should contain one SNP ID per row.')
parser.add_argument('--keep', default=None, type=str,
    help='File with individuals to include in LD Score estimation. '
    'The file should contain one individual ID per row.')
parser.add_argument('--ld-wind-snps', default=None, type=int,
    help='Specify the window size to be used for estimating LD Scores in units of '
    '# of SNPs. You can only specify one --ld-wind-* option.')
parser.add_argument('--ld-wind-kb', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'kilobase-pairs (kb). You can only specify one --ld-wind-* option.')
parser.add_argument('--ld-wind-cm', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'centiMorgans (cM). You can only specify one --ld-wind-* option.')
parser.add_argument('--print-snps', default=None, type=str,
    help='This flag tells LDSC to only print LD Scores for the SNPs listed '
    '(one ID per row) in PRINT_SNPS. The sum r^2 will still include SNPs not in '
    'PRINT_SNPs. This is useful for reducing the number of LD Scores that have to be '
    'read into memory when estimating h2 or rg.' )
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
parser.add_argument('--cts-breaks', default=None, type=str,
    help='Use this flag to specify names for the continuous variables cut into bins '
    'with --cts-bin. For each continuous variable, specify breaks as a comma-separated '
    'list of breakpoints, and separate the breakpoints for each variable with an x. '
    'For example, if binning on MAF and distance to gene (in kb), '
    'you might set --cts-breaks 0.1,0.25,0.4x10,100,1000 ')
parser.add_argument('--cts-names', default=None, type=str,
    help='Use this flag to specify names for the continuous variables cut into bins '
    'with --cts-bin. The argument to this flag should be a comma-separated list of '
    'names. For example, if binning on DAF and distance to gene, you might set '
    '--cts-bin DAF,DIST_TO_GENE ')
parser.add_argument('--per-allele', default=False, action='store_true',
    help='Setting this flag causes LDSC to compute per-allele LD Scores, '
    'i.e., \ell_j := \sum_k p_k(1-p_k)r^2_{jk}, where p_k denotes the MAF '
    'of SNP j. ')
parser.add_argument('--pq-exp', default=None, type=float,
    help='Setting this flag causes LDSC to compute LD Scores with the given scale factor, '
    'i.e., \ell_j := \sum_k (p_k(1-p_k))^a r^2_{jk}, where p_k denotes the MAF '
    'of SNP j and a is the argument to --pq-exp. ')
parser.add_argument('--no-print-annot', default=False, action='store_true',
    help='By defualt, seting --cts-bin or --cts-bin-add causes LDSC to print '
    'the resulting annot matrix. Setting --no-print-annot tells LDSC not '
    'to print the annot matrix. ')
parser.add_argument('--maf', default=None, type=float,
    help='Minor allele frequency lower bound. Default is MAF > 0.')
# Basic Flags for Working with Variance Components
parser.add_argument('--h2', default=None, type=str,
    help='Filename prefix for a .chisq file for one-phenotype LD Score regression. '
    'LDSC will automatically append .chisq or .chisq.gz to the filename prefix.'
    '--h2 requires at minimum also setting the --ref-ld and --w-ld flags.')
parser.add_argument('--rg', default=None, type=str,
    help='Comma-separated list of prefixes of .chisq filed for genetic correlation estimation.')
parser.add_argument('--ref-ld', default=None, type=str,
    help='Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD '
    'Score regression. '
    'LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix.')
parser.add_argument('--ref-ld-chr', default=None, type=str,
    help='Same as --ref-ld, but will automatically concatenate .l2.ldscore files split '
    'across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz '
    'to the filename prefix. If the filename prefix contains the symbol @, LDSC will '
    'replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome '
    'numbers to the end of the filename prefix.'
    'Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz'
    'Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz')
parser.add_argument('--w-ld', default=None, type=str,
    help='Filename prefix for file with LD Scores with sum r^2 taken over SNPs included '
    'in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz.')
parser.add_argument('--w-ld-chr', default=None, type=str,
    help='Same as --w-ld, but will read files split into 22 chromosomes in the same '
    'manner as --ref-ld-chr.')
parser.add_argument('--overlap-annot', default=False, action='store_true',
    help='This flag informs LDSC that the partitioned LD Scores were generates using an '
    'annot matrix with overlapping categories (i.e., not all row sums equal 1), '
    'and prevents LDSC from displaying output that is meaningless with overlapping categories.')
parser.add_argument('--no-intercept', action='store_true',
    help = 'If used with --h2, this constrains the LD Score regression intercept to equal '
    '1. If used with --rg, this constrains the LD Score regression intercepts for the h2 '
    'estimates to be one and the intercept for the genetic covariance estimate to be zero.')
parser.add_argument('--intercept-h2', action='store', default=None,
    help = 'Intercepts for constrained-intercept single-trait LD Score regression.')
parser.add_argument('--intercept-gencov', action='store', default=None,
    help = 'Intercepts for constrained-intercept cross-trait LD Score regression.'
    ' Must have same length as --rg. The first entry is ignored.')
parser.add_argument('--M', default=None, type=str,
    help='# of SNPs (if you don\'t want to use the .l2.M files that came with your .l2.ldscore.gz files)')
parser.add_argument('--two-step', default=None, type=float,
    help='Test statistic bound for use with the two-step estimator. Not compatible with --no-intercept and --constrain-intercept.')
parser.add_argument('--chisq-max', default=None, type=float,
    help='Max chi^2.')

# Flags for both LD Score estimation and h2/gencor estimation
# Flags you should almost never use
parser.add_argument('--chunk-size', default=50, type=int,
    help='Chunk size for LD Score calculation. Use the default.')
parser.add_argument('--pickle', default=False, action='store_true',
    help='Store .l2.ldscore files as pickles instead of gzipped tab-delimited text.')
parser.add_argument('--yes-really', default=False, action='store_true',
    help='Yes, I really want to compute whole-chromosome LD Score.')
parser.add_argument('--invert-anyway', default=False, action='store_true',
    help="Force LDSC to attempt to invert ill-conditioned matrices.")
parser.add_argument('--n-blocks', default=200, type=int,
    help='Number of block jackknife blocks.')
parser.add_argument('--not-M-5-50', default=False, action='store_true',
    help='This flag tells LDSC to use the .l2.M file instead of the .l2.M_5_50 file.')
parser.add_argument('--return-silly-things', default=False, action='store_true',
    help='Force ldsc to return silly genetic correlation estimates.')
parser.add_argument('--no-check-alleles', default=False, action='store_true',
    help='For rg estimation, skip checking whether the alleles match. This check is '
    'redundant for pairs of chisq files generated using munge_sumstats.py and the '
    'same argument to the --merge-alleles flag.')
parser.add_argument('--print-coefficients',default=False,action='store_true',
    help='when categories are overlapping, print coefficients as well as heritabilities.')
parser.add_argument('--frqfile', type=str,
    help='For use with --overlap-annot. Provides allele frequencies to prune to common '
    'snps if --not-M-5-50 is not set.')
parser.add_argument('--frqfile-chr', type=str,
    help='Prefix for --frqfile files split over chromosome.')

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
