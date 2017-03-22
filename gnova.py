#!/usr/bin/env python
import argparse, collections, os.path, sys
from prep import prep
from ldsc_thin import ldscore
from calculate import calculate
import pandas as pd


# returns whether the parent directory of path exists
def parent_dir_exists(path):
    return os.path.exists(os.path.abspath(os.path.join(path, os.pardir)))


def pipeline(args):
    pd.options.mode.chained_assignment = None

    # Sanity check args
    if args.save_ld is not None and args.use_ld is not None:
        raise ValueError('Both the --save-ld and --use-ld flags were set. '
                         'Please use only one of them.')
    if args.save_ld is not None:
        if not parent_dir_exists(args.save_ld):
            raise ValueError('--save-ld flag points to an invalid path.')
    if args.use_ld is not None:
        if not os.path.exists(args.use_ld + '.csv.gz'):
            raise ValueError('--use-ld flag points to an invalid path.')
    if not parent_dir_exists(args.out):
        raise ValueError('--out flag points to an invalid path.')

    print('Preparing files for analysis...')
    gwas_snps, N1, N2, annots = prep(args.bfile,
                                     args.annot,
                                     args.sumstats1,
                                     args.sumstats2)
    if args.N1 is not None:
        N1 = args.N1
    if args.N2 is not None:
        N2 = args.N2
    if args.use_ld is not None:
        print('Loading LD scores from {}'.format(args.use_ld))
        ld_scores = pd.read_csv(args.use_ld + '.csv.gz', sep=' ')
    else:
        print('Calculating LD scores...')
        ld_scores = ldscore(args.bfile, annots, gwas_snps, args.save_ld)
    print('Calculating correlation...')
    results = calculate(gwas_snps, ld_scores, annots, N1, N2)
    results_to_print = collections.OrderedDict(
        [('rho', results['rho']),
         ('rho_corrected', results['rho_corrected']),
         ('pvalue', results['p_value']),
         ('pvalue_corrected', results['p_value_corrected']),
         ('corr', results['corr']),
         ('corr_corrected', results['corr_corrected']),
         ('h2_1', results['h2_1']),
         ('h2_2', results['h2_2']),
         ('p', results['p']),
         ('p0', results['p0'])
        ]
    )
    out = pd.DataFrame(results_to_print)
    print('Final results:\n{}\n'.format(out))
    print('\nView ldsc.log for verbose output.')
    out.insert(0, 'annot_name', out.index)
    out.to_csv(args.out, sep=' ', index=False)


parser = argparse.ArgumentParser()

parser.add_argument('sumstats1',
    help='The first sumstats file.')
parser.add_argument('sumstats2',
    help='The second sumstats file.')

parser.add_argument('--bfile', required=True, type=str,
    help='Prefix for Plink .bed/.bim/.fam file.')
parser.add_argument('--annot', type=str,
    help='Filename prefix for annotation file for partitioned LD Score estimation. '
    'LDSC will automatically append .annot or .annot.gz to the filename prefix. '
    'See docs/file_formats_ld for a definition of the .annot format.')
parser.add_argument('--N1', type=int,
    help='N of the sumstats1 file. If not provided, this value will be inferred '
    'from the sumstats1 arg.')
parser.add_argument('--N2', type=int,
    help='N of the sumstats2 file. If not provided, this value will be inferred '
    'from the sumstats2 arg.')

parser.add_argument('--out', required=True, type=str,
    help='Location to output results.')
parser.add_argument('--save-ld', type=str,
    help='Prefix of the location to save LD score calculations as gzipped .csv '
         'files. If not set, then then no intermediate calculations will be saved.')
parser.add_argument('--use-ld', type=str,
    help='Prefix of the location to load LD score calculations from.')

if __name__ == '__main__':
    if sys.version_info[0] != 2:
        print('ERROR: GNOVA does not run on Python 3. Please run it on Python 2.7.x.')
        sys.exit(1)
    pipeline(parser.parse_args())
