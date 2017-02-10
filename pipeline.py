import argparse
from prep import prep
from ldsc_thin import ldscore
from calculate import calculate
import pandas as pd

def pipeline(args):
    pd.options.mode.chained_assignment = None
    gwas_snps, N1, N2, annots = prep(args.bfile,
                                args.annot,
                                args.sumstats1,
                                args.already_munged1,
                                args.N1,
                                args.sumstats2,
                                args.already_munged2,
                                args.N2)
    ld_scores = ldscore(args.bfile, annots, gwas_snps)
    results = calculate(gwas_snps, ld_scores, annots, N1, N2)
    print(results)


parser = argparse.ArgumentParser()

parser.add_argument('sumstats1',
    help='The first sumstats file.')
parser.add_argument('sumstats2',
    help='The second sumstats file.')

parser.add_argument('--already_munged1', action='store_true',
    help='Denotes that the first sumstats file has not already been run '
         'through munge_sumstats.')
parser.add_argument('--already_munged2', action='store_true',
    help='Denotes that the second sumstats file has not already been run '
         'through munge_sumstats.')

parser.add_argument('--N1', default=None, type=int,
    help='N for sumstats1 if there is no N column.')
parser.add_argument('--N2', default=None, type=int,
    help='N for sumstats2 if there is no N column.')

parser.add_argument('--bfile', required=True, type=str,
    help='Prefix for Plink .bed/.bim/.fam file.')
parser.add_argument('--annot', type=str,
    help='Filename prefix for annotation file for partitioned LD Score estimation. '
    'LDSC will automatically append .annot or .annot.gz to the filename prefix. '
    'See docs/file_formats_ld for a definition of the .annot format.')

# parser.add_argument('--out', required=True, type=str,
#     help='Location to output results.')


if __name__ == '__main__':
    pipeline(parser.parse_args())
