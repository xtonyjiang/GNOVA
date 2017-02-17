#!/usr/bin/env python
import argparse, collections, cPickle
from prep import prep
from ldsc_thin import ldscore
from calculate import calculate
import pandas as pd

def pipeline(args):
    pd.options.mode.chained_assignment = None
    gwas_snps, N1, N2, annots = prep(args.bfile,
                                     args.annot,
                                     args.sumstats1,
                                     args.sumstats2)
    if args.use_ld is not None:
        ld_scores = cPickle.load(open(args.use_ld, 'rb'))
    else:
        ld_scores = ldscore(args.bfile, annots, gwas_snps, args.save_ld)
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

parser.add_argument('--out', required=True, type=str,
    help='Location to output results.')
parser.add_argument('--save-ld', type=str,
    help='Location to save LD score calculations as pickle files to. If not '
         'set, then then no intermediate calculations will be saved.')
parser.add_argument('--use-ld', type=str,
    help='Location to load LD score calculations from.')

if __name__ == '__main__':
    pipeline(parser.parse_args())
