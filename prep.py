#!/usr/bin/python

import argparse, glob, os, sys

import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ldsc'))
import munge_sumstats


def _allign_alleles(df):
    """Look for reversed alleles and inverts the z-score for one of them.

    Here, we take advantage of numpy's vectorized functions for performance.
    """
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    alleles = []
    for colname in ['A1_x', 'A2_x', 'A1_y', 'A2_y']:
        tmp = np.empty(len(df[colname]), dtype=int)
        for k, v in d.items():
            tmp[np.array(df[colname]) == k] = v
        alleles.append(tmp)
    reversed_alleles = np.logical_and(alleles[0] == alleles[3],
        alleles[1] == alleles[2])
    reversed_strand_flip_alleles = np.logical_and(alleles[0] == 3 - alleles[3],
        alleles[1] == 3 - alleles[2])
    to_flip = np.logical_or(reversed_alleles, reversed_strand_flip_alleles)
    df['Z_y'] *= -2 * to_flip + 1


def pre_function():
    munge_sumstats.parser.add_argument('sumstats1',
        help='the first sumstats file')
    munge_sumstats.parser.add_argument('sumstats2',
        help='the second sumstats file')
    munge_sumstats.parser.add_argument('--bimfile', default=None, type=str,
        required=True, help='bim filename. replace chrosome number with @ if \
            there are multiple.')
    munge_sumstats.parser.add_argument('--N1', default=None, type=int,
        help='N for sumstats1 if there is no N column')
    munge_sumstats.parser.add_argument('--N2', default=None, type=int,
        help='N for sumstats2 if there is no N column')

    args = munge_sumstats.parser.parse_args()
    # read in bim files
    all_bim_dfs = (pd.read_csv(f,
                               header=0,
                               names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                               delim_whitespace=True)
                   for f in glob.glob(args.bimfile.replace('@', '*')))
    bim = pd.concat(all_bim_dfs, ignore_index=True)

    # call munge_sumstats on the two files
    args.out = 'ldsc'  # we set this because it is required by munge_sumstats,
                        # but it is not used for our purposes.
    dfs = []
    for file, n in [(args.sumstats1, args.N1), (args.sumstats2, args.N2)]:
        print '=== MUNGING SUMSTATS FOR {} ==='.format(file)
        args.sumstats = file
        args.N = n
        dfs.append(munge_sumstats.munge_sumstats(args, p=False))

    # rename cols
    bim.rename(columns={'A1': 'A1_ref', 'A2': 'A2_ref'}, inplace=True)
    dfs[0].rename(columns={'A1': 'A1_x', 'A2': 'A2_x', 'N': 'N_x'},
        inplace=True)
    dfs[1].rename(columns={'A1': 'A1_y', 'A2': 'A2_y', 'N': 'N_y'},
        inplace=True)

    # take overlap between output and ref genotype files
    df = pd.merge(bim, dfs[1], on=['SNP']).merge(dfs[0], on=['SNP'])

    # flip sign of z-score for allele reversals
    _allign_alleles(df)

    # take maximum value of N
    df['N_x'] = np.maximum(df['N_x'], 0)
    df['N_y'] = np.maximum(df['N_y'], 0)
    return df
    # pickle.dump(df, open('df.p', 'wb'))


if __name__ == "__main__":
    # parse args
    df = pre_function()