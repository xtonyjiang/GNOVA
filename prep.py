#!/usr/bin/python

import argparse, glob, os, re, sys

import numpy as np
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ldsc'))
import munge_sumstats


def allign_alleles(df):
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


def matched_or_reversed(df):
    """Returns boolean array signifying whether rows have matched or reversed
    alleles.
    """
    d = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    a = []  # array of alleles
    for colname in ['A1_x', 'A2_x', 'A1_y', 'A2_y']:
        tmp = np.empty(len(df[colname]), dtype=int)
        for k, v in d.items():
            tmp[np.array(df[colname]) == k] = v
        a.append(tmp)
    matched_alleles = (((a[0] == a[2]) & (a[1] == a[3])) |
        ((a[0] == 3 - a[2]) | (a[1] == 3 - a[3])))
    reversed_alleles = (((a[0] == a[3]) & (a[1] == a[2])) |
        ((a[0] == 3 - a[0]) | (a[1] == 3 - a[2])))
    return matched_alleles | reversed_alleles


def pre_function(args):
    # read in bim files
    if '@' in args.bimfile:
        bim_dfs = (pd.read_csv(args.bimfile.replace('@', str(i)),
                               header=None,
                               names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                               delim_whitespace=True)
                   for i in range(1, 23)
                   if os.path.isfile(args.bimfile.replace('@', str(i))))
        bim = pd.concat(bim_dfs, ignore_index=True)
    else:
        bim = pd.read_csv(args.bimfile,
                          header=None,
                          names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                          delim_whitespace=True)

    # call munge_sumstats on the two sumstats files
    dfs = []
    ms_args = munge_sumstats.parser.parse_args([])
    ms_args.out = 'ldsc'  # we set this because it is required by munge_sumstats
                          # but it is not used for our purposes.
    for file, n, should_munge in [(args.sumstats1, args.N1, args.munge1),
                                  (args.sumstats2, args.N2, args.munge2)]:
        if should_munge:
            print '=== MUNGING SUMSTATS FOR {} ==='.format(file)
            ms_args.sumstats = file
            ms_args.N = n
            dfs.append(munge_sumstats.munge_sumstats(ms_args, p=False))
        else:
            dfs.append(pd.read_csv(file, delim_whitespace=True))

    # rename cols
    bim.rename(columns={'A1': 'A1_ref', 'A2': 'A2_ref'}, inplace=True)
    dfs[0].rename(columns={'A1': 'A1_x', 'A2': 'A2_x', 'N': 'N_x', 'Z': 'Z_x'},
        inplace=True)
    dfs[1].rename(columns={'A1': 'A1_y', 'A2': 'A2_y', 'N': 'N_y', 'Z': 'Z_y'},
        inplace=True)

    # take overlap between output and ref genotype files
    df = pd.merge(bim, dfs[1], on=['SNP']).merge(dfs[0], on=['SNP'])

    # flip sign of z-score for allele reversals
    allign_alleles(df)
    df = df[matched_or_reversed(df)]

    # take maximum value of N
    df['N_x'] = np.maximum(df['N_x'], 0)
    df['N_y'] = np.maximum(df['N_y'], 0)

    return df[['CHR', 'SNP', 'N_x', 'Z_x', 'N_y', 'Z_y']]


if __name__ == "__main__":
    # parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('sumstats1',
        help='the first sumstats file')
    parser.add_argument('sumstats2',
        help='the second sumstats file')
    parser.add_argument('--bimfile', default=None, type=str, required=True,
        help='bim filename. replace chromosome number with @ if \
            there are multiple.')
    parser.add_argument('--N1', default=None, type=int,
        help='N for sumstats1 if there is no N column')
    parser.add_argument('--N2', default=None, type=int,
        help='N for sumstats2 if there is no N column')
    df = pre_function(parser.parse_args())
