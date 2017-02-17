#!/usr/bin/python
from __future__ import division
from itertools import product
import numpy as np
import pandas as pd
from sklearn import linear_model
from scipy.stats import norm
from numpy.linalg import inv

def calculate(gwas_snps, ld_scores, annots, N1, N2):
    ###   Clean up data   ###
    if annots is None:
        annot = ld_scores[['SNP']]
        annot['ALL_'] = 1
    else:
        annot = pd.concat(annots)

    ld_snps = set(ld_scores['SNP'])
    annot = annot.loc[annot['SNP'].isin(ld_snps)].reset_index(drop=True)

    ld_scores = ld_scores.drop(['CHR', 'BP', 'CM', 'MAF'], axis=1, errors='ignore').reset_index(drop=True)
    annot = annot.drop(['BP', 'SNP', 'CHR', 'CM'], axis=1, errors='ignore')
    gwas_snps.drop(['idx'], axis=1, errors='ignore', inplace=True)

    num_annotations = len(annot.columns)

    merged = pd.merge(gwas_snps,
                      pd.concat([ld_scores, annot], axis=1),
                      on=['SNP'])

    ld_score_all = merged.iloc[:,4]
    if num_annotations == 1:  # non-stratified analysis
        ld_scores = merged.iloc[:,4:5]
        annot = merged.iloc[:,5:6]

    else:  # we added in an all 1's column in prep step, so exclude that
        ld_scores = merged.iloc[:,5:4 + num_annotations]
        annot = merged.iloc[:,5 + num_annotations: 4 + 2 * num_annotations]
        num_annotations -= 1

    ###   Calculate genetic correlation   ###
    # Calculate S and W matrix
    P = annot.sum()
    p0 = len(ld_scores)

    S = np.empty([num_annotations, num_annotations])
    for i, j in product(range(num_annotations), range(num_annotations)):
        S[i][j] = np.sum(ld_scores[annot.iloc[:,i] == 1].iloc[:,j]) / (P[i] * P[j])

    W = np.empty([num_annotations, num_annotations])
    for i, j in product(range(num_annotations), range(num_annotations)):
        W[i][j] = np.sum((annot.iloc[:,i]==1) & (annot.iloc[:,j]==1)) / np.sum(annot.iloc[:,j] == 1)

    # Calculate heritability
    Z_x, Z_y = merged['Z_x'], merged['Z_y']

    if annots is None:
        h2_1 = np.array([p0 * (np.mean(Z_x ** 2) - 1) / (N1 * np.mean(ld_score_all))])
        h2_2 = np.array([p0 * (np.mean(Z_y ** 2) - 1) / (N2 * np.mean(ld_score_all))])
    else:
        tau1 = (np.mean((Z_x) ** 2) - 1)/(N1 * np.mean(ld_score_all))
        tau2 = (np.mean((Z_y) ** 2) - 1)/(N2 * np.mean(ld_score_all))
        w1 = 1 /(ld_score_all * (1 + N1 * tau1 * ld_score_all) ** 2)
        w2 = 1 /(ld_score_all * (1 + N2 * tau2 * ld_score_all) ** 2)
        w1[(w1 < 0) | (w1 == np.inf) | (w1 == -np.inf)] = 0
        w2[(w2 < 0) | (w2 == np.inf) | (w2 == -np.inf)] = 0
        m1 = linear_model.LinearRegression().fit(ld_scores, pd.DataFrame((Z_x) ** 2), sample_weight=w1)
        m2 = linear_model.LinearRegression().fit(ld_scores, pd.DataFrame((Z_y) ** 2), sample_weight=w2)
        h2_1 = np.dot(W, m1.coef_.T * pd.DataFrame(P) / N1)
        h2_2 = np.dot(W, m2.coef_.T * pd.DataFrame(P) / N2)

    # Calculate sample overlap correction
    if annots is None:
        w1 = 1 + N1 * (h2_1 * ld_score_all / len(ld_score_all))
        w2 = 1 + N2 * (h2_2 * ld_score_all / len(ld_score_all))
    else:
        w1 = 1 + p0 * (np.mean(Z_x ** 2) - 1) / np.mean(ld_score_all) * ld_score_all / len(ld_score_all)
        w2 = 1 + p0 * (np.mean(Z_y ** 2) - 1) / np.mean(ld_score_all) * ld_score_all / len(ld_score_all)

    w3 = np.mean(Z_x * Z_y) * ld_score_all
    w = 1 / (w1 * w2 + w3 * w3)
    m = linear_model.LinearRegression().fit(pd.DataFrame(ld_score_all), pd.DataFrame(Z_x * Z_y), sample_weight=w)
    corr_pheno = m.intercept_[0]

    # Calculate Jackknife variance estimate
    nblock = 200
    q_block = np.empty([num_annotations, nblock])

    for i in range(num_annotations):
        df_x = Z_x[annot.iloc[:,i] == 1]
        df_y = Z_y[annot.iloc[:,i] == 1]
        tot = np.dot(df_x, df_y)
        for j, (b_x, b_y) in enumerate(zip(np.array_split(df_x, nblock), np.array_split(df_y, nblock))):
            q_block[i][j] = (tot - np.dot(b_x, b_y)) / ((len(df_x) - len(b_x) - corr_pheno) * ((N1 * N2) ** 0.5))

    q = np.mean(q_block, axis=1)
    cov_q = np.cov(q_block, bias=True) * (nblock - 1)

    # rho
    rho = W.dot(inv(S)).dot(q)
    rho_corrected = W.dot(inv(S)).dot(q - corr_pheno / ((N1 * N2) ** 0.5))

    # covariance of rho
    cov_rho = W.dot(inv(S)).dot(cov_q).dot(inv(S)).dot(W.T)

    # genetic correlation
    corr = rho / ((h2_1 * h2_2) ** 0.5).T
    corr_corrected = rho_corrected / ((h2_1 * h2_2) ** 0.5).T

    # p-value
    p_value = norm.sf(abs(rho / (cov_rho.diagonal() ** 0.5))) * 2
    p_value_corrected = norm.sf(abs(rho_corrected / (cov_rho.diagonal() ** 0.5))) * 2

    return {
                'rho': rho,
                'rho_corrected': rho_corrected,
                'p_value': p_value,
                'p_value_corrected': p_value_corrected,
                'var_rho': cov_rho.diagonal(),
                'corr': corr[0],
                'corr_corrected': corr_corrected[0],
                'h2_1': h2_1.T[0],
                'h2_2': h2_2.T[0],
                'p': P,
                'p0': p0
           }
