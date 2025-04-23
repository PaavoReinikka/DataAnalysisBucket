# imports for the global null experiments
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats import multitest
from statsmodels.stats.multitest import fdrcorrection_twostage, fdrcorrection
import copy

import sys
sys.path.insert(1, "../code/utils")
from celldata import *
from experiment_utils import *

import multiprocessing
from joblib import Parallel, delayed


# READY TO USE
def global_null_experimentV1(X,n_iter, alpha=0.05, filter=None, equal_var=True, keep_k=0):
    # Filtering is applied to the data before the t-test, only one filter is applied
    # negative keep_k means that all p-values are kept, zero means none are kept, and positive means k are kept
    
    if filter is None:
        filter = lambda x: x
        
    FD_bh = np.zeros(n_iter,dtype=int)
    FD_by = np.zeros(n_iter,dtype=int)
    FD_bh_2step = np.zeros(n_iter,dtype=int)
    FD_bky_2step = np.zeros(n_iter,dtype=int)
    FD_bonf = np.zeros(n_iter,dtype=int)
    FD_unadj = np.zeros(n_iter,dtype=int)
    
    flag=False
    if keep_k > 50:
        flag=True
        pvalues = -1.0*np.ones((n_iter, keep_k))
    elif keep_k > 0 and keep_k <= 50:
        pvalues = np.zeros((n_iter, keep_k))
    else:
        pvalues=None
        
    counts = np.zeros(n_iter,dtype=int)
    m0_estimates_bh = np.zeros(n_iter)
    m0_estimates_bky = np.zeros(n_iter)
    
    for i in range(n_iter):
        x_perm = X.flatten()
        x_perm = np.random.permutation(x_perm).reshape(X.shape)
        x_perm = filter(x_perm)
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        
        pvals = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var).pvalue
        pvals = np.sort(pvals)
        counts[i] = len(pvals)
        
        if pvalues is not None:
            k = np.min([len(pvals), pvalues.shape[1]])
            pvalues[i,:k] = pvals[:k]
        
        FD_bonf[i] = multitest.multipletests(pvals, alpha=alpha, method='bonferroni', is_sorted=True)[0].sum()
        FD_unadj[i] = (pvals < alpha).sum()
        FD_bh[i] = multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)[0].sum()
        FD_by[i] = multitest.multipletests(pvals, alpha=alpha, method='fdr_by', is_sorted=True)[0].sum()
        
        m0_estimate_bh, alpha_adjusted = estimate_m0_bh(pvals, alpha=alpha, is_sorted=True)
        FD_bh_2step[i] = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=True)[0].sum()
        m0_estimates_bh[i] = m0_estimate_bh
        
        m0_estimate_bky, alpha_adjusted = estimate_m0_bky(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)
        FD_bky_2step[i] = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=True)[0].sum()
        m0_estimates_bky[i] = m0_estimate_bky
        #FD_bh_2step[i] = fdrcorrection_twostage(pvals, alpha=alpha, method='bh')[0].sum()
        #FD_bky_2step[i] = fdrcorrection_twostage(pvals, alpha=alpha, method='bky', iter=True)[0].sum()
        
    if flag:
        # flatten and remove -1.0
        pvalues = pvalues.flatten()
        pvalues = pvalues[pvalues != -1.0]    
    
    # return everything in a dictionary
    return {'FD_bh': FD_bh,
            'FD_by': FD_by,
            'FD_bh_2step': FD_bh_2step,
            'FD_bky_2step': FD_bky_2step,
            'FD_bonf': FD_bonf,
            'FD_unadj': FD_unadj,
            'pvalues': pvalues,
            'counts': counts,
            'm0_estimates_bh': m0_estimates_bh,
            'm0_estimates_bky': m0_estimates_bky}

# READY TO USE
def global_null_experimentV2(X, GAPS, n_iter=1e4, alpha=0.05, k=10, filter=None, equal_var=True, keep_k=0):
    # 2 filters are applied to the data before the t-test, one based on data and one based on reliability (gaps)
    
    if filter is None:
        filter = lambda x: (x,np.arange(x.shape[0]))
        
    FD_bh = np.zeros(n_iter,dtype=int)
    FD_by = np.zeros(n_iter,dtype=int)
    FD_bh_2step = np.zeros(n_iter,dtype=int)
    FD_bky_2step = np.zeros(n_iter,dtype=int)
    FD_bonf = np.zeros(n_iter,dtype=int)
    FD_unadj = np.zeros(n_iter,dtype=int)
    
    flag=False
    if keep_k > 50:
        flag=True
        pvalues = -1.0*np.ones((n_iter, keep_k))
    elif keep_k > 0 and keep_k <= 50:
        pvalues = np.zeros((n_iter, keep_k))
    else:
        pvalues=None
    
    counts = np.zeros(n_iter,dtype=int)
    m0_estimates_bh = np.zeros(n_iter)
    m0_estimates_bky = np.zeros(n_iter)
    
    for i in range(n_iter):
        x_perm = X.flatten()
        gaps_perm = GAPS.flatten()
        random_perm = np.random.permutation(len(x_perm))
        x_perm = x_perm[random_perm].reshape(X.shape)
        gaps_perm = gaps_perm[random_perm].reshape(GAPS.shape)
        
        x_perm, inds = filter(x_perm)
        gaps_perm = gaps_perm[inds]
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        gaps_perm = compile_pairwise_comparisons_legacy(gaps_perm)
        inds = reliability(gaps_perm[:,:10], gaps_perm[:,10:], 0.8, [0,8,64,128])[0]
        x_perm = x_perm[inds]
        
        pvals = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var).pvalue
        pvals = np.sort(pvals)
        counts[i] = len(pvals)
        
        if pvalues is not None:
            k = np.min([len(pvals), pvalues.shape[1]])
            pvalues[i,:k] = pvals[:k]
        
        FD_bonf[i] = multitest.multipletests(pvals, alpha=alpha, method='bonferroni', is_sorted=True)[0].sum()
        FD_unadj[i] = (pvals < alpha).sum()
        FD_bh[i] = multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)[0].sum()
        FD_by[i] = multitest.multipletests(pvals, alpha=alpha, method='fdr_by', is_sorted=True)[0].sum()
        
        m0_estimate_bh, alpha_adjusted = estimate_m0_bh(pvals, alpha=alpha, is_sorted=True)
        FD_bh_2step[i] = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=True)[0].sum()
        m0_estimates_bh[i] = m0_estimate_bh
        
        m0_estimate_bky, alpha_adjusted = estimate_m0_bky(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)
        FD_bky_2step[i] = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=True)[0].sum()
        m0_estimates_bky[i] = m0_estimate_bky
        
        #FD_bh_2step[i] = fdrcorrection_twostage(pvals, alpha=alpha, method='bh')[0].sum()
        #FD_bky_2step[i] = fdrcorrection_twostage(pvals, alpha=alpha, method='bky')[0].sum()
        
    if flag:
        # flatten and remove -1.0
        pvalues = pvalues.flatten()
        pvalues = pvalues[pvalues != -1.0]
    
    return {'FD_bh': FD_bh, \
            'FD_by': FD_by, \
            'FD_bh_2step': FD_bh_2step, \
            'FD_bky_2step': FD_bky_2step, \
            'FD_bonf': FD_bonf, \
            'FD_unadj': FD_unadj, \
            'pvalues': pvalues, \
            'counts': counts, \
            'm0_estimates_bh': m0_estimates_bh, \
            'm0_estimates_bky': m0_estimates_bky}

# READY TO USE
def global_null_experimentV3(X, GAPS, n_iter=1e4, alpha=0.05, k=10, filter=None, equal_var=True, keep_k=0):
    # 3 filters are applied to the data before the t-test, one based on data, one based on reliability (gaps), and one based on the Fisher's exact test
    # Fisher's test is approximated here by a very conservative 10-to-0 or 0-to-10 test

    if filter is None:
        filter = lambda x: (x,np.arange(x.shape[0]))
        
    FD_bh = np.zeros(n_iter,dtype=int)
    FD_by = np.zeros(n_iter,dtype=int)
    FD_bh_2step = np.zeros(n_iter,dtype=int)
    FD_bky_2step = np.zeros(n_iter,dtype=int)
    FD_bonf = np.zeros(n_iter,dtype=int)
    FD_unadj = np.zeros(n_iter,dtype=int)
    FD_fisher = np.zeros(n_iter,dtype=int)
    
    flag=False
    if keep_k > 50:
        flag=True
        pvalues = -1.0*np.ones((n_iter, keep_k))
    elif keep_k > 0 and keep_k <= 50:
        pvalues = np.zeros((n_iter, keep_k))
    else:
        pvalues=None
    
    counts = np.zeros(n_iter,dtype=int)
    m0_estimates_bh = np.zeros(n_iter)
    m0_estimates_bky = np.zeros(n_iter)
    
    for i in range(n_iter):
        x_perm = X.flatten()
        gaps_perm = GAPS.flatten()
        random_perm = np.random.permutation(len(x_perm))
        x_perm = x_perm[random_perm].reshape(X.shape)
        gaps_perm = gaps_perm[random_perm].reshape(GAPS.shape)
        
        x_perm, inds = filter(x_perm)
        gaps_perm = gaps_perm[inds]
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        gaps_perm = compile_pairwise_comparisons_legacy(gaps_perm)
        inds, counts1, counts2 = reliability(gaps_perm[:,:10], gaps_perm[:,10:], 0.8, [0,8,64,128])
        inds = np.logical_and(counts1 ==10, counts2 ==0) | np.logical_and(counts1 == 0, counts2 == 10)
        fd_fisher = inds.sum()
        inds = ~inds
        x_perm = x_perm[inds]
        
        pvals = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var).pvalue
        pvals = np.sort(pvals)
        counts[i] = len(pvals)
        
        if pvalues is not None:
            k = np.min([len(pvals), pvalues.shape[1]])
            pvalues[i,:k] = pvals[:k]
        
        FD_fisher[i] = fd_fisher
        FD_bonf[i] = multitest.multipletests(pvals, alpha=alpha, method='bonferroni', is_sorted=True)[0].sum()
        FD_unadj[i] = (pvals < alpha).sum()
        FD_bh[i] = multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)[0].sum()
        FD_by[i] = multitest.multipletests(pvals, alpha=alpha, method='fdr_by', is_sorted=True)[0].sum()
        
        m0_estimate_bh, alpha_adjusted = estimate_m0_bh(pvals, alpha=alpha, is_sorted=True)
        FD_bh_2step[i] = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=True)[0].sum()
        m0_estimates_bh[i] = m0_estimate_bh
        
        m0_estimate_bky, alpha_adjusted = estimate_m0_bky(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)
        FD_bky_2step[i] = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=True)[0].sum()
        m0_estimates_bky[i] = m0_estimate_bky
        
        #FD_bh_2step[i] = fdrcorrection_twostage(pvals, alpha=alpha, method='bh')[0].sum()
        #FD_bky_2step[i] = fdrcorrection_twostage(pvals, alpha=alpha, method='bky')[0].sum()
        
    if flag:
        # flatten and remove -1.0
        pvalues = pvalues.flatten()
        pvalues = pvalues[pvalues != -1.0]
    
    return {'FD_bh': FD_bh, \
            'FD_by': FD_by, \
            'FD_bh_2step': FD_bh_2step, \
            'FD_bky_2step': FD_bky_2step, \
            'FD_bonf': FD_bonf, \
            'FD_unadj': FD_unadj, \
            'FD_fisher': FD_fisher, \
            'pvalues': pvalues, \
            'counts': counts, \
            'm0_estimates_bh': m0_estimates_bh, \
            'm0_estimates_bky': m0_estimates_bky}
    
    
# READY TO USE
def global_null_resampling_experimentV1(X, n_iter, n_iter_inner=1000, alpha=0.05, filter=None, equal_var=True, full_resampling=False):

    if filter is None:
        filter = lambda x: x
        
    counts = np.zeros(n_iter,dtype=int)
    FD_wy = np.zeros(n_iter,dtype=int)
    FD_reiner_bh = np.zeros(n_iter,dtype=int)
    FD_reiner_by = np.zeros(n_iter,dtype=int)
    
    for i in range(n_iter):
        x_perm = X.flatten()
        x_perm = np.random.permutation(x_perm).reshape(X.shape)
        x_perm = filter(x_perm)
        
        # full resampling means that we permute across all samples and groups
        if full_resampling:
            # resampling procedure has to use same ordering of pairwise comparisons as the following t-tests
            _, __, resampled_stats, max_stats = westfall_young_resampled_pvalues_full(x_perm, n_iter=n_iter_inner, equal_var=equal_var)
            # compile pairwise comparisons using same methods as is done within westfall_young_resampled_pvalues
            x_perm = compile_pairwise_comparisons_legacy(x_perm)
        else:
            # compile pairwise comparisons before resampling
            x_perm = compile_pairwise_comparisons_legacy(x_perm)
            _, __, resampled_stats, max_stats = westfall_young_resampled_pvalues(x_perm, n_iter=n_iter_inner, equal_var=equal_var)        
        
        #pval_wy = np.sort(min_pvals)[int(alpha*n_iter_inner)]
        stat_wy = np.sort(max_stats)[::-1][int(alpha*n_iter_inner)]
        result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var)
        pvals = result.pvalue
        statistics = result.statistic
        #adjusted_pvals_reiner, adjusted_pvals_marginal = resampling_adjusted_pvalues(pvals, resampled_pvals)
        adjusted_pvals_reiner = resampling_adjusted_pvalues_from_stats(statistics, resampled_stats)
        
        
        counts[i] = len(pvals)
        #FD_wy[i] = (pvals < pval_wy).sum()
        FD_wy[i] = (np.abs(statistics) > stat_wy).sum()
        FD_reiner_bh[i] = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_bh')[0].sum()
        FD_reiner_by[i] = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_by')[0].sum()
        
    return {'FD_wy': FD_wy, \
            'FD_reiner_bh': FD_reiner_bh, \
            'FD_reiner_by': FD_reiner_by, \
            'counts': counts}

# READY TO USE
def global_null_resampling_experimentV2(X, GAPS, n_iter=1e4, n_iter_inner=1000, alpha=0.05, k=10, filter=None, equal_var=True):
    
    if filter is None:
        filter = lambda x: (x,np.arange(x.shape[0]))
            
    counts = np.zeros(n_iter,dtype=int)
    FD_wy = np.zeros(n_iter,dtype=int)
    FD_reiner_bh = np.zeros(n_iter,dtype=int)
    FD_reiner_by = np.zeros(n_iter,dtype=int)
        
    for i in range(n_iter):
        x_perm = X.flatten()
        gaps_perm = GAPS.flatten()
        random_perm = np.random.permutation(len(x_perm))
        x_perm = x_perm[random_perm].reshape(X.shape)
        gaps_perm = gaps_perm[random_perm].reshape(GAPS.shape)
        
        x_perm, inds = filter(x_perm)
        gaps_perm = gaps_perm[inds]
        gaps_perm = compile_pairwise_comparisons_legacy(gaps_perm)            
        inds = reliability(gaps_perm[:,:10], gaps_perm[:,10:], 0.8, [0,8,64,128])[0]
            
        # compile pairwise comparisons before resampling, and using only the reliable tests
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        x_perm = x_perm[inds]
        #resampled_pvals, min_pvals, _ = westfall_young_resampled_pvalues(x_perm, n_iter=n_iter_inner)
        _, __, resampled_stats, max_stat = westfall_young_resampled_pvalues(x_perm, n_iter=n_iter_inner)
            
        #pval_wy = np.sort(min_pvals)[int(alpha*n_iter_inner)]
        stat_wy = np.sort(max_stat)[::-1][int(alpha*n_iter_inner)]
        result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var)
        pvals = result.pvalue
        statistics = result.statistic
        #adjusted_pvals_reiner, adjusted_pvals_marginal = resampling_adjusted_pvalues(pvals, resampled_pvals)
        adjusted_pvals_reiner = resampling_adjusted_pvalues_from_stats(statistics, resampled_stats)
         
        counts[i] = len(pvals)
        #FD_wy[i] = (pvals < pval_wy).sum()
        FD_wy[i] = (np.abs(statistics) > stat_wy).sum()
        FD_reiner_bh[i] = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_bh')[0].sum()
        FD_reiner_by[i] = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_by')[0].sum()
            
    return {'FD_wy': FD_wy, \
            'FD_reiner_bh': FD_reiner_bh, \
            'FD_reiner_by': FD_reiner_by, \
            'counts': counts}
            
# READY TO USE
def global_null_resampling_experimentV3(X, GAPS, n_iter=1e4, n_iter_inner=1000, alpha=0.05, k=10, filter=None, equal_var=True):
    # like the previous, but now we also include the Fisher's exact test
    if filter is None:
        filter = lambda x: (x,np.arange(x.shape[0]))
        
    counts = np.zeros(n_iter,dtype=int)
    FD_wy = np.zeros(n_iter,dtype=int)
    FD_reiner_bh = np.zeros(n_iter,dtype=int)
    FD_reiner_by = np.zeros(n_iter,dtype=int)
    FD_fisher = np.zeros(n_iter,dtype=int)
    
    for i in range(n_iter):
        x_perm = X.flatten()
        gaps_perm = GAPS.flatten()
        random_perm = np.random.permutation(len(x_perm))
        x_perm = x_perm[random_perm].reshape(X.shape)
        gaps_perm = gaps_perm[random_perm].reshape(GAPS.shape)
        
        x_perm, inds = filter(x_perm)
        gaps_perm = gaps_perm[inds]
        gaps_perm = compile_pairwise_comparisons_legacy(gaps_perm)
        inds, counts1, counts2 = reliability(gaps_perm[:,:10], gaps_perm[:,10:], 0.8, [0,8,64,128])
        inds = np.logical_and(counts1 ==10, counts2 ==0) | np.logical_and(counts1 == 0, counts2 == 10)
        fd_fisher = inds.sum()
        inds = ~inds
        
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        x_perm = x_perm[inds]
        #resampled_pvals, min_pvals, _ = westfall_young_resampled_pvalues(x_perm, n_iter=n_iter_inner)
        _, __, resampled_stats, max_stat = westfall_young_resampled_pvalues(x_perm, n_iter=n_iter_inner)
        
        #pval_wy = np.sort(min_pvals)[int(alpha*n_iter_inner)]
        stat_wy = np.sort(max_stat)[::-1][int(alpha*n_iter_inner)]
        result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var)
        pvals = result.pvalue
        statistics = result.statistic
        #adjusted_pvals_reiner, adjusted_pvals_marginal = resampling_adjusted_pvalues(pvals, resampled_pvals)
        adjusted_pvals_reiner = resampling_adjusted_pvalues_from_stats(statistics, resampled_stats)
        
        counts[i] = len(pvals)
        #FD_wy[i] = (pvals < pval_wy).sum()
        FD_wy[i] = (np.abs(statistics) > stat_wy).sum()
        FD_reiner_bh[i] = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_bh')[0].sum()
        FD_reiner_by[i] = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_by')[0].sum()
        FD_fisher[i] = fd_fisher
        
    return {'FD_wy': FD_wy, \
            'FD_reiner_bh': FD_reiner_bh, \
            'FD_reiner_by': FD_reiner_by, \
            'FD_fisher': FD_fisher, \
            'counts': counts}
    
    