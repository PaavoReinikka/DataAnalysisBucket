import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats import multitest
import copy

#NOTE: This code is not meant to be run as a script. It is meant to be imported into
#      another script. You need to import the below utilities before importing this

#import sys
#sys.path.insert(1, "../utils")
#sys.path.insert(1, "../experiment_code")
#from experiment_utils import *
#from celldata import *
#from effects import *
#from n_effective import *
#from experiment_dumps import *


def FDR_experimentV1(X,n_iter, alpha=0.05, filter=None, equal_var=True, \
                    mask_model = lambda x: np.zeros((x.shape[0],4)), effect_model=lambda x,y: x,\
                        scale_model = None, TAG='UNKNOWN'):
    # data X should have no significant effects (either X is residuals or X has sigificant effects removed)

    # filter is a function that takes in data and returns a subset of the data
    if filter is None:
        filter = lambda x: (x, np.arange(x.shape[0]))
    if scale_model is None:
        scale_model = lambda x: x
    elif scale_model is not None:
        assert callable(scale_model), "scale_model must be a function"
        assert np.all(np.isclose(X[:,:10].mean(axis=1),X[:,10:20].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,10:20].mean(axis=1),X[:,20:30].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,20:30].mean(axis=1),X[:,30:].mean(axis=1))), "groups must be centered"
                                   
    FN_bh = np.zeros(n_iter)
    FP_bh = np.zeros(n_iter)
    rejections_bh = np.zeros(n_iter,dtype=int)
    
    FN_by = np.zeros(n_iter)
    FP_by = np.zeros(n_iter)
    rejections_by = np.zeros(n_iter,dtype=int)
    
    FN_bonf = np.zeros(n_iter,dtype=int)
    FP_bonf = np.zeros(n_iter,dtype=int)
    rejections_bonf = np.zeros(n_iter,dtype=int)
    
    FN_bh_2step = np.zeros(n_iter,dtype=int)
    FP_bh_2step = np.zeros(n_iter,dtype=int)
    rejections_bh_2step = np.zeros(n_iter,dtype=int)

    FN_bky_2step = np.zeros(n_iter,dtype=int)
    FP_bky_2step = np.zeros(n_iter,dtype=int)
    rejections_bky_2step = np.zeros(n_iter,dtype=int)
    
    counts = np.zeros(n_iter,dtype=int)
    n_effects = np.zeros(n_iter,dtype=int)
    
    m0_estimates_bh = np.zeros(n_iter)
    m0_estimates_bky = np.zeros(n_iter)
    
    for i in range(n_iter):
        rand_inds = np.random.permutation(X.shape[1])
        x_perm = X[:,rand_inds]
        
        # make sure that mask function is randomizing the locations of the effects
        x_perm = scale_model(x_perm)
        mask = mask_model(x_perm)
        x_perm = effect_model(x_perm, mask) 
        
        x_perm, inds = filter(x_perm)
        mask = mask[inds]
        
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        true_rejections = compile_pairwise_effect_indicators_legacy(mask)
        n_effects[i] = np.sum(true_rejections)
        
        pvals = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var).pvalue
        counts[i] = len(pvals)
        
        h = multitest.multipletests(pvals, alpha=alpha, method='fdr_bh')[0]
        rejections_bh[i] = np.sum(h)
        FN_bh[i] = false_negatives(h,true_rejections)
        FP_bh[i] = false_positives(h,true_rejections)
        
        h = multitest.multipletests(pvals, alpha=alpha, method='fdr_by')[0]
        rejections_by[i] = np.sum(h)
        FN_by[i] = false_negatives(h,true_rejections)
        FP_by[i] = false_positives(h,true_rejections)
        
        h = multitest.multipletests(pvals, alpha=alpha, method='bonferroni')[0]
        rejections_bonf[i] = np.sum(h)
        FN_bonf[i] = false_negatives(h,true_rejections)
        FP_bonf[i] = false_positives(h,true_rejections)
        
        m0_estimate_bh, alpha_adjusted = estimate_m0_bh(pvals, alpha=alpha, is_sorted=False)
        h = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=False)[0]
        rejections_bh_2step[i] = np.sum(h)
        FN_bh_2step[i] = false_negatives(h,true_rejections)
        FP_bh_2step[i] = false_positives(h,true_rejections)        
        
        m0_estimate_bky, alpha_adjusted = estimate_m0_bky(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
        h = multitest.multipletests(pvals, alpha=alpha_adjusted, method='fdr_bh', is_sorted=False)[0]
        rejections_bky_2step[i] = np.sum(h)
        FN_bky_2step[i] = false_negatives(h,true_rejections)
        FP_bky_2step[i] = false_positives(h,true_rejections)
        
        m0_estimates_bh[i] = m0_estimate_bh
        m0_estimates_bky[i] = m0_estimate_bky
        
    
    return {'TAG':TAG,'BH':{'FN': FN_bh, 'FP': FP_bh, 'rejections': rejections_bh}, \
            'BY':{'FN': FN_by, 'FP': FP_by, 'rejections': rejections_by}, \
            'BONF':{'FN': FN_bonf, 'FP': FP_bonf, 'rejections': rejections_bonf}, \
            'BH_2STEP':{'FN': FN_bh_2step, 'FP': FP_bh_2step, 'rejections': rejections_bh_2step}, \
            'BKY_2STEP':{'FN': FN_bky_2step, 'FP': FP_bky_2step, 'rejections': rejections_bky_2step}, \
            'n_effects': n_effects, 'n_tests': counts, 'm0_estimates_bh': m0_estimates_bh, 'm0_estimates_bky': m0_estimates_bky}
 

def FDR_experimentV2(X,n_iter, alpha=0.05, filter=None, equal_var=True, \
                    mask_model = lambda x: np.zeros((x.shape[0],4)), effect_model=lambda x,y: x,\
                        scale_model = None, method = 'fdr_bh', n_effective_method = 'n_effective_full', th=0.1, TAG='UNKNOWN'):
    # data X should have no significant effects (either X is residuals or X has sigificant effects removed)

    # filter is a function that takes in data and returns a subset of the data
    if filter is None:
        filter = lambda x: (x, np.arange(x.shape[0]))
    if scale_model is None:
        scale_model = lambda x: x
    elif scale_model is not None:
        assert callable(scale_model), "scale_model must be a function"
        assert np.all(np.isclose(X[:,:10].mean(axis=1),X[:,10:20].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,10:20].mean(axis=1),X[:,20:30].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,20:30].mean(axis=1),X[:,30:].mean(axis=1))), "groups must be centered"
    
    FN_fdr_spectral = np.zeros(n_iter)
    FP_fdr_spectral = np.zeros(n_iter)
    rejections_fdr_spectral = np.zeros(n_iter,dtype=int)
    
    FN_fdr_cluster = np.zeros(n_iter)
    FP_fdr_cluster = np.zeros(n_iter)
    rejections_fdr_cluster = np.zeros(n_iter,dtype=int)
                                   
    FN_bonf_spectral = np.zeros(n_iter,dtype=int)
    FP_bonf_spectral = np.zeros(n_iter,dtype=int)
    rejections_bonf_spectral = np.zeros(n_iter,dtype=int)
    
    FN_bonf_cluster = np.zeros(n_iter,dtype=int)
    FP_bonf_cluster = np.zeros(n_iter,dtype=int)
    rejections_bonf_cluster = np.zeros(n_iter,dtype=int)
    
    counts = np.zeros(n_iter,dtype=int)
    n_effects = np.zeros(n_iter,dtype=int)
    M_eff_spectral = np.zeros(n_iter)
    M_eff_cluster = np.zeros(n_iter)
    
    for i in range(n_iter):
        rand_inds = np.random.permutation(X.shape[1])
        x_perm = X[:,rand_inds]
        
        # make sure that mask function is randomizing the locations of the effects
        x_perm = scale_model(x_perm)
        mask = mask_model(x_perm)
        x_perm = effect_model(x_perm, mask)
        
        x_perm, inds = filter(x_perm)
        mask = mask[inds]
        
        if n_effective_method == 'n_effective_full':
            cor = np.corrcoef(x_perm)
            M_eff_spectral[i] = n_eff_spectral_fast(cor, factor=6)
            M_eff_cluster[i] = n_eff_cluster_fast(cor, threshold=th, factor=6)
        
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        true_rejections = compile_pairwise_effect_indicators_legacy(mask)
        n_effects[i] = np.sum(true_rejections)
        
        if n_effective_method == 'n_effective':
            cor = np.corrcoef(x_perm)
            M_eff_spectral[i] = n_eff_spectral_fast(cor, factor=1)
            M_eff_cluster[i] = n_eff_cluster_fast(cor, threshold=th, factor=1)
        
        pvals = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var).pvalue
        counts[i] = len(pvals)
        
        h = (M_eff_spectral[i]*pvals < alpha)        
        rejections_bonf_spectral[i] = np.sum(h)
        FN_bonf_spectral[i] = false_negatives(h,true_rejections)
        FP_bonf_spectral[i] = false_positives(h,true_rejections)
        
        h = (M_eff_cluster[i]*pvals < alpha)
        rejections_bonf_cluster[i] = np.sum(h)
        FN_bonf_cluster[i] = false_negatives(h,true_rejections)
        FP_bonf_cluster[i] = false_positives(h,true_rejections)
         
        h = multitest.multipletests(pvals, alpha=alpha*len(pvals)/M_eff_spectral[i], method=method)[0]
        rejections_fdr_spectral[i] = np.sum(h)
        FN_fdr_spectral[i] = false_negatives(h,true_rejections)
        FP_fdr_spectral[i] = false_positives(h,true_rejections)
        
        h = multitest.multipletests(pvals, alpha=alpha*len(pvals)/M_eff_cluster[i], method=method)[0]
        rejections_fdr_cluster[i] = np.sum(h)
        FN_fdr_cluster[i] = false_negatives(h,true_rejections)
        FP_fdr_cluster[i] = false_positives(h,true_rejections)
        
    return {'TAG':TAG,'FDR_SPECTRAL':{'FN': FN_fdr_spectral, 'FP': FP_fdr_spectral, 'rejections': rejections_fdr_spectral}, \
            'FDR_CLUSTER':{'FN': FN_fdr_cluster, 'FP': FP_fdr_cluster, 'rejections': rejections_fdr_cluster}, \
            'BONF_SPECTRAL':{'FN': FN_bonf_spectral, 'FP': FP_bonf_spectral, 'rejections': rejections_bonf_spectral}, \
            'BONF_CLUSTER':{'FN': FN_bonf_cluster, 'FP': FP_bonf_cluster, 'rejections': rejections_bonf_cluster}, \
            'n_effects': n_effects, 'n_tests': counts, 'M_eff_spectral': M_eff_spectral, 'M_eff_cluster': M_eff_cluster}        


def FDR_experimentV3(X,n_iter, alpha=0.05, filter=None, equal_var=True, \
                    mask_model = lambda x: np.zeros((x.shape[0],4)), effect_model=lambda x,y: x,\
                        scale_model = None, method = 'fdr_bh', combine_method = 'min', th=0.1, TAG='UNKNOWN'):
    # data X should have no significant effects (either X is residuals or X has sigificant effects removed)

    # filter is a function that takes in data and returns a subset of the data
    if filter is None:
        filter = lambda x: (x, np.arange(x.shape[0]))
    if scale_model is None:
        scale_model = lambda x: x
    elif scale_model is not None:
        assert callable(scale_model), "scale_model must be a function"
        assert np.all(np.isclose(X[:,:10].mean(axis=1),X[:,10:20].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,10:20].mean(axis=1),X[:,20:30].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,20:30].mean(axis=1),X[:,30:].mean(axis=1))), "groups must be centered"

    FN = np.zeros(n_iter)
    FP = np.zeros(n_iter)
    rejections = np.zeros(n_iter,dtype=int)
    
    counts = np.zeros(n_iter,dtype=int)
    n_effects = np.zeros(n_iter,dtype=int)
    M_eff = np.zeros(n_iter)
    
    pvalues = np.zeros((n_iter,X.shape[0]*6))
    
    for i in range(n_iter):
        rand_inds = np.random.permutation(X.shape[1])
        x_perm = X[:,rand_inds]
        
        # make sure that mask function is randomizing the locations of the effects
        x_perm = scale_model(x_perm)
        mask = mask_model(x_perm)
        x_perm = effect_model(x_perm, mask)
        
        x_perm, inds = filter(x_perm)
        mask = mask[inds]
        
        x_perm = compile_pairwise_comparisons_legacy(x_perm)
        true_rejections = compile_pairwise_effect_indicators_legacy(mask)
        n_effects[i] = np.sum(true_rejections)
        
        cor = np.corrcoef(x_perm)
        cor = n_eff_cluster_fast(cor, threshold = th, method = 'complete', return_labels = True, factor = 1)
        M_eff[i] = cor[0]
        labels = cor[1]
                
        pvals = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var).pvalue
        counts[i] = len(pvals)
        pvals = combine_pvalues(pvals[:,np.newaxis], labels, method = combine_method)
        pvalues[i,:] = pvals.flatten()
        
        h = multitest.multipletests(pvals, alpha=alpha, method=method)[0]
        rejections[i] = np.sum(h)
        FN[i] = false_negatives(h,true_rejections)
        FP[i] = false_positives(h,true_rejections)
        
    return {'TAG':TAG,method.upper():{'FN': FN, 'FP': FP, 'rejections': rejections}, \
            'n_effects': n_effects, 'n_tests': counts, 'M_eff': M_eff, 'pvalues': pvalues}
        
        
    
def FDR_resampling_experimentV1(X,n_iter, n_iter_inner=1000, alpha=0.05, filter=None, equal_var=True,\
    mask_model = lambda x: np.zeros((x.shape[0],4)), effect_model=lambda x, y: x,\
                        scale_model = None, full_resampling=False, TAG='UNKNOWN'):
    # data X should have no significant effects (either X is residuals or X has sigificant effects removed)
    
    if filter is None:
        filter = lambda x: x
    if scale_model is None:
        scale_model = lambda x: x
    elif scale_model is not None:
        assert callable(scale_model), "scale_model must be a function"
        assert np.all(np.isclose(X[:,:10].mean(axis=1),X[:,10:20].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,10:20].mean(axis=1),X[:,20:30].mean(axis=1))), "groups must be centered"
        assert np.all(np.isclose(X[:,20:30].mean(axis=1),X[:,30:].mean(axis=1))), "groups must be centered"
        
    FN_wy = np.zeros(n_iter)
    FP_wy = np.zeros(n_iter)
    
    FN_reiner_bh = np.zeros(n_iter)
    FP_reiner_bh = np.zeros(n_iter)
    
    FN_reiner_by = np.zeros(n_iter)
    FP_reiner_by = np.zeros(n_iter)
    
    rejections_wy = np.zeros(n_iter,dtype=int)
    rejections_reiner_bh = np.zeros(n_iter,dtype=int)
    rejections_reiner_by = np.zeros(n_iter,dtype=int)
    
    counts = np.zeros(n_iter,dtype=int)
    n_effects = np.zeros(n_iter,dtype=int)
    
    for i in range(n_iter):
        # permute columns of X to form null distribution
        rand_inds = np.random.permutation(X.shape[1])
        x_perm = X[:,rand_inds]

        # make sure that mask function is randomizing the locations of the effects
        x_perm = scale_model(x_perm)
        mask = mask_model(x_perm)
        x_perm = effect_model(x_perm, mask)

        # filter out rows from data and mask
        x_perm, inds = filter(x_perm)
        mask = mask[inds,:]
        
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
        
        true_rejections = compile_pairwise_effect_indicators_legacy(mask)
        n_effects[i] = np.sum(true_rejections)
        
        #pval_wy = np.sort(min_pvals)[int(alpha*n_iter_inner)]
        stat_wy = np.sort(max_stats)[::-1][int(alpha*n_iter_inner)]
        result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], axis=1, equal_var=equal_var)
        pvals = result.pvalue
        statistics = result.statistic
        #adjusted_pvals_reiner, adjusted_pvals_marginal = resampling_adjusted_pvalues(pvals, resampled_pvals)
        adjusted_pvals_reiner = resampling_adjusted_pvalues_from_stats(statistics, resampled_stats)
        
        counts[i] = len(pvals)
        #h_wy[i] = (pvals < pval_wy)
        h_wy = (np.abs(statistics) > stat_wy)
        h_reiner_bh = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_bh', is_sorted=False)[0]
        h_reiner_by = multitest.multipletests(adjusted_pvals_reiner, alpha=alpha, method='fdr_by', is_sorted=False)[0]
        
        rejections_wy[i] = np.sum(h_wy)
        rejections_reiner_bh[i] = np.sum(h_reiner_bh)
        rejections_reiner_by[i] = np.sum(h_reiner_by)
        
        FN_wy[i] = false_negatives(h_wy,true_rejections)
        FP_wy[i] = false_positives(h_wy,true_rejections)

        FN_reiner_bh[i] = false_negatives(h_reiner_bh,true_rejections)
        FP_reiner_bh[i] = false_positives(h_reiner_bh,true_rejections)
        
        FN_reiner_by[i] = false_negatives(h_reiner_by,true_rejections)
        FP_reiner_by[i] = false_positives(h_reiner_by,true_rejections)
        
    # return everything in a dictionary
    return {'TAG':TAG,'WY':{'FN': FN_wy, 'FP': FP_wy, 'rejections': rejections_wy}, \
            'REINER_BH':{'FN': FN_reiner_bh, 'FP': FP_reiner_bh, 'rejections': rejections_reiner_bh},\
            'REINER_BY':{'FN': FN_reiner_by, 'FP': FP_reiner_by, 'rejections': rejections_reiner_by},\
            'n_effects': n_effects, 'n_tests': counts}
        
    
