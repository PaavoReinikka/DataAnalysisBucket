import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats import multitest
from statsmodels.stats.multitest import fdrcorrection_twostage, fdrcorrection

import pickle
import sys
sys.path.insert(1, "../code/utils")
from celldata import *

 
            
            
def variance_filter(val=0, method='prop'):
    # method: 'prop' or 'value'
    if val == 0:
        return lambda X: X
    elif method == 'value':
        return lambda X: X[X.var(axis=1) > val]
    elif method == 'prop':
        # return 1-val proportion of the data
        # by dropping the lowest variance rows
        return lambda X: X[X.var(axis=1) > np.quantile(X.var(axis=1), val)]
    else:
        raise ValueError("method must be 'prop' or 'value'")

def variance_filter_inds(val=0, method='prop'):
    # method: 'prop' or 'value'
    if val==0:
        return lambda X: (X, np.arange(X.shape[0]))
    elif method == 'value':
        return lambda X: (X[X.var(axis=1) > val], np.where(X.var(axis=1) > val)[0])
    elif method == 'prop':
        # return 1-val proportion of the data
        # by dropping the lowest variance rows
        return lambda X: (X[X.var(axis=1) > np.quantile(X.var(axis=1), val)], np.where(X.var(axis=1) > np.quantile(X.var(axis=1), val))[0])
    else:
        raise ValueError("method must be 'prop' or 'value'")
        
        
def mean_filter(val=0, method='prop'):
    # method: 'prop' or 'value'
    if val == 0:
        return lambda X: X
    elif method == 'value':
        return lambda X: X[np.abs(X.mean(axis=1)) > val]
    elif method == 'prop':
        # return 1-val proportion of the data
        # by dropping the lowest variance rows
        return lambda X: X[np.abs(X.mean(axis=1)) > np.quantile(np.abs(X.mean(axis=1)), val)]
    else:
        raise ValueError("method must be 'prop' or 'value'")

def mean_filter_inds(val=0, method='prop'):
    # method: 'prop' or 'value'
    if val == 0:
        return lambda X: (X, np.arange(X.shape[0]))
    elif method == 'value':
        return lambda X: (X[np.abs(X.mean(axis=1)) > val], np.where(np.abs(X.mean(axis=1)) > val)[0])
    elif method == 'prop':
        # return 1-val proportion of the data
        # by dropping the lowest variance rows
        return lambda X: (X[np.abs(X.mean(axis=1)) > np.quantile(np.abs(X.mean(axis=1)), val)], np.where(np.abs(X.mean(axis=1)) > np.quantile(np.abs(X.mean(axis=1)), val))[0])
    else:
        raise ValueError("method must be 'prop' or 'value'")


def westfall_young_resampled_pvalues_full(x, n_iter=1e4, equal_var=True):
    # westfall young resampling procedure over all samples and groups
    x_perm = compile_pairwise_comparisons_legacy(x)
    pvals = np.zeros((n_iter, x_perm.shape[0]))
    statistics = np.zeros((n_iter, x_perm.shape[0]))
    min_pvals = np.zeros(n_iter)
    max_stats = np.zeros(n_iter)
    result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], equal_var=equal_var, axis=1)
    pvals[0,:] = result.pvalue
    statistics[0,:] = result.statistic
    
    for i in range(1,n_iter):
        rand_inds = np.random.permutation(x.shape[1])
        x_perm = compile_pairwise_comparisons_legacy(x[:,rand_inds])
        result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], equal_var=equal_var, axis=1)
        pvals[i,:] = result.pvalue
        statistics[i,:] = result.statistic
        min_pvals[i] = np.min(pvals[i,:])
        max_stats[i] = np.max(np.abs(statistics[i,:]))
    
    return pvals, min_pvals, statistics, max_stats

def westfall_young_resampled_pvalues(x, n_iter=1e4, equal_var=True):
    # westfall young resampling procedure over all samples and two groups
    pvals = np.zeros((n_iter, x.shape[0]))
    statistics = np.zeros((n_iter, x.shape[0]))
    min_pvals = np.zeros(n_iter)
    max_stats = np.zeros(n_iter)
    result = stats.ttest_ind(x[:,:10], x[:,10:], equal_var=equal_var, axis=1)
    pvals[0,:] = result.pvalue
    statistics[0,:] = result.statistic
    
    for i in range(1,n_iter):
        rand_inds = np.random.permutation(x.shape[1])
        x_perm = x[:,rand_inds]
        result = stats.ttest_ind(x_perm[:,:10], x_perm[:,10:], equal_var=equal_var, axis=1)
        pvals[i,:] = result.pvalue
        statistics[i,:] = result.statistic
        min_pvals[i] = np.min(pvals[i,:])
        max_stats[i] = np.max(np.abs(statistics[i,:]))
    
    return pvals, min_pvals, statistics, max_stats



def resampling_adjusted_pvalues(pvalues_original, resampled_pvalues):
    # pvalues_original: p-values from the original data
    # resampled_pvalues: p-values from resampled data
    # The method by Reiner et al. in "Resampling-based FDR control in microarray analysis", (2002)
    # The marginal method does not collapse the resampled p-values
    n_iter = resampled_pvalues.shape[0]
    pvalues_adjusted_reiner = np.zeros(pvalues_original.shape[0])
    #pvalues_adjusted_marginal = np.zeros(pvalues_original.shape[0])
    
    for i in range(pvalues_original.shape[0]):
        s=0
        for j in range(n_iter):
            s += (resampled_pvalues[j,:] <= pvalues_original[i]).mean()
        pvalues_adjusted_reiner[i] = s/n_iter
        #pvalues_adjusted_marginal[i] = (resampled_pvalues[:,i] <= pvalues_original[i]).mean()
        
    return pvalues_adjusted_reiner#, pvalues_adjusted_marginal

def resampling_adjusted_pvalues_from_stats(stas_original, resampled_stats):
    # like the above, but now doing it for the statistics (2-sided)
    n_iter = resampled_stats.shape[0]
    pvalues_adjusted_reiner = np.zeros(stas_original.shape[0])
    #pvalues_adjusted_marginal = np.zeros(stas_original.shape[0])
    
    for i in range(stas_original.shape[0]):
        s=0
        for j in range(n_iter):
            s += (np.abs(resampled_stats[j,:]) >= np.abs(stas_original[i])).mean()
        pvalues_adjusted_reiner[i] = s/n_iter
        #pvalues_adjusted_marginal[i] = (np.abs(resampled_stats[:,i]) >= np.abs(stas_original[i])).mean()
        
    return pvalues_adjusted_reiner#, pvalues_adjusted_marginal

def estimate_m0_bky(pvalues, alpha=0.05, method='fdr_bh', is_sorted=False):
    # estimate the proportion of true null hypotheses
    # using the Benjamini-Krieger-Yekutieli method
    # from: "Adaptive linear step-up procedures that control the false discovery rate", (2006)
    # definition 6
    # returns a conservative estimate of the number of true null hypotheses and the adjusted alpha
    m = len(pvalues)
    q = alpha / (1 + alpha)
    r = multitest.multipletests(pvalues, alpha=q, method=method, is_sorted=is_sorted)[0].sum()
    m0 = m - r
    
    return m0, m/m0 * alpha

def estimate_m0_bh(pvalues, alpha=0.05, is_sorted=False):
    # estimate the proportion of true null hypotheses
    # using the Benjamini-Hochberg method
    # returns a conservative estimate of the number of true null hypotheses
    # and the adjusted alpha
    # This is from the same paper as the Benjamini-Krieger-Yekutieli method
    # definition 3
    # Step 1: Use the linear step-up procedure at level alpha
    # Step 2: Estimate m0(k) by (m+1-k)/(1-p(k))
    # Step 3: starting with k=2 stop when for the first time m0(k) > m0(k-1)
    # Step 4: estimate m0 = min(m0(k), m) rounding up to the nearest integer
    # Step 5: return m0 and alpha/m0
    if is_sorted:
        pvals = pvalues
    else:
        pvals = np.sort(pvalues)
    
    m = pvalues.shape[0]
    if multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=True)[0].sum() == 0:
        return m, alpha
    
    m0 = np.zeros(m)
    m0[0] = m/(1-pvals[0])
    for k in range(2,m+1):
        m0[k-1] = (m+1-k)/(1-pvals[k-1])
        if m0[k-1] > m0[k-2]:
            break
    m0 = np.min([m0[k-1], m])
    return m0, alpha*m/m0
    

#############################################

def true_positives(rejections, true_effects):
    # rejections is a boolean array of length n
    # true_effects is a boolean array of length n
    return np.sum(rejections & true_effects)

def false_positives(rejections, true_effects):
    # rejections is a boolean array of length n
    # true_effects is a boolean array of length n
    return np.sum(rejections & ~true_effects)

def true_negatives(rejections, true_effects):
    # rejections is a boolean array of length n
    # true_effects is a boolean array of length n
    return np.sum(~rejections & ~true_effects)

def false_negatives(rejections, true_effects):
    # rejections is a boolean array of length n
    # true_effects is a boolean array of length n
    return np.sum(~rejections & true_effects)

def power(rejections, true_effects):
    # rejections is a boolean array of length n
    # true_effects is a boolean array of length n
    return true_positives(rejections, true_effects) / np.sum(true_effects)


def print_FDR_and_POWER(fp,rejections, n_effects):
    FDP =[]
    for i in range(len(fp)):
        #denom = result['FP'][i]+result['TP'][i]
        denom = rejections[i] 
        if denom == 0:
            FDP.append(0)
            continue
        FDP.append(fp[i]/denom)

    tp = rejections - fp
    power = tp/n_effects
    
    print(np.mean(FDP), np.mean(power))


def print_ERROR_and_POWER(fp, rejections, n_effects):
    FDP =[]
    for i in range(len(fp)):
        denom = rejections[i] 
        if rejections[i]==0:
            FDP.append(0)
            continue
        FDP.append(fp[i]/rejections[i])

    tp = rejections - fp
    POWER = str(np.mean(tp/n_effects))
    FDR = str(np.mean(FDP))
    FWER = str(np.mean(fp>=1))

    print('FWER: ' + FWER + ', FDR: ' + FDR + ', POWER: ' + POWER)



#############################################





def calculate_relevant_measures(result):
    # calculate relevant measures for each method:
    # Sensitivity, Specificity, FDR, FWER, average power, average proportion of rejections to true effects.
    # additionally, for any key with 'M_eff' in it, calculate the average M_eff, and same for 'n_effects',
    # and anything with 'm0' in it, calculate the average m0 estimate (use same keys as in the result dictionary)
    # This is the structure of the result dictionary
    '''
    {'TAG':TAG,'WY':{'FN': FN_wy, 'FP': FP_wy, 'rejections': rejections_wy}, \
            'REINER_BH':{'FN': FN_reiner_bh, 'FP': FP_reiner_bh, 'rejections': rejections_reiner_bh},\
            'REINER_BY':{'FN': FN_reiner_by, 'FP': FP_reiner_by, 'rejections': rejections_reiner_by},\
            'n_effects': n_effects, 'n_tests': counts}
    '''
    measures = {}
    m1 = result['n_effects']
    m = result['n_tests']
    m0 = m - m1
    for key in result.keys():
        if key == 'TAG':
            # add tag to measures
            measures['TAG'] = result['TAG']
        # if the value is a dictionary, then we need to calculate the measures
        # We know there are always FN, FP, and rejections
        elif isinstance(result[key],dict):
            FN = result[key]['FN']
            V = FP = result[key]['FP']
            R = result[key]['rejections']
            S = TP = R - FP
            TN = m - R - FN
            # make a dictionary to store the measures with the same keys as the result dictionary
            measures[key] = {}
            measures[key]['Sensitivity'] = np.mean(TP/(TP+FN))  
            measures[key]['Specificity'] = np.mean(TN/(TN+FP))
            FDP = [v/r if r > 0 else 0 for v,r in zip(V,R)]
            measures[key]['FDR'] = np.mean(FDP)
            measures[key]['FDP_var'] = np.var(FDP)
            measures[key]['FWER'] = np.mean(V>=1)
            measures[key]['Average_Power'] = np.mean(S/m1)
            # accuracy is proportion of correctly assigned tests
            measures[key]['Accuracy'] = np.mean((TP+TN)/m)
            measures[key]['Average_Proportion_Rejections'] = np.mean(R/m1)

        # if the value is an array, then we need to calculate the average
        elif isinstance(result[key],np.ndarray):
            # if the key has 'M_eff' in it, then calculate the average M_eff
            if 'M_eff' in key:
                new_key = key + '_average_proportion'
                measures[new_key] = np.mean(result[key]/m)
            if 'm0' in key:
                new_key = key + '_average_proportion'
                measures[new_key] = np.mean(result[key]/m0)
                        
    return measures

def print_relevant_measures(measures):
    for key in measures.keys():
        if isinstance(measures[key],dict):
            print(key)
            for k in measures[key].keys():
                print(k,measures[key][k])
        else:
            print(key,measures[key])
            if key == 'TAG':
                print()
            
            

def save_relevant_measures(measures, filename):
    with open(filename, 'w') as f:
        for key in measures.keys():
            if isinstance(measures[key],dict):
                f.write(key + '\n')
                for k in measures[key].keys():
                    f.write(k + ' ' + str(measures[key][k]) + '\n')
            else:
                f.write(key + ' ' + str(measures[key]) + '\n')



def pickle_results(results, filename):
    with open(filename, 'wb') as f:
        pickle.dump(results, f)
        
def load_results(filename):
    with open(filename, 'rb') as f:
        results = pickle.load(f)
    return results
