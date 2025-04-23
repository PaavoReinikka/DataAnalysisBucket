import numpy as np
from scipy import stats
from scipy.spatial.distance import squareform
from statsmodels.stats import multitest
from joblib import Parallel, delayed
import multiprocessing
import time

# test the significance of each correlation coefficient
def test_correlation(X, alpha=0.05, keep_dims=False):
    # This function tests the significance of each correlation coefficient between the features in X
    # X is a matrix of size n_features x n_samples
    
    pvalues = np.zeros((X.shape[0], X.shape[0]))

    # time the loop
    start = time.time()
    for i in range(X.shape[0]):
        for j in range(i+1, X.shape[0]):
            pvalues[i,j] = stats.pearsonr(X[i,:], X[j,:])[1]
    end = time.time()
    print("Time elapsed: ", end-start)
    pvalues_vector = pvalues[np.triu_indices(pvalues.shape[0], k=1)]

    rejected, p_FDR, alphacSidak, alphacBonf = multitest.multipletests(pvalues_vector, alpha=alpha, method='fdr_by')
    
    # make a dict of the results
    result = {}
    if keep_dims:
        result['pvalues'] = squareform(pvalues_vector, force='tomatrix')
        result['rejected'] = squareform(rejected, force='tomatrix')
        result['p_FDR'] = squareform(p_FDR, force='tomatrix')
        result['corrcoef'] = np.corrcoef(X)
    else:
        result['pvalues'] = pvalues_vector
        result['rejected'] = rejected
        result['p_FDR'] = p_FDR
        result['corrcoef'] = np.corrcoef(X)[np.triu_indices(X.shape[0], k=1)]
            
    result['alphacSidak'] = alphacSidak
    result['alphacBonf'] = alphacBonf
    
    
    return result

def test_correlation_parallel(X, alpha=0.05, keep_dims=False):
    # Like test_correlation but uses parallel processing
    # This function tests the significance of each correlation coefficient between the features in X
    # X is a matrix of size n_features x n_samples
    
    # parallel processing
    start = time.time()
    num_cores = multiprocessing.cpu_count()
    print("Running on {} cores".format(num_cores))
    results = Parallel(n_jobs=num_cores)(delayed(stats.pearsonr)(X[i,:], X[j,:]) for i in range(X.shape[0]) for j in range(i+1, X.shape[0]))
    pvalues_vector = np.asarray(results)[:,1]
    print()
    end = time.time()
    print("Time elapsed: ", end-start)
    
    rejected, p_FDR, alphacSidak, alphacBonf = multitest.multipletests(pvalues_vector, alpha=alpha, method='fdr_by')
    
    # make a dict of the results
    result = {}
    if keep_dims:
        result['pvalues'] = squareform(pvalues_vector, force='tomatrix')
        result['rejected'] = squareform(rejected, force='tomatrix')
        result['p_FDR'] = squareform(p_FDR, force='tomatrix')
        result['corrcoef'] = np.corrcoef(X)
    else:
        result['pvalues'] = pvalues_vector
        result['rejected'] = rejected
        result['p_FDR'] = p_FDR
        result['corrcoef'] = np.corrcoef(X)[np.triu_indices(X.shape[0], k=1)]
    result['alphacSidak'] = alphacSidak
    result['alphacBonf'] = alphacBonf
    
    return result


def confidence_interval_correlation(COR, alpha=0.05, threshold=0.95):
    # This function calculates confidence intervals for the correlation coefficients
    # using the Fisher transformation
    # https://en.wikipedia.org/wiki/Fisher_transformation
    # https://www.statisticshowto.com/probability-and-statistics/correlation-coefficient-formula/
    # https://www.statisticshowto.com/probability-and-statistics/z-score/

    cor_vector = COR[np.triu_indices(COR.shape[0], k=1)]
    clipped_cor_vector = np.clip(cor_vector, np.finfo(float).eps-1, 1-np.finfo(float).eps)

    r_transformed = np.arctanh(clipped_cor_vector)
    SE = 1/np.sqrt(COR.shape[0]-3)
    # calculate the z-score for the 95% confidence interval, with multiple testing correction
    z = stats.norm.ppf(1-alpha/(2*len(cor_vector)))
    # calculate the confidence intervals
    CI_upper = np.tanh(r_transformed + 1* z * SE)
    CI_lower = np.tanh(r_transformed - 1* z * SE)

    # Check if the confidence intervals do not overlap with +/- 0.95
    CI_mask_vector = np.logical_or(CI_lower > threshold, CI_upper < -threshold)
    # convert the vector into a matrix
    CI_mask = np.zeros((COR.shape[0], COR.shape[0]), dtype=bool)
    CI_mask[np.triu_indices(COR.shape[0], k=1)] = CI_mask_vector
    CI_mask += CI_mask.T
    
    # make dictionary with the results
    result = {}
    result['COR'] = COR
    result['CI_mask'] = CI_mask
    result['cor_vector'] = cor_vector
    result['CI_mask_vector'] = CI_mask_vector
    result['CI_lower'] = CI_lower
    result['CI_upper'] = CI_upper
    
    return result
    
