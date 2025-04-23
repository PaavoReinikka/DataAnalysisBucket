import numpy as np
import scipy.stats as stats
import time
import multiprocessing
from joblib import Parallel, delayed


def MWSL(dataX, dataY, alpha = 0.05, K = None):
    """
    from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03975-2
    
    MWSL: Metabolome-Wide Significance Level
    Solves the multiple testing problem by permutation testing.
    Accounts for the correlation between the tests by permuting the labels of the samples.
    
    Idea:
    • Step (1): Shuffle i.e. re-sample without replacement, the outcome variable Y.
    In this way, the n subjects are re-sampled under the null hypothesis of no association.
    
    • Step (2): To estimate the relationship between the outcome and the set of features
    while accounting for possible confounding effects, run M t-tests.
    
    • Step (3): Extract the minimum of the set of M p values as this indicates the highest
    threshold value which would reject all M null hypotheses.
    
    • Step (4): Repeat Step (1)-(3) for K times, where K is at least n/2 times M. The K
    minimum p values are the elements of the new vector q.
    
    • Step (5): Sort the elements of q, and take the ( alpha*K)-value of this vector. This value
    is the MWSL estimate. An approximate confidence interval can be obtained by
    treating the true position of the MWSL estimate as a Binomial random variable
    with parameters K and alpha . Then, using the Normal approximation to the Binomial,
    we obtain the z_(1-alpha)% confidence limits by extracting the elements of q in positions
    (alpha*K) ± (1 - alpha)*sqrt(alpha*K(1 - alpha).
    
    • Step (6): Compute the effective number of tests (ENT) defined as the number of
    independent tests that would be required to obtain the same significance level using
    the Bonferroni correction ENT = alpha/MWSL . The ENT estimate measures the extent that
    the M markers are non-redundant. Therefore, the ratio R = ENT/M % of the effective
    and the actual number of tests (ANT or M) is a measure of the dependence among
    features, which is expected to be closer to 0% when highly correlated features are
    considered
    
    X: data matrix
    it: number of permutations
    alpha: significance level
    """
    assert alpha > 0 and alpha < 1, "alpha must be between 0 and 1"
    assert dataX.shape[0] == dataY.shape[0], "dataX and dataY must have the same number of rows"
    
    data = np.hstack((dataX, dataY))
    nX = dataX.shape[1]
    n=data.shape[1]
    M=data.shape[0]
    
    if K is None:
        K = int(np.ceil(n/2*M))
        print("K not specified, using K = n/2*M = ", K)
    else:
        assert K > n/2*M, "K must be greater than n/2*M"
    
    start = time.time()
    
    q=np.zeros(K)
    all_pvalues = np.zeros((K, M))
    for i in range(K):
        rand_ind = np.random.permutation(n)
        x, y = data[:,rand_ind[:nX]], data[:,rand_ind[nX:]]
        result = stats.ttest_ind(x, y, axis=1, equal_var=False, nan_policy='raise')
        q[i] = result.pvalue.min()
        all_pvalues[i,:] = result.pvalue
        
    q.sort()
    MWSL = q[int(alpha*K)]
    ENT = alpha/MWSL
    R = ENT/M
    MWSL_CI = (q[int(alpha*K - 1 - (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))], q[int(alpha*K - 1 + (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))])
    
    # make a dictionary with the results
    result = {}
    result['MWSL'] = MWSL
    result['ENT'] = ENT
    result['R'] = R
    result['MWSL_CI'] = MWSL_CI
    result['all_pvalues'] = all_pvalues
    
    end = time.time()
    print("Time elapsed: ", end-start)
    
    return result

def MWSL_parallel(dataX, dataY, alpha = 0.05, K = None):
    """
    from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03975-2
    
    MWSL: Metabolome-Wide Significance Level
    Solves the multiple testing problem by permutation testing.
    Accounts for the correlation between the tests by permuting the labels of the samples.
    
    Idea:
    • Step (1): Shuffle i.e. re-sample without replacement, the outcome variable Y.
    In this way, the n subjects are re-sampled under the null hypothesis of no association.
    
    • Step (2): To estimate the relationship between the outcome and the set of features
    while accounting for possible confounding effects, run M t-tests.
    
    • Step (3): Extract the minimum of the set of M p values as this indicates the highest
    threshold value which would reject all M null hypotheses.
    
    • Step (4): Repeat Step (1)-(3) for K times, where K is at least n/2 times M. The K
    minimum p values are the elements of the new vector q.
    
    • Step (5): Sort the elements of q, and take the ( alpha*K)-value of this vector. This value
    is the MWSL estimate. An approximate confidence interval can be obtained by
    treating the true position of the MWSL estimate as a Binomial random variable
    with parameters K and alpha . Then, using the Normal approximation to the Binomial,
    we obtain the z_(1-alpha)% confidence limits by extracting the elements of q in positions
    (alpha*K) ± (1 - alpha)*sqrt(alpha*K(1 - alpha).
    
    • Step (6): Compute the effective number of tests (ENT) defined as the number of
    independent tests that would be required to obtain the same significance level using
    the Bonferroni correction ENT = alpha/MWSL . The ENT estimate measures the extent that
    the M markers are non-redundant. Therefore, the ratio R = ENT/M % of the effective
    and the actual number of tests (ANT or M) is a measure of the dependence among
    features, which is expected to be closer to 0% when highly correlated features are
    considered
    
    X: data matrix
    K: number of permutations
    alpha: significance level
    """
    assert alpha > 0 and alpha < 1, "alpha must be between 0 and 1"
    assert dataX.shape[0] == dataY.shape[0], "dataX and dataY must have the same number of rows"
    
    data = np.hstack((dataX, dataY))
    nX = dataX.shape[1]
    n=data.shape[1]
    M=data.shape[0]
    
    if K is None:
        K = int(np.ceil(n/2*M))
        print("K not specified, using K = n/2*M = ", K)
    else:
        assert K > n/2*M, "K must be greater than n/2*M"
    
    q=np.zeros(K)
    all_pvalues = np.zeros((K, M))
    
    # parallel processing
    start = time.time()
    num_cores = multiprocessing.cpu_count()
    print("Running on {} cores".format(num_cores))
    
    # Generate the random indices
    rand_inds = np.stack([np.random.permutation(n) for i in range(K)], axis=0)
    results = Parallel(n_jobs=num_cores)(delayed(stats.ttest_ind)(data[:,rand_inds[i,:nX]],\
                                                                    data[:,rand_inds[i,nX:]],\
                                                                        axis=1, equal_var=False, nan_policy='raise') for i in range(K))
    
    
    q = np.asarray([np.min(result.pvalue) for result in results])
    all_pvalues = np.asarray([result.pvalue for result in results])
    
    q.sort()
    MWSL = q[int(alpha*K)]
    ENT = alpha/MWSL
    R = ENT/M
    MWSL_CI = (q[int(alpha*K - 1 - (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))], q[int(alpha*K - 1 + (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))])
    
    # make a dictionary with the results
    result = {}
    result['MWSL'] = MWSL
    result['ENT'] = ENT
    result['R'] = R
    result['MWSL_CI'] = MWSL_CI
    result['all_pvalues'] = all_pvalues
    
    end = time.time()
    print("Time elapsed: ", end-start)
    
    return result


def MWSL_independent(dataX, dataY, alpha=0.05, K=None):
    # Like the above functions but imposes independence between the tests
    # by permuting the labels of the samples independently for each test
    # instead of permuting the labels of the samples once for all tests
    
    assert alpha > 0 and alpha < 1, "alpha must be between 0 and 1"
    assert dataX.shape[0] == dataY.shape[0], "dataX and dataY must have the same number of rows"
    
    data = np.hstack((dataX, dataY))
    nX = dataX.shape[1]
    n=data.shape[1]
    M=data.shape[0]
    
    if K is None:
        K = int(np.ceil(n/2*M))
        print("K not specified, using K = n/2*M = ", K)
    else:
        assert K > n/2*M, "K must be greater than n/2*M"
    
    start = time.time()
    
    q=np.zeros(K)
    all_pvalues = np.zeros((K, M))
    
    for i in range(K):
        rand_ind = np.stack([np.random.permutation(n) for i in range(M)], axis=0)
        
        x, y = np.ones((M, nX)), np.ones((M, nX))
        for j in range(M):
            x[j,:] = data[j,rand_ind[j,:nX]]
            y[j,:] = data[j,rand_ind[j,nX:]]
        
        
        all_pvalues[i,:] = stats.ttest_ind(x, y, axis=1, equal_var=False, nan_policy='raise').pvalue
        q[i] = all_pvalues[i,:].min()
        
    q.sort()
    MWSL = q[int(alpha*K)]
    ENT = alpha/MWSL
    R = ENT/M
    MWSL_CI = (q[int(alpha*K - 1 - (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))], q[int(alpha*K - 1 + (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))])
    
    # make a dictionary with the results
    result = {}
    result['MWSL'] = MWSL
    result['ENT'] = ENT
    result['R'] = R
    result['MWSL_CI'] = MWSL_CI
    result['all_pvalues'] = all_pvalues
    
    end = time.time()
    print("Time elapsed: ", end-start)
    
    return result

def MWSL_independent_parallel(dataX, dataY, alpha=0.05, K=None):
    # Like the above function but uses parallel processing
    
    assert alpha > 0 and alpha < 1, "alpha must be between 0 and 1"
    assert dataX.shape[0] == dataY.shape[0], "dataX and dataY must have the same number of rows"
    
    data = np.hstack((dataX, dataY))
    nX = dataX.shape[1]
    n=data.shape[1]
    M=data.shape[0]
    
    if K is None:
        K = int(np.ceil(n/2*M))
        print("K not specified, using K = n/2*M = ", K)
    else:
        assert K > n/2*M, "K must be greater than n/2*M"
        
    start = time.time()
    
    q=np.zeros(K)
    all_pvalues = np.zeros((K, M))
    
    # parallel processing
    num_cores = multiprocessing.cpu_count()
    print("Running on {} cores".format(num_cores))
    # Generate the random indices of size K x M x n
    rand_inds = np.stack([np.stack([np.random.permutation(n) for i in range(M)], axis=0) for j in range(K)], axis=0)
    
    @staticmethod
    def ttest_ind_parallel(data, r_inds):
        x, y = np.ones((M, nX)), np.ones((M, nX))
        for j in range(M):
            x[j,:] = data[j,r_inds[j,:nX]]
            y[j,:] = data[j,r_inds[j,nX:]]
        return stats.ttest_ind(x, y, axis=1, equal_var=False, nan_policy='raise')

    #data_expanded = data.reshape((1, M, n))
    #data_expanded = np.tile(data_expanded, (K, 1, 1))
    #results = Parallel(n_jobs=num_cores)(delayed(ttest_ind_parallel)(data_expanded[i], rand_inds[i]) for i in range(K))
    results = Parallel(n_jobs=num_cores)(delayed(ttest_ind_parallel)(data, rand_inds[i]) for i in range(K))
    
    
    all_pvalues = np.asarray([result.pvalue for result in results]).reshape((K, M))
    q = np.asarray([np.min(result.pvalue) for result in results])   
    
    q.sort()
    MWSL = q[int(alpha*K)]
    ENT = alpha/MWSL
    R = ENT/M
    MWSL_CI = (q[int(alpha*K - 1 - (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))], q[int(alpha*K - 1 + (1 - alpha)*np.sqrt(alpha*K*(1 - alpha)))])
    
    # make a dictionary with the results
    result = {}
    result['MWSL'] = MWSL
    result['ENT'] = ENT
    result['R'] = R
    result['MWSL_CI'] = MWSL_CI
    result['all_pvalues'] = all_pvalues
    
    end = time.time()
    print("Time elapsed: ", end-start)
    
    return result

