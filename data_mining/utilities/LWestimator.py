import numpy as np

def mod_frobenius_prod(A, B):
    '''
    computes the modified frobenius product between two matrices A and B
    as the sum of squared differences divided by the number of features
    '''
    assert A.shape == B.shape, "Matrices must have the same shape"
    diff = A - B
    squared_diff = diff ** 2
    sum_squared_diff = np.sum(squared_diff)
    num_features = A.shape[1]
    mod_frobenius = sum_squared_diff / num_features
    return mod_frobenius

def mod_frobenius_norm(A):
    '''
    computes the modified frobenius norm of matrix A
    as the square root of the sum of squared elements divided by the number of elements
    '''
    return np.sqrt(mod_frobenius_prod(A, A))

# NOTE: Check the implementation against sklearn.covariance.LedoitWolf for correctness
# This is very hastily put together and may contain errors/fibs...
def ledoit_wolf_shrinkage(X):
    '''
    applies the Ledoit-Wolf shrinkage to the given covariance matrix.
    
    Parameters:
        cov_matrix: numpy array, the sample covariance matrix
    Returns:
        shrunk_cov: numpy array, the shrunk covariance matrix
    '''
    cov_matrix = np.cov(X, rowvar=False)
    n = cov_matrix.shape[0]
    m = mod_frobenius_prod(cov_matrix, np.identity(n))
    d = mod_frobenius_prod(cov_matrix - (m * np.identity(n)), cov_matrix - (m * np.identity(n)))
    b2 = d**2
    b2hat = np.zeros_like(cov_matrix)
    for i in range(X.shape[0]):
        xi = X[i, :].reshape(-1, 1)
        outer_prod = np.dot(xi, xi.T)
        diff = outer_prod - cov_matrix
        b2hat += mod_frobenius_prod(diff, diff)
    b2hat /= X.shape[0]**2
    b2 = min(b2, b2hat)
    a2 = d**2 - b2
    
    cov_ledoit_wolf = (b2 * m * np.identity(n) + a2 * cov_matrix) / d**2
    return cov_ledoit_wolf
