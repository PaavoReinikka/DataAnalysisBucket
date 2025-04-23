import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

def n_eff_spectral(C, alpha = 0.05, method='1', factor = 1):
    """
    Calculate the effective number of tests from a spectral decomposition
    of the correlation matrix. The function also returns Sidak and Bonferroni 
    corrected significance level.

    methods:
    1. https://www.nature.com/articles/6889010
    2. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1181954/
    3. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03975-2
    
    Parameters:
    C: numpy array of shape (features, features)

    """
    # Get the eigenvalues of the correlation matrix
    eigvals = np.linalg.eigvals(C)
    M = len(eigvals)
    # Calculate the effective number of tests
    if method == '1':
        M_eff = M*(1 - (M - 1)*np.var(eigvals) / M**2)
    elif method == '2':
        M_eff = 1 + (M - 1)*(1 - np.var(eigvals) / M)
    elif method == '3':
        # assert that the eigenvalues are positive and real
        assert np.all(eigvals >= 0) and np.all(np.isreal(eigvals)), 'Eigenvalues must be positive and real'
        num = (np.sqrt(eigvals).sum()/np.log(eigvals.max()))**2
        den = eigvals.sum()/eigvals.max() + np.sqrt(eigvals.max())
        M_eff = num/den
    else:
        raise ValueError('Invalid method')
    
    M_eff = int(M_eff*factor)
    result = dict()
    result['M_eff'] = M_eff
    result['Sidak'] = 1 - (1 - alpha)**(1/M_eff)
    result['Bonferroni'] = alpha/M_eff

    return result

def n_eff_cluster(COR, threshold = 0.1, alpha = 0.05, method = 'complete', plot = False, return_labels = False, factor = 1):
    """
    Calculate the effective number of tests by performing hierarchical clustering.
    The function also returns Sidak and Bonferroni corrected significance level.

    Parameters:
    X: numpy array of features (obs, features) if metric = 'correlation'
    or otherwise a distance/kernel matrix (features, features)
    """

    CD = 1 - COR
    CD -= np.diag(np.diag(CD))
    CD=np.clip(CD,0,2)
    # Calculate the distance matrix
    D = squareform(CD,force='tovector')
    # Perform hierarchical clustering
    Z = hierarchy.linkage(D, method='complete')
    # Get the clusters
    clusters = hierarchy.fcluster(Z, t=threshold, criterion='distance')
    # Get the number of clusters
    M_eff = len(np.unique(clusters))
    M_eff = int(M_eff*factor)
    
    if plot:
        hierarchy.dendrogram(Z, truncate_mode='none', labels=np.asarray(np.arange(CD.shape[0])))
        plt.axhline(y=threshold, c='grey', lw=1, linestyle='dashed')
        plt.yticks(np.arange(0,2.1,0.2))
        plt.title("{} clusters (out of {}) with threshold of {}".format(np.unique(clusters).shape[0], CD.shape[0], threshold))
        plt.show()
    
    
    result = dict()
    result['M_eff'] = M_eff
    result['Sidak'] = 1 - (1 - alpha)**(1/M_eff)
    result['Bonferroni'] = alpha/M_eff
    
    if return_labels:
        result['labels'] = clusters

    return result
