import numpy as np
import pandas as pd
from scipy.cluster.vq import vq
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import distance_matrix
from scipy.stats import mode
import matplotlib.pyplot as plt

def matching_dist(X1, X2, aggregator = np.mean):
    n1, n2 = X1.shape[0], X2.shape[0]
    result = np.empty((n1, n2))
    for i in range(n1):
        for j in range(n2):
            result[i,j]=aggregator(X1[i,:]!=X2[j,:])
            
    return result

def euclidean_dist(X1, X2):
    #TODO-alternative: implement pairwise distance function(s) wich return
    # a matrix of shape (n, k) for inputs with following sizes: 
    #  X1.shape=(n,d), X2.shape=(k,d)
    return distance_matrix(X1, X2, p=2)
    
    
def manhattan_dist(X1, X2):
    #TODO-alternative: implement pairwise distance function(s) wich return
    # a matrix of shape (n, k) for inputs with following sizes: 
    #  X1.shape=(n,d), X2.shape=(k,d)
    return distance_matrix(X1, X2, p=1)
    

def assign(D, k):#metric):
    '''Takes the distances to k current cluster representatives (last ,
       and returns the cluster assignments and the distance to
       the closest representative (for each datapoint)
    '''
    ind = np.argmin(D, axis=1)
    bool_ind = (np.arange(k).reshape(-1,1)==ind).T
    
    return ind, D[bool_ind]

def update_representatives(X, assignments, representatives, method):
    '''Based on provided assignments, iterates through every cluster
       and calculates the new representatives using method (e.g., np.mean/median or scipy.stats.mode).
    '''
    for i in np.unique(assignments):
        subX = X[np.where(assignments==i),:][0]
        representatives[i] = method(subX, axis=0)
    
    return representatives


def KRepresentatives(X, k, x_init, n_iters = 100, dist_func = euclidean_dist, rep_func=np.mean):
    
    representatives = x_init
    
    #scores_arr = np.empty((n_iters,1))
    #representatives_arr = np.empty((k,n_iters, X.shape[1]))
    
    for i in range(n_iters):
        D = dist_func(X, representatives)
        assignments, dist = assign(D,k)
        representatives = update_representatives(X, assignments, representatives, rep_func)
        
        #these are optionally used to get the whole optimization path
        #score_arr[i] = score(dist)
        #representatives_arr[:,i,:] = representatives
        
    return assignments#, representatives_arr, scores_arr
    
    

def my_mode(x, axis):
    #scipy's mode returns a tupple, and therefore needs a wrapper
    return mode(x, axis=axis)[0]

