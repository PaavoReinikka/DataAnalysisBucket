import numpy as np
from numba import njit

@njit(cache=True)
def np_all_axis0(x):
    """Numba compatible version of np.all(x, axis=0)."""
    out = np.ones(x.shape[1], dtype=np.bool8)
    for i in range(x.shape[0]):
        out = np.logical_and(out, x[i, :])
    return out

@njit(cache=True)
def np_all_axis1(x):
    """Numba compatible version of np.all(x, axis=1)."""
    out = np.ones(x.shape[0], dtype=np.bool8)
    for i in range(x.shape[1]):
        out = np.logical_and(out, x[:, i])
    return out

@njit(cache=True)
def np_any_axis0(x):
    """Numba compatible version of np.any(x, axis=0)."""
    out = np.zeros(x.shape[1], dtype=np.bool8)
    for i in range(x.shape[0]):
        out = np.logical_or(out, x[i, :])
    return out

@njit(cache=True)
def np_any_axis1(x):
    """Numba compatible version of np.any(x, axis=1)."""
    out = np.zeros(x.shape[0], dtype=np.bool8)
    for i in range(x.shape[1]):
        out = np.logical_or(out, x[:, i])
    return out

# TODO: add a function to replace un-used centroids, e.g., by random re-initialization

@njit
def kmeans(X, k, n_iter, init_centroids):
    #Fast parallel kmeans
    #Original implementation from old numba examples (here slightly modified)
    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids

    for l in range(n_iter):
        dist = np.array([[np.sqrt(np.sum((X[i, :] - centroids[j, :])**2))
                          for j in range(k)] for i in range(N)])

        predictions = np.array([dist[i, :].argmin() for i in range(N)])

        centroids = np.array([[np.sum(X[predictions == i, j])/np.sum(predictions == i)
                               for j in range(D)] for i in range(k)])

    return centroids, dist, predictions


@njit
def kmeans_spherical_v1(X, k, n_iter, init_centroids):

    # This implementation collapses empty clusters to origin.
    # The origin tends to pull other clusters towards it.
    # Not good for high values of k.
    
    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids

    for i in range(centroids.shape[0]):
        centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])

    for l in range(n_iter):
        dist = 1 - X@centroids.T

        predictions = np.array([dist[i, :].argmin() for i in range(N)])
        
        centroids = np.array([[np.sum(X[predictions == i, j])/max(1,np.sum(predictions == i))
                               for j in range(D)] for i in range(k)])

        for i in range(centroids.shape[0]):
            centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])

    return centroids, dist, predictions

@njit
def kmeans_spherical_v2(X, k, n_iter, init_centroids):
    
    # This implementation drops empty clusters

    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids

    for i in range(centroids.shape[0]):
        centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])
        
    for l in range(n_iter):
        dist = 1 - X@centroids.T

        predictions = np.array([dist[i, :].argmin() for i in range(N)])

        centroids = np.array([[np.sum(X[predictions == i, j])/np.sum(predictions == i)
                               for j in range(D)] for i in np.unique(predictions)])
        
        for i in range(centroids.shape[0]):
            centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])
            
    return centroids, dist, predictions

@njit
def kmeans_spherical_v3(X, k, n_iter, init_centroids):
    
    # This keeps empty clusters and does not pull other clusters towards the origin.
    # It is the same as v2 but with a different way of handling empty clusters.
    # You can also change the way empty clusters are handled -- e.g., by random re-initialization.
    
    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids

    for i in range(centroids.shape[0]):
        centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])

    for l in range(n_iter):
        dist = 1 - X@centroids.T

        predictions = np.array([dist[i, :].argmin() for i in range(N)])

        centroids_new = np.array([[np.sum(X[predictions == i, j])/max(1,np.sum(predictions == i))
                               for j in range(D)] for i in range(k)])
        
        mask = np_all_axis1(centroids_new == 0)
        centroids_new[mask,:] = centroids[mask,:]
        
        centroids = centroids_new
        
        for i in range(centroids.shape[0]):
            centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])

    return centroids, dist, predictions

@njit
def kmeans_spherical_v4(X, k, n_iter, init_centroids, discard_freq=10):

    # This keeps empty clusters and does not pull other clusters towards the origin.
    # It is the same as v2 but with a different way of handling empty clusters.
    # You can also change the way empty clusters are handled -- e.g., by random re-initialization.

    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids

    for i in range(centroids.shape[0]):
        centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])

    discarded_centroids = np.zeros((k,1))
    for l in range(n_iter):
        dist = 1 - X@centroids.T

        predictions = np.array([dist[i, :].argmin() for i in range(N)])

        centroids_new = np.array([[np.sum(X[predictions == i, j])/max(1,np.sum(predictions == i))
                                   for j in range(D)] for i in range(k)])

        # TODO: add a function to replace un-used centroids, e.g., by random re-initialization
        mask = np_all_axis1(centroids_new == 0)
        centroids_new[mask,:] = centroids[mask,:]
        discarded_centroids[mask] += 1
        
        if np_any_axis0(discarded_centroids >= discard_freq):
            centroids_new[discarded_centroids >= discard_freq,:] = np.random.randn(np.sum(discarded_centroids >= discard_freq), D)
            discarded_centroids[discarded_centroids >= discard_freq] = 0

        centroids = centroids_new

        for i in range(centroids.shape[0]):
            centroids[i,:] = centroids[i,:]/np.sqrt(centroids[i,:]@centroids[i,:])

    return centroids, dist, predictions


