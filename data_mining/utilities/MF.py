import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.preprocessing import StandardScaler

def SPCA(X, n_components, n_iter=100, alpha=0, l1_ratio=0.5):
    """
    Alternating Sparse PCA

    Parameters
    ----------
    X : array_like, shape (n_samples, n_features)
        Training data.
    alpha : float, optional (default=0)
        Regularization parameter.

    """
    #loss = lambda X, A, B: np.linalg.norm(X@A - X@B, ord='fro')**2 + alpha*((1-l1_ratio)*np.linalg.norm(B, ord=2) + l1_ratio*np.linalg.norm(A, ord=2)
    N, D = X.shape
    C = np.cov(X.T)
    eigvals, eigvecs = np.linalg.eig(C)
    ind = np.argsort(eigvals)[::-1]
    eigvals = eigvals[ind][:n_components]
    eigvecs = eigvecs[:, ind][:, :n_components]
    
    A = eigvecs
    B = np.zeros((D, n_components))
    x_star = np.sqrt(N*C)
    
    for i in range(n_iter):
        y_star = x_star @ A
        for j in range(n_components):
            B[:, j] = ElasticNet(alpha=alpha, l1_ratio=l1_ratio).fit(x_star, y_star[:,j]).coef_
        U, S, Vt = np.linalg.svd(X.T@X@B, full_matrices=True)
        A = U@Vt
        
    return A, B




def random_initialization(D,rank):
    number_of_documents = D.shape[0]
    number_of_terms = D.shape[1]
    U = np.random.uniform(1,2,(number_of_documents,rank))
    V = np.random.uniform(1,2,(rank,number_of_terms))
    return U,V.T

# NNDSVD method, Boutsidis & Gallopoulos, 2007
def nnd_svd_initialization(D,rank):
    u,s,v=np.linalg.svd(D,full_matrices=False)
    v=v.T
    w=np.zeros((D.shape[0],rank))
    h=np.zeros((rank,D.shape[1]))

    w[:,0]=np.sqrt(s[0])*np.abs(u[:,0])
    h[0,:]=np.sqrt(s[0])*np.abs(v[:,0].T)

    for i in range(1,rank):
        
        ui=u[:,i]
        vi=v[:,i]
        ui_pos=(ui>=0)*ui
        ui_neg=(ui<0)*-ui
        vi_pos=(vi>=0)*vi
        vi_neg=(vi<0)*-vi
        
        ui_pos_norm=np.linalg.norm(ui_pos,2)
        ui_neg_norm=np.linalg.norm(ui_neg,2)
        vi_pos_norm=np.linalg.norm(vi_pos,2)
        vi_neg_norm=np.linalg.norm(vi_neg,2)
        
        norm_pos=ui_pos_norm*vi_pos_norm
        norm_neg=ui_neg_norm*vi_neg_norm
        
        if norm_pos>=norm_neg:
            w[:,i]=np.sqrt(s[i]*norm_pos)/ui_pos_norm*ui_pos
            h[i,:]=np.sqrt(s[i]*norm_pos)/vi_pos_norm*vi_pos.T
        else:
            w[:,i]=np.sqrt(s[i]*norm_neg)/ui_neg_norm*ui_neg
            h[i,:]=np.sqrt(s[i]*norm_neg)/vi_neg_norm*vi_neg.T

    return w,h.T


# This does not enforce non-negativity (just a direct minimization of reconstruction
# error, possibly with L2 regularization term
# Needs to be cythonized
def matrix_factorization_slow(ratings, U_init, V_init, latent_dims, n_iter=5000, alpha=0.0002, beta=0.02, tol=1e-3):
    '''
    D: rating matrix
    U: #users x latent_dims
    V: #items x latent_dims
    alpha: learning rate
    beta: regularization parameter'''
    D = ratings
    U = U_init
    V = V_init.T
    n_u, n_i = D.shape
    
    for step in range(n_iter):
        for i in range(n_u):
            for j in range(n_i):
                if D[i,j] > 0:
                    # calculate error for prediction
                    e_ij = D[i,j] - np.dot(U[i,:],V[:,j])

                    for k in range(latent_dims):
                        # update with gradient information
                        U[i,k] = U[i,k] + alpha * (2 * e_ij * V[k,j] - beta * U[i,k])
                        V[k,j] = V[k,j] + alpha * (2 * e_ij * U[i,k] - beta * V[k,j])

        e = (D - U@V)[np.where(D>0)].reshape(-1)
        if(np.linalg.norm(e)<tol):
            break
    print("Number of steps: {}".format(step+1))
    return U, V.T


#Using multiplicative updates. Enforces non-negativity (e.g., Aggarwals book)
#This is a very "banana"-model - adding regularization should improve slightly.
#Mostly interested in trying this with hybrid recommender systems,  where
# initial sparse data is "densified" with content based method and only then
# using NMF/MF for final recommendations (or svd-low-rank, collaborative filtering, etc.)
def nn_matrix_factorization(X, U_init, V_init, n_iter=100, tol=1e-5):
    D = X
    U = U_init.copy()
    V = V_init.copy()
    n_u, n_i = D.shape
    K = V.shape[1]
    err = np.ones(n_iter)*np.inf
    for step in range(n_iter):
        #update V
        DtU = D.T@U
        VUtU = V@U.T@U
        V *= DtU / VUtU
        #update U
        DV = D@V
        UVtV = U@V.T@V
        U *= DV / UVtV
        # compute error
        #err[step] = np.linalg.norm(D-U@V.T, 'fro')
        err[step] = np.sqrt(((D-U@V.T)**2).sum())
        if err[step-1] - err[step] < tol:
            err = err[:step]
            print("Stop at step {}".format(step))
            break
    return U,V, err




def matrix_factorization(ratings, U_init, V_init, n_iter=100, n_iter_inner=10, alpha=0.0002, beta=0.02, tol=1e-5):
    '''
    D: rating matrix
    U: #users x latent_dims
    V: #items x latent_dims
    alpha: learning rate
    beta: regularization parameter'''
    D = ratings
    U = U_init.copy()
    V = V_init.copy().T
    n_u, n_i = D.shape
    err = np.ones(n_iter)*np.inf
    
    for step in range(n_iter):
        
        # This is a modified version of alternating least squares
        # where we do multiple steps of gradient descent for each
        # update of U and V
        
        # update U
        for inner_step in range(n_iter_inner):
            e_i = D - np.dot(U,V)
            U = U + alpha * (2 * e_i @ V.T - beta * U)
            #for i in range(n_u):
            #    e_i = D[i,:] - np.dot(U[i,:],V)
            #    U[i,:] = U[i,:] + alpha * (2 * e_i[np.newaxis,:] @ V.T - beta * U[i,:])
        
        # update V
        for inner_step in range(n_iter_inner):
            e_j = D - np.dot(U,V)
            V = V + alpha * (2 * U.T @ e_j - beta * V)
            #for j in range(n_i):
            #    e_j = D[:,j] - np.dot(U,V[:,j])
            #    V[:,j] = V[:,j] + alpha * (2 * e_j @ U - beta * V[:,j])
        
        # compute error
        err[step] = np.sqrt(((D-U@V)**2).sum())
        if np.abs(err[step-1] - err[step]) < tol:
            err = err[:step]
            print("Stop at step {}".format(step))
            break
        
    return U, V.T, err

    
    
