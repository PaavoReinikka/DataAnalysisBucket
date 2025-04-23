import numpy as np

def get_ind(A, I, J):
    max_i, max_j = -1, -1
    max_val = -1e6
    for i in I:
        for j in J:
            val = A[i,j]
            if(val>max_val):
                max_val = val
                max_i, max_j = i, j
    return max_i, max_j

def VAT(D):
    n = D.shape[0]
    K = set(range(n))
    P = np.zeros((n,), dtype=int)
   
    #init
    i,_ = get_ind(D,K,K)
    P[0] = i
    I, J = set([i]), K - set([i])
    for r in range(1, n):
        i,j = get_ind(-D,I,J)
        P[r] = j
        I = I | set([j])
        J = J - set([j])
   
    D_new = np.zeros_like(D)
    for i in range(n):
        for j in range(n):
            D_new[i,j] = D[P[i],P[j]]
    return D_new
   
