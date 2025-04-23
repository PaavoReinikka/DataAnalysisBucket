import numpy as np

def reliable(x, tol, accept):
    '''
    checks gap status for a single test over one group
    by counting the number of acceptable status codes in the group
    and comparing the result against tol.
    '''
    s = 0
    for elem in x:
        s += elem in accept
    flag = s/len(x) >= tol
    n_unreliable = len(x) - s
    return flag, n_unreliable
            
    
def reliability(gapX, gapY, tol=0.8, accept=[0, 8, 64, 128]):
    '''
    checks gap status the entire data
    by counting the number of acceptable status codes in each group
    and comparing the result against tol 
    
    returns:
        flags: numpy array of bools, True if the test X-Y is reliable
        countsX: numpy array of ints, number of unreliable tests in X
        countsY: numpy array of ints, number of unreliable tests in Y
    '''
    assert gapX.shape == gapY.shape
    
    flags = np.zeros(gapX.shape[0], dtype=bool)
    countsX = np.zeros(gapX.shape[0], dtype=int)
    countsY = np.zeros(gapX.shape[0], dtype=int)
    
    for i in range(gapX.shape[0]):
        flagX, countX = reliable(gapX[i,:], tol, accept)
        flagY, countY = reliable(gapY[i,:], tol, accept)
        flags[i] = flagX | flagY
        countsX[i] = countX
        countsY[i] = countY
        
    return flags, countsX, countsY