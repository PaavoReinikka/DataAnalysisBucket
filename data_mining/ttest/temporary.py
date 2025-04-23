import numpy as np
from numpy.matlib import repmat
import pandas as pd

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

def compile_ttest_features(path: str, fname: str):
	df = pd.read_csv(path + fname,sep=';', header=None)

	# Compile the data for testing
	AREA = df.iloc[3:,8:48].to_numpy(dtype=float)
	GROUPS = df.iloc[0,8:48].to_numpy(dtype=str)
	GAPS = df.iloc[3:,88:].to_numpy(dtype=int)
	FEATURE_INFO = df.iloc[3:,:8].to_numpy(dtype=str)

	# return a dictionary
	return dict(AREA=AREA, GROUPS=GROUPS, GAPS=GAPS, FEATURE_INFO=FEATURE_INFO)

# this function takes a dictionary and a function that filters the features.
# The filtering function should return boolean values for each feature.
def filter_features(D: dict, filter: callable, verbose=False, **kwargs):
    if verbose:
        print('Filtering features')
        print('Features before:', D['AREA'].shape[0])
    # Filter out features that are not significant
    AREA = D['AREA']
    GAPS = D['GAPS']
    FEATURE_INFO = D['FEATURE_INFO']
    
    keep, filtering_measure = filter(AREA)
    filtering_measure = filtering_measure[keep][:,np.newaxis]
    AREA = AREA[keep,:]
    GAPS = GAPS[keep,:]
    FEATURE_INFO = FEATURE_INFO[keep,:]
    
    if verbose:
        print('Features removed:', np.sum(~keep))
        print('Features after:', np.sum(keep))
        if 'method' in kwargs:
            print('Filtering method:', kwargs['method'])
            
    # return a dictionary
    return dict(AREA=AREA, GAPS=GAPS, FEATURE_INFO=FEATURE_INFO, filtering_measure=filtering_measure)

def compile_ttest_groups(D, verbose=False):
    # Compile the data for testing
    patterns = ["aSYN--UT", "comb.--UT", "INFg--UT", "aSYN--INFg","comb.--INFg.", "aSYN--comb."]
    n_patterns = len(patterns)
    
    AREA = D['AREA']
    GAPS = D['GAPS']
    FEATURE_INFO = D['FEATURE_INFO']
    filtering_measure = D['filtering_measure']
    
    # Compile the test identifiers, that map test/row to test identifier string
    I = np.ones((AREA.shape[0],1), dtype=int)
    pattern_inds = np.vstack((0*I, 1*I, 2*I, 3*I, 4*I, 5*I)).flatten()
    FULL_FEATURE_INFO = repmat(FEATURE_INFO, n_patterns, 1)
    FULL_FILTERING_MEASURE = repmat(filtering_measure, n_patterns, 1)
    
    # Compile the data for testing
    aSYN = AREA[:,0:10]
    comb = AREA[:,10:20]
    IFNg = AREA[:,20:30]
    UT = AREA[:,30:40]

    dataX = np.vstack((aSYN, comb, IFNg, aSYN, comb, aSYN))
    dataY = np.vstack((UT,   UT,   UT,   IFNg, IFNg, comb))
    
    GAP_aSYN = GAPS[:,0:10]
    GAP_comb = GAPS[:,10:20]
    GAP_IFNg = GAPS[:,20:30]
    GAP_UT = GAPS[:,30:40]

    gapX = np.vstack((GAP_aSYN, GAP_comb, GAP_IFNg, GAP_aSYN, GAP_comb, GAP_aSYN))
    gapY = np.vstack((GAP_UT,   GAP_UT,   GAP_UT,   GAP_IFNg, GAP_IFNg, GAP_comb))
    
    if verbose:
        print('Number of tests:', dataX.shape[0])
    
    return dict(dataX=dataX, dataY=dataY, \
                gapX=gapX, gapY=gapY, \
                pattern_inds=pattern_inds, \
                FULL_FEATURE_INFO=FULL_FEATURE_INFO, \
                FULL_FILTERING_MEASURE=FULL_FILTERING_MEASURE)
    
def remove_unreliable_tests(D, tol=0.8, accept=[0, 8, 64, 128], verbose=False):
    if verbose:
        print('Removing unreliable tests')
        print('Tolerance:', tol)
        print('Acceptable status codes:', accept)
        print('Tests before:', D['dataX'].shape[0])
    
    dataX = D['dataX']
    dataY = D['dataY']
    gapX = D['gapX']
    gapY = D['gapY']
    pattern_inds = D['pattern_inds']
    FULL_FEATURE_INFO = D['FULL_FEATURE_INFO']
    FULL_FILTERING_MEASURE = D['FULL_FILTERING_MEASURE']
    
    keep, countX, countY = reliability(gapX, gapY, tol, accept)
    if verbose:
        print('Tests removed:', np.sum(~keep))
        print('Tests after:', np.sum(keep))
    
    dataX = dataX[keep,:]
    dataY = dataY[keep,:]
    # gap data is not used anymore, only the counts
    FULL_FEATURE_INFO = FULL_FEATURE_INFO[keep,:]
    FULL_FILTERING_MEASURE = FULL_FILTERING_MEASURE[keep]
    pattern_inds = pattern_inds[keep]
    countX = countX[keep]
    countY = countY[keep]
    
    return dict(dataX=dataX, dataY=dataY, \
                countX=countX, countY=countY, \
                pattern_inds=pattern_inds, \
                FULL_FEATURE_INFO=FULL_FEATURE_INFO, \
                FULL_FILTERING_MEASURE=FULL_FILTERING_MEASURE)
    
    