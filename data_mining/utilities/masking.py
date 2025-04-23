import numpy as np
import pandas as pd


def make_proximity_mask(arr, tol, thresholding='hard', multiplier=None):
    
    arr = arr.reshape(len(arr),1)
    D = np.abs(arr - arr.T)
    mask = D <= tol
    if thresholding == 'soft':
        if(multiplier is None):
            multiplier = 1
        
        relu_ = np.vectorize(np.maximum)
        relu = lambda x: relu_(x,0)
        
        mask -= multiplier * np.clip(D,0,tol)
        mask = relu(mask)
        
    return mask 