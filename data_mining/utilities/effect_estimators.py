import numpy as np
import pandas as pd

def bootstrapped_effect(group1, group2, it=1000, replace=True):
    
    effect = lambda x1,x2: np.mean(x1) - np.mean(x2) 
    n1, n2 = len(group1), len(group2)
    effects = np.zeros((it,))
    
    for i in range(it):
        s1 = np.random.choice(group1, n1, replace=replace)
        s2 = np.random.choice(group2, n2, replace=replace)
        effects[i] = effect(s1, s2)
        
    return np.mean(effects), effects


def sampled_effect(group1, group2, it=1000):
    
    n1, n2 = len(group1), len(group2)
    effects = np.zeros((it,))
    
    for i in range(it):
        s1 = np.random.choice(group1, 1)
        s2 = np.random.choice(group2, 1)
        effects[i] = s1 - s2
        
    return np.mean(effects), effects
