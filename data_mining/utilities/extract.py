import numpy as np
import pandas as pd

def display_pvalues(MTs, Ps, Ts, Es, title=None):
    if(title is not None):
        print(title)
    for i, tup in enumerate(zip(Ps, Ts, Es)):
        print("-----------------------")
        print("Peak:", MTs[i])
        print("P-values:", tup[0])
        print("Tests:", tup[1])
        print("Effects:", tup[2])
        

def extract_pvalues(df, MTs, mt_col=' masstime', pvalue_col=" p_FDR", test_col=' test', effect_col=' FC', head=''):
    pvalues=[]
    tests=[]
    effects=[]
    for mt in MTs:
        p, t, e = extract_pvalues_(df, mt, mt_col, pvalue_col, test_col, effect_col, head)
        pvalues.append(p)
        tests.append(t)
        effects.append(e)
    return MTs, pvalues, tests, effects

def extract_pvalues_(df,MT, mt_col=' masstime', pvalue_col=" p_FDR", test_col=' test', effect_col=' FC', head=''):
    pvalues=[]
    tests=[]
    effects=[]
    inds=np.where(df[mt_col]==head + MT)[0]
    for ind in inds:
        pvalues.append(df[pvalue_col].iloc[ind])
        tests.append(df[test_col].iloc[ind])
        effects.append(df[effect_col].iloc[ind])
        
    return pvalues, tests, effects

def get_feature_matrix(masstime, data, significant=True):
    
    # from "data" dataframe, extract significant features
    # based on mass and rtime in "masstime" (series or array)
    # ignore first 3 rows, first 8 columns and return the 
    # float values in a numpy array
    
    if not significant:
        return data.iloc[3:,8:].to_numpy(dtype=float).T
    inds=[]
    for i in range(3,data.shape[0]):
        mt=data.iloc[i,4] + '@' + data.iloc[i,5]
        if mt in masstime:
            inds.append(i)
    inds=np.unique(np.asarray(inds, dtype=int))
    return data.iloc[inds,8:].to_numpy(dtype=float).T

def get_features(masstime, data, significant=True, skip=[0,1,2], head=''):

    # same as get_feature_matrix, except keep the same data format
    # as the input dataframe
    
    if not significant:
        return data
    inds=skip
    for i in range(len(skip),data.shape[0]):
        mt=head + data.iloc[i,4] + '@' + data.iloc[i,5]
        if mt in masstime:
            inds.append(i)
    inds=np.unique(np.asarray(inds, dtype=int))
    return data.iloc[inds,:]

