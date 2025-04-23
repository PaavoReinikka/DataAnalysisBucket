import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


def print_duplicates(path, filenames, column_str, tol, split_token = '_'):
    dm, dt = tol[0], tol[1]
    token = split_token
    n_duplicates=0
    for i in range(1,9):
        for j in range(i+1,10):
            data1 = pd.read_csv(path + filenames[i], sep=';')
            data2 = pd.read_csv(path + filenames[j], sep=';')
            MT1 = np.unique(data1[column_str])
            MT2 = np.unique(data2[column_str])
            n_duplicates+=print_equals(MT1,MT2,filenames[i].split(token)[0],filenames[j].split(token)[0], dm, dt)
    print("Number of pairs: {}".format(n_duplicates))
    

def save_duplicates(path, filenames, column_str, tol, split_token = '_'):
    dm, dt = tol[0], tol[1]
    token = split_token
    n_duplicates=0
    fname_results = 'postTestDuplicates.txt'
    if os.path.exists(fname_results):
        os.remove(fname_results)
    for i in range(1,9):
        for j in range(i+1,10):
            data1 = pd.read_csv(path + filenames[i], sep=';')
            data2 = pd.read_csv(path + filenames[j], sep=';')
            MT1 = np.unique(data1[column_str])
            MT2 = np.unique(data2[column_str])
            n_duplicates+=write_equals(fname_results,MT1,MT2,filenames[i].split(token)[0],filenames[j].split(token)[0], dm, dt)
    print("Number of duplicates: {}".format(n_duplicates))


def equality(mt1, mt2, dm, dt):
    m1, t1 = np.float64(mt1.split('@'))
    m2, t2 = np.float64(mt2.split('@'))
    bm = np.abs(m1-m2)<=dm
    bt = np.abs(t1-t2)<=dt
    
    return bm & bt

def equality_full(mt1, mt2, dm, dt):
    m1, t1 = np.float64(mt1.split('@'))
    m2, t2 = np.float64(mt2.split('@'))
    dm_actual = np.abs(m1-m2)
    dt_actual = np.abs(t1-t2)
    bm = dm_actual<=dm
    bt = dt_actual<=dt
    
    return bm & bt, dm_actual, dt_actual

def print_equals(arr1, arr2, f1, f2, dm, dt):
    n_duplicates=0
    for mt1 in arr1:
        for mt2 in arr2:
            if(equality(mt1, mt2, dm, dt)):
                n_duplicates+=1
                print('{}:{}\n{}:{}\n-----------'.format(f1,mt1,f2,mt2))
    
    return n_duplicates

def print_equals_full(arr1, arr2, f1, f2, dm, dt):
    n_duplicates=0
    for mt1 in arr1:
        for mt2 in arr2:
            b, dm_actual, dt_actual = equality(mt1, mt2, dm, dt)
            if(b):
                n_duplicates+=1
                print('{}:{}\n{}:{}    dm={}, dt={}\n-----------'.format(f1,mt1,f2,mt2, dm_actual, dt_actual))
    
    return n_duplicates
    
def write_equals(fname, arr1, arr2, f1, f2, dm, dt):
    n_duplicates=0
    with open(fname,'a') as f:
        for mt1 in arr1:
            for mt2 in arr2:
                if(equality(mt1, mt2, dm, dt)):
                    n_duplicates+=1
                    f.write('{}:{}\n{}:{}\n-------------------\n'.format(f1,mt1,f2,mt2))
                    
    return n_duplicates
    

