import pandas as pd
import numpy as np
import PySimpleGUI as psg

def get_MTime_duplicates(data, dm, dt):
    indices = []
    mt_pairs = []
    similarities = ['' for i in range(40)]
    mass = np.asarray(data[4], dtype=float)
    time = np.asarray(data[5], dtype=float)
    for i in range(data.shape[0]-1):
        for j in range(i+1,data.shape[0]):
            if(np.abs(mass[i]-mass[j])<=dm and np.abs(time[i]-time[j])<=dt):
                mt1 = '{}@{}'.format(mass[i],time[i])
                mt2 = '{}@{}'.format(mass[j],time[j])
                mt_pairs.append((mt1,mt2))
                indices.append((i,j))
    return indices, mt_pairs, similarities


def get_value_duplicates(data, tol):
    indices = []
    mt_pairs = []
    similarities = []
    
    size = data.shape[0]-1
    counter=0
    for i in range(data.shape[0]-1):
        counter += 1
        if not psg.one_line_progress_meter('Progress Meter', counter, size,'Duplicate finder', no_button=False):
            exit()
        
        a = np.asarray(data.iloc[i,4:], dtype=float)
        for j in range(i+1,data.shape[0]):
            b = np.asarray(data.iloc[j,4:], dtype=float)
            similarity = np.sum(a==b)/len(a)
            if(similarity >= tol):
                indices.append((i,j))
                similarities.append('  ({}%)'.format(similarity))
                mt1 = '' + str(data.iloc[i,2]) + '@' + str(data.iloc[i,3])
                mt2 = '' + str(data.iloc[j,2]) + '@' + str(data.iloc[j,3])
                mt_pairs.append((mt1,mt2))
    
    return indices, mt_pairs, similarities
    
    
def get_VAL_duplicates(data, tol):
    indices = []
    mt_pairs = []
    
    for i in range(data.shape[0]-1):
        a = np.asarray(data.iloc[i,4:], dtype=float)
        for j in range(i+1,data.shape[0]):
            b = np.asarray(data.iloc[j,4:], dtype=float)
            if(np.sum(a==b)/len(a) >= tol):
                indices.append((i,j))
                mt1 = '' + str(data.iloc[i,2]) + '@' + str(data.iloc[i,3])
                mt2 = '' + str(data.iloc[j,2]) + '@' + str(data.iloc[j,3])
                mt_pairs.append((mt1,mt2))
    
    return indices, mt_pairs


