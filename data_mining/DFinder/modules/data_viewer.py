import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def pair_string(a,b,inds,data):
    ind0,ind1=inds[0],inds[1]
    n_duplicates = np.sum(a==b)
    s = 'SIMILARITY: {} ({} out of 40)\n'.format(n_duplicates/len(a), n_duplicates)
    s += '----------------\n'
    for i in range(4):
        s += '{}\n'.format(data.iloc[ind0,i]) 
    s += '----------------\n'
    for i in range(4):
        s += '{}\n'.format(data.iloc[ind1,i]) 
    s += '----------------\n'
    return s

def hist_pair(a,b,inds,data):
    ind0,ind1=inds[0],inds[1]
    plt.hist(a)
    plt.hist(b)
    plt.legend([data.iloc[ind0,:1][0], data.iloc[ind1,:1][0]])
    plt.show()


def report_pair(a,b,inds,data):
    print(pair_string(a,b,inds,data))
    hist_pair(a,b,inds,data)
    
def group_hist(a,b):
    groups = [np.arange(0,40), np.arange(0,10), np.arange(10,20), np.arange(20,30), np.arange(30,40)]
    
    xmin, xmax = np.min([np.min(a),np.min(b)]), np.max([np.max(a),np.max(b)])
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].hist(a[groups[1]])
    axs[0, 0].hist(b[groups[1]])
    axs[0, 1].hist(a[groups[2]])
    axs[0, 1].hist(b[groups[2]])
    axs[1, 0].hist(a[groups[3]])
    axs[1, 0].hist(b[groups[3]])
    axs[1, 1].hist(a[groups[4]])
    axs[1, 1].hist(b[groups[4]])
    for i in range(2):
        for j in range(2):
            axs[i,j].set_xlim([xmin,xmax])
    
    plt.show()


