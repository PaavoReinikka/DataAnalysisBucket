import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import argparse
import os
import yaml
import sys
sys.path.insert(0, './utilities/')
from MF import *

def get_features(masstime, data, significant=True):
    if not significant:
        return data
    inds=[0,1,2]
    for i in range(3,data.shape[0]):
        mt=data.iloc[i,4] + '@' + data.iloc[i,5]
        if mt in masstime:
            inds.append(i)
    inds=np.unique(np.asarray(inds, dtype=int))
    return data.iloc[inds,:]


def get_ids(data):
    ID=[]
    for line in data.iloc[2,8:]:
        ID.append(line.split('_')[-1].split('.')[0])
    return np.array(ID)


# parse the arguments. Only one is required, the input yaml file
parser = argparse.ArgumentParser(description='Plot PCA')
parser.add_argument('input', type=str, help='yaml config file', default=None)

# read the input file with yaml
args = yaml.load(open(parser.parse_args().input), Loader=yaml.FullLoader)
filename_significant = args['filename_significant']
filename_data = args['filename_data']
folder_out = args['folder_out']
plot_format = args['plot_format']
title = args['title']
figsize = args['figsize']
method = args['method']
INCLUDE_ID = args['includeID']
SAVE = folder_out is not None

# check if output file is given that the folder exists
# if not create it
if folder_out is not None:
    if not os.path.exists(os.path.dirname(folder_out)):
        os.makedirs(os.path.dirname(folder_out))

if figsize is None:
    figsize = (13,13)
else:
    figsize = (int(figsize), int(figsize))

if plot_format is None:
    plot_format = 'eps'

# Parse title
if title == 'None':
    # check if system is using '/' or '\'
    if '\\' in filename_significant:
        tle = filename_significant.split('\\')[-1]
    else:
        tle = filename_significant.split('/')[-1]
    tle = tle.split('ALPHA')[0]
else:
    tle = title


# read the data
data = pd.read_csv(filename_data, sep=';', header=None)
df_significance = pd.read_csv(filename_significant, sep=';')

masstime = np.unique(df_significance[' masstime'].to_numpy())
significant_features = get_features(masstime, data, True)
x1 = significant_features.iloc[3:,8:48].to_numpy(dtype=float).T
ID=get_ids(data.iloc[:,:48])

# print information
print()
print("input file: {}".format(filename_data))
print("significant features file: {}".format(filename_significant))
print("plot title: {}".format(tle))
print("{} significant features".format(len(masstime)))
print()


U_init, V_init = nnd_svd_initialization(x1, 3)
# ratings, U_init, V_init, n_iter=100, n_iter_inner=10, alpha=0.0002, beta=0.02, tol=1e-5
if method.lower() == 'mf':
    U, V, err = matrix_factorization(x1, U_init, V_init, n_iter=1000, n_iter_inner=10, alpha=0.0002, beta=0.02, tol=1e-5)
elif method.lower() == 'nnmf':
    U, V, err = nn_matrix_factorization(x1, U_init, V_init, n_iter=10000, tol=1e-5)
elif method.lower() == 'pca':
    pca = PCA().fit_transform(x1)
    U = pca[:,:3];err=None

if method != 'PCA':
    # check if the algorithm converged by comparing n_iter and err
    print('n_iter: {}'.format(len(err)))
    if len(err) < 1000:
        print('Converged')

assert(x1.shape[0]==40==U.shape[0])

# plot the first tree PCs
for k in range(2):
    for kk in range(k+1,3):
        fig, ax = plt.subplots(figsize=figsize)
        ax.scatter(U[:10,k],U[:10,kk], marker='^', c='red', s=100)
        ax.scatter(U[10:20,k],U[10:20,kk], marker='o', c='blue', s=100)
        ax.scatter(U[20:30,k],U[20:30,kk], marker='v', c='green', s=100)
        ax.scatter(U[30:,k],U[30:,kk], marker='s', c='gold', s=100)
        ax.legend(['aSYN','comb.','IFNg','UT'], fontsize=15)
        
        if INCLUDE_ID:
            for i, txt in enumerate(ID):
                ax.annotate(txt, (U[i,k], U[i,kk]),fontsize=12)
        
        plt.xlabel('LAT {}'.format(k+1), fontsize=25)    
        plt.ylabel('LAT {}'.format(kk+1), fontsize=25)
        plt.title(tle, fontsize=20)
        
        if SAVE:
            tle2 = tle + '_LAT{}vs{}'.format(k+1,kk+1)
            filename_output2 = folder_out + tle2
            print('Saving plot to {}'.format(filename_output2))
            # pick the format
            if plot_format == 'eps':
                fig.savefig(filename_output2 + '.eps')
            elif plot_format == 'tif':
                fig.savefig(filename_output2 + '.tif')
            elif plot_format == 'jpeg':
                fig.savefig(filename_output2 + '.jpg')
            else:
                print('Unknown format, use eps, tif or jpeg')
        else:
            plt.show()
            
            
    

