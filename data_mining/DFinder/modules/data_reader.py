import pandas as pd
import numpy as np

def read_file(file, path=None):
	if path is None:
		data = pd.read_csv(file, sep=';', header=None).drop([0,1,2],axis=0).drop([2,3,6,7],axis=1)
	else:
		data = pd.read_csv(path + file, sep=';', header=None).drop([0,1,2],axis=0).drop([2,3,6,7],axis=1)
	data=data.iloc[np.argsort(np.asarray(data[4],dtype=float)),:]
	return data
	
def array_pair(data, inds, start = 4):
	a = np.asarray(data.iloc[inds[0],start:], dtype=float)
	b = np.asarray(data.iloc[inds[1],start:], dtype=float)
	return a, b

