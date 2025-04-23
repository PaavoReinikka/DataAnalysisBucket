import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import chardet
import string
from collections import defaultdict


from nltk.corpus import stopwords as sw
from nltk.corpus import wordnet as wn
from nltk import wordpunct_tokenize
from nltk import WordNetLemmatizer
from nltk import sent_tokenize
from nltk import pos_tag
from nltk.probability import FreqDist
from nltk.stem import SnowballStemmer, PorterStemmer
import nltk
nltk.download('averaged_perceptron_tagger')
nltk.download('wordnet')
nltk.download('punkt') 


from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_extraction.text import TfidfVectorizer, TfidfTransformer, CountVectorizer
from sklearn.pipeline import Pipeline
from sklearn.cluster import KMeans, AgglomerativeClustering, \
                            SpectralClustering 
from sklearn.decomposition import TruncatedSVD, PCA
from sklearn.manifold import SpectralEmbedding, TSNE
from sklearn.preprocessing import LabelEncoder, normalize
from sklearn.metrics import pairwise_kernels,\
                            normalized_mutual_info_score

from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform, cdist


class Preprocessor(BaseEstimator, TransformerMixin):

    def __init__(self, stopwords=set(), punctuation=set(), min_length=1,
                 lower=True, strip=True, utf=True, 
                 numencode=False, stem=True, nonumbers=False, verbose=False):
                
        self.min_length = min_length
        self.utf = utf
        self.numencode = numencode
        self.stem = stem
        self.nonumbers = nonumbers
        self.lower      = lower
        self.strip      = strip
        self.stopwords  = stopwords | set(sw.words('english'))
        self.punct      = punctuation | set(string.punctuation)
        self.lemmatizer = WordNetLemmatizer()
        self.stemmer = SnowballStemmer('english')# PorterStemmer()
        
        if(verbose):
            print("Preprocessor...\nlower: {}".format(self.lower))
            print("strip: {}".format(self.strip))
            print("min_length: {}".format(self.min_length))
            print("utf: {}".format(self.utf))
            print("numencode: {}".format(self.numencode))
            print("stem: {}".format(self.stem))
            print("lemmatize: {}".format(not self.stem))
            print("nonumbers: {}".format(self.nonumbers))
        

    def fit(self, X, y=None):
        return self

    def inverse_transform(self, X):
        return [" ".join(doc) for doc in X]

    def transform(self, X):
        tmp = [
            list(self.tokenize(doc)) for doc in X
        ]
        return [" ".join(x) for x in tmp]

    def tokenize(self, document):
        
        for sent in sent_tokenize(document):
            # Every sent to token and pos tag (regardless of whether or not it is used)
            for token, tag in pos_tag(wordpunct_tokenize(sent)):
                if(self.utf):
                    if(len(token)!=len(token.encode('utf-8'))):
                        continue
                # Simple preprocessing of the tokens
                token = token.lower() if self.lower else token
                token = token.strip() if self.strip else token
                token = token.strip('_') if self.strip else token
                token = token.strip('*') if self.strip else token

                # ignore stopword
                if token in self.stopwords:
                    continue

                # ignore punctuation (change to all if necessary)
                if any(char in self.punct for char in token):
                    continue
                
                #choose stem or lemma
                if(self.stem):
                    result = self.stemmer.stem(token)
                else:
                    tag = {
                            'N': wn.NOUN,
                            'V': wn.VERB,
                            'R': wn.ADV,
                            'J': wn.ADJ
                        }.get(tag[0], wn.NOUN)
                    result = self.lemmatizer.lemmatize(token, tag)
                    #result = self.lemmatize(token, tag)
                
                # ignore stopword
                if result in self.stopwords:
                    continue
                
                #drop singular tokens
                if(len(result)<=self.min_length):
                    continue
                
                #drop any token including digit(s) (also drops all digits)
                if(self.nonumbers):
                    if any(char.isdigit() for char in result):
                        continue
                        
                #set numberflag if necessary
                if(self.numencode and result.isdigit()):
                    if(self.numencode=='drop'):
                        continue
                    result = 'numerolippu'
                
                yield result

    def lemmatize(self, token, tag):
        tag = {
            'N': wn.NOUN,
            'V': wn.VERB,
            'R': wn.ADV,
            'J': wn.ADJ
        }.get(tag[0], wn.NOUN)

        return self.lemmatizer.lemmatize(token, tag)

    
def top_cluster_features(countVec, labels, k=10):
    n = len(np.unique(labels))
    result = np.zeros((n,k), dtype=int)
    
    for label in range(n):
        x = countVec[np.where(labels==label),:][0]
        counts = np.sum(x,axis=0)
        ind = np.argsort(counts, axis=1)[:,-k:]
        result[label,:] = ind
    
    return result[:,::-1]
      
    
def plot_top_words(X, k=10, title=None):
    #FreqDist plot of top k words in X (e.g., a cluster)
    all_words = []
    
    for item in X:
        for word in item.split():
            all_words.append(word)
            
    fdist = FreqDist(all_words)
    if(title is not None):
        plt.title(title)
    fdist.plot(k,cumulative=False)

def plot_all_clusters(X, labels, k=10):
    for label in np.unique(labels):
        plot_top_words(np.asarray(X)[np.where(labels==label)],k,'Cluster {}'.format(label))

def table_top_words(cv, top_words):
    result=pd.DataFrame()
    for i in range(5):
        result['Cluster {}'.format(i)]= np.array(cv.get_feature_names())[top_words[i,:]]
    return result

    
def prune_dict(D, keys_to_prune):
    for key in keys_to_prune:
        del D[key]
    return D


def get_word_freqs(X, pruning=0):
    freq = defaultdict(int)
    for text in X:
        for token in text.split():
            freq[token] += 1
    
    n = len(X)
    if(pruning>0):
        if(pruning<1):
            pruning *= n
        keys = []
        for key in freq.keys():
            if(freq[key]<pruning):
                keys.append(key)
        prune_dict(freq, keys)
    
    return freq


def laplacian(S, normalize='rw'):
    if(normalize=='rw'):
        return np.diag(1./np.sum(S,axis=1)) @ (np.diag(np.sum(S,axis=1)) - S)
    elif(normalize=='sym'):
        return np.sqrt(np.diag(1./np.sum(S,axis=1))) @ (np.diag(np.sum(S,axis=1)) - S) @\
        np.sqrt(np.diag(1./np.sum(S,axis=1)))
    else:
        return np.diag(np.sum(S,axis=1)) - S

    
def spect_embedding(L, n_components=1, normalized=False):
    #Assumes a fully connected graph (i.e., drops the first [minor] eigen vec by default)
    evals, evecs = np.linalg.eig(L)
    index = sorted(range(evals.shape[0]), key=lambda k: evals[k] )
    res = evecs[:,index[1:n_components+1]]
    if(normalized):
        return normalize(res, axis=1)
    return res


def kmeans_objective(D, labels):
    obj = 0
    for label in np.unique(labels):
        D2 = D[np.where(labels==label)]**2
        obj += np.sum(D2)
    return obj/D.shape[0]


#THESE ARE FOR THE PARALLEL K-MEANS FUNCTIONS
#from numba.extending import overload <-> to extend to cluster digest (TODO)
'''
from numba import njit

@njit(parallel=True)
def kmeans(X, k, n_iter, init_centroids):
    #Fast parallel kmeans
    #Original implementation from old numba examples (here slightly modified)
    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids
    
    for l in range(n_iter):
        dist = np.array([[np.sqrt(np.sum((X[i, :] - centroids[j, :])**2))
                                for j in range(k)] for i in range(N)])

        predictions = np.array([dist[i, :].argmin() for i in range(N)])
            
        centroids = np.array([[np.sum(X[predictions == i, j])/np.sum(predictions == i)
                                 for j in range(D)] for i in range(k)])

    return centroids, dist, predictions



@njit(parallel=True)
def kmeans_spherical(X, k, n_iter, init_centroids):
    #Fast parallel spherical kmeans
    #Original implementation from old numba examples (here slightly modified)
    N = X.shape[0]
    D = X.shape[1]
    centroids = init_centroids
    
    for l in range(n_iter):
        dist = X@centroids.T
        
        predictions = np.array([dist[i, :].argmin() for i in range(N)])
            
        centroids = np.array([[np.sum(X[predictions == i, j])/np.sum(predictions == i)
                                 for j in range(D)] for i in range(k)])

    return centroids, dist, predictions
'''


def run_kmeans_experiments(X, y_true, repeats=10, verbose=False):
    seed_set = np.random.randint(0,1e5, repeats)
    max_so_far = 0
    ind=-1
    k=0
    all_nmis = []

    for seed in seed_set:
        kmeans = KMeans(n_clusters=5,init='k-means++',random_state=seed).fit(X)
        nmi=normalized_mutual_info_score(y_true,kmeans.labels_)
        all_nmis.append(nmi)
        
        if nmi>max_so_far:
            ind=k
            max_so_far=nmi
        k+=1
        
    average = np.mean(all_nmis)
    #stdev = np.std(all_nmis)
    seed = seed_set[ind]
    
    if(verbose):
        print('max: {}'.format(max_so_far))
        print('average: {}'.format(average))
        #print('std: {}'.format(stdev))
        print('seed: {}'.format(seed))
    
    return max_so_far, average, seed


def buckshot_init(data, k=5, sample_size=None, link='ward', 
                  metric='euclidean', labels=None):
    N = data.shape[0]
    
    if(sample_size is None):
        sample_size = int(np.sqrt(k*N))
    
    ind = np.random.choice(range(N), sample_size, replace=False)
    X = data[ind,:]
    cd_matrix = pdist(X, metric=metric)
    
    link = hierarchy.linkage(cd_matrix, method=link,  optimal_ordering=True)
    plt.figure(figsize=(6,6))
    hierarchy.dendrogram(link, truncate_mode='lastp', p=5)
    predictions = hierarchy.fcluster(link, t=5, criterion='maxclust')
    
    if(labels is not None):
        print('NMI of buckshot: {}'.format(normalized_mutual_info_score(labels[ind],
                                                                        predictions)))
    
    return X[sample_cluster_centers(predictions),:]


def sample_cluster_centers(labels):
    return  np.array([
            np.random.choice(np.where(labels==1)[0],1, replace=False)[0],\
            np.random.choice(np.where(labels==2)[0],1, replace=False)[0],\
            np.random.choice(np.where(labels==3)[0],1, replace=False)[0],\
            np.random.choice(np.where(labels==4)[0],1, replace=False)[0],\
            np.random.choice(np.where(labels==5)[0],1, replace=False)[0]
    ])



