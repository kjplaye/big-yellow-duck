import numpy as np
import scipy.stats
from data_visualization import ggobi
from math import log
from matplotlib import pyplot as plt

# Regression model using (sqrt(m), sqrt(n), 1) for the max sv
# of a random (m,n) matrix with iid entries from N(0,1)
_mean_regress = np.array([1.02247417,  1.0224104 , -0.93166985])

# Pval bits for max_sv for a random (m,n) matrix with iid entries from N(0,1)
def _max_sv_bits(m, n, max_sv):
    mean = _mean_regress @ np.array([m**0.5, n**0.5, 1])
    mp_bound = m**0.5 + n**0.5
    if mp_bound >= 8:
    	sigma = 0.25
    elif mp_bound >= 5:
    	sigma = 0.2
    else:
    	sigma = 0.1
    zscore = (max_sv - mean) / sigma
    return -scipy.stats.norm.logsf(zscore) / log(2.0), mean

def _pca_bits(m,n,svs):
    ans = []
    svs_array = np.array(svs)
    for k in range(len(svs)):
        rm_sigma = np.sum(svs_array[k:]**2 / (m*n))**0.5
        (bits, mean0) = _max_sv_bits(m-k,n-k,svs[k] / rm_sigma)
        if k==0:
            mean = mean0
        ans.append(bits)
    return np.array(ans), mean

class pca:
    """PCA
    Let M be a matrix with shape (m,n).  m is the amount of data
    and n is the number of features.  Usually m > n.

    USAGE EXAMPLE
    >>> M = np.random.normal(0,1,[128,400])
    >>> M[0][0] = 100
    >>> P = pca(M, remove_mean = False)
    >>> P.plot_diag()
    >>> cluster = P.ggobi()
    >>> np.allclose(P.U @ np.diag(P.D) @ P.V, M)
    >>> print(P)
    PCA(128, 400)
    estimated rank = 1, mean =  31.08
    ----------------------------------------
    s_0 102.54 46158.92 bits
        ............
    s_1  30.53   0.55 bits
    s_2  30.22   0.72 bits
    s_3  29.96   1.12 bits
    s_4  29.65   1.31 bits
    """
    def __init__(self, M, remove_mean = False):
    	self.M = M
    	if remove_mean:
        	self.M -= np.mean(self.M,0)
    	[self.U, self.D, self.V] = np.linalg.svd(self.M,0)
    	[self.m, self.n] = self.M.shape
    	self.pval, self.max_sv_mean = _pca_bits(self.m, self.n, self.D)
    def estimated_rank(self, threshold_bits = 20.0):
    	return np.argmin(self.pval >= threshold_bits)
    def plot_diag(self):
    	plt.plot(self.D,'o', color = 'blue')
    	plt.plot(self.D,'.', color = 'red')
    	plt.plot([self.max_sv_mean]*len(self.D), color = 'green')
    	plt.ylabel("Singular Values")
    	plt.xlabel("Rank")
    def ggobi(self, input_clusters = None, dims = 20):
    	if input_clusters is None:
        	return ggobi(self.U[:,:dims])
    	else:
        	return ggobi(self.U[:,:dims], input_clusters)
    def __repr__(self):
        header_bar = '-' * 40
        middle_bar = '    ............'
        L = ['s_%d %6.2f %6.2f bits' % (k, self.D[k], self.pval[k]) for k in range(min(10,len(self.D)))]
        est_rank = self.estimated_rank()
        L = [header_bar] + L[:est_rank] + [middle_bar] + L[est_rank:]
        return f'PCA({self.m}, {self.n})\nestimated rank = {self.estimated_rank()}, mean = {"%6.2f" % self.max_sv_mean}\n' + '\n'.join(L)
   	 
   	 
