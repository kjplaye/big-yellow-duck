import numpy as np
import scipy.stats
from data_visualization import ggobi
from math import log
from matplotlib import pyplot as plt
from basic_util import col
from sklearn.decomposition import FastICA

try:
    from mojave_eda import mojave
except:
    print("No Mojave")

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
    return -scipy.stats.norm.logsf(zscore) / log(2.0), mean, zscore

def _pca_bits(m,n,svs):
    ans = []
    zscore_ans = []
    svs_array = np.array(svs)
    for k in range(len(svs)):
        rm_sigma = np.sum(svs_array[k:]**2 / (m*n))**0.5
        (bits, mean0, zscr0) = _max_sv_bits(m-k,n-k,svs[k] / rm_sigma)
        if k==0:
            mean = mean0
        ans.append(bits)
        zscore_ans.append(zscr0)
    return np.array(ans), mean, np.array(zscore_ans)

class ICA(np.ndarray):
    def __new__(self, input_array):
    	"""Make sure to PCA before this call and truncate to desired size.

        Matrix decomposition: 
            self.M = I @ I.new_V
        """
    	FI = FastICA()
    	sources = FI.fit_transform(np.array(input_array))
    	obj = np.asarray(sources).view(self)
    	obj.A = FI.mixing_   	 
    	return obj
    def mojave(self, input_clusters = None, dims = 20):
        if input_clusters is None:
            return mojave(self[:,:dims])
        else:
            return mojave(self[:,:dims], input_clusters)


class PCA:
    """PCA
    Let M be a matrix with shape (m,n).  m is the amount of data
    and n is the number of features.  Usually m > n. Matrix 
    decomposition:

        P.M = P.U @ P.D @ P.V

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
    
    ind  MP  PV     sv   zscr     bits
    ---  --  --     --   ----     ----
    0 (M) (P) 102.53 251.72 45717.25
    1          29.89  -3.51     0.00
    2          29.68  -3.01     0.00
    3          29.37  -2.90     0.00
    4          29.18  -2.34     0.01
    5          29.00  -1.72     0.06
    6          28.77  -1.28     0.15
    7          28.59  -0.68     0.41
    8          28.22  -0.85     0.32
    9          27.92  -0.75     0.37
    """
    def __init__(self, M, remove_mean = False):
        M = np.array(M)
        self.M = M
        if remove_mean:
                self.M -= np.mean(self.M,0)
        [self.U, self.D, self.V] = np.linalg.svd(self.M,0)
        [self.m, self.n] = self.M.shape
        self.pval, self.max_sv_mean, self.zscore = _pca_bits(self.m, self.n, self.D)
    def get_mp_threshold(self):
    	return self.max_sv_mean * np.sum(self.D**2 / (self.m * self.n))**0.5
    def estimated_rank(self, threshold_bits = 20.0, use_mp_thresh = True):
        if use_mp_thresh:
            return np.sum(self.D > self.get_mp_threshold())
        else:
            return np.argmin(self.pval >= threshold_bits)
    def plot_diag(self):
    	plt.plot(self.D,'o', color = 'blue')
    	plt.plot(self.D,'.', color = 'red')
    	plt.plot([self.get_mp_threshold()]*len(self.D), color = 'green')
    	plt.ylabel("Singular Values")
    	plt.xlabel("Rank")
    def ggobi(self, input_clusters = None, dims = 20):
    	if input_clusters is None:
            return ggobi(self.U[:,:dims])
    	else:
            return ggobi(self.U[:,:dims], input_clusters)
    def mojave(self, input_clusters = None, dims = 20):
    	if input_clusters is None:
            return mojave(self.U[:,:dims])
    	else:
            return mojave(self.U[:,:dims], input_clusters)
    def ica(self, dims):
        I = ICA(self.U[:,:dims])
        I.new_V = I.A.T @ (np.diag(self.D) @ self.V)[:len(I.A)]
        return I
    def __repr__(self, threshold_bits = 20.0):
        mp_thresh = self.get_mp_threshold()
        header = ['ind','MP','PV','sv','zscr','bits']
        header_lines = ['---','--','--','--','----','----']
        L = [[k,
              '(M)' if self.D[k] > mp_thresh else '   ',
              '(P)' if self.pval[k] > threshold_bits else '   ',
              "%6.2f" % self.D[k],
              "%6.2f" % self.zscore[k],
              "%6.2f" % self.pval[k]] for k in range(min(10,len(self.D)))]
        est_rank = self.estimated_rank()
        LH = [header] + [header_lines] + L
        begin = f'PCA({self.m}, {self.n})\nestimated rank = ' \
            f'{self.estimated_rank()}, mean = {"%6.2f" % self.max_sv_mean}\n\n'
        return begin + str(col(LH))
