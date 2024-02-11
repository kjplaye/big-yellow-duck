from collections import Counter
from cols import transpose, col
import math
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd
import pickle
import pickle
import random
import scipy.stats
import sys
import tempfile


def dim(x):
    if isinstance(x,str):
        return 0
    try:
        return dim(x[0]) + 1
    except:
        return 0

def blockup(x, n, pad = False, padv = None):
    if pad:
    	jag = [x[i:i+n] for i in range(0,len(x),n)]
    	jag[-1] = list(jag[-1]) + [padv] * (n - len(jag[-1]))
    	return jag
    return [x[i:i+n] for i in range(0,len(x)-n+1,n)]

def rotate(x,n):
    return x[n:] + x[:n]

def reverse_it(x):
    if isinstance(x,str):
        return ''.join(reverse_it(list(x)))
    else:
        L = x[:]
        L.reverse()
        return L

def delta(x,d = 1):
    return [x[i+d] - x[i] for i in range(len(x) - d)]

def flatten(x):
    if dim(x)==2:
        return [e for l in x for e in l]
    elif dim(x) > 2:
        return flatten([flatten(e) for e in x])
    else:
        return x

def sortby(x,L):
    E = sorted([[x[i],L[i]] for i in range(len(x))],cmp = lambda a,b: 1-2*(a[1] < b[1]))
    return([e[0] for e in E])

def shuffle_it(x):
    y = x[:]
    random.shuffle(y)
    return y

def dither(x,e = 0.5):
  a = np.array(x)
  n = np.random.random(a.shape) * e
  return a + n

def pval(x):
    if hasattr(x,'pval'):
        return x.pval()
    if hasattr(x,'bits'):
        return 2**-x.bits()

def bits(x):
    if hasattr(x,'bits'):
        return x.bits()
    if hasattr(x,'pval'):
        return -math.log(x.pval(),2)


def entropy(x):
    c = dict()
    T = 0
    for e in x:
        c[e] = c.get(e,0) + 1
        T+=1
    P = [float(c[k])/T for k in c]
    return -sum([math.log(p,2) * p for p in P])

class form_table_result(col):
    def pval(self):
        return scipy.stats.chi2_contingency(self)[1]
    def bits(self):
        x, p, dof, ex = scipy.stats.chi2_contingency(self)
        return -scipy.stats.chi2.logsf(x,dof) / math.log(2)
    def chi_sqr_components(self):
        x, p, dof, ex = scipy.stats.chi2_contingency(self)
        a = np.array(self)
        return (a - ex)**2 / ex
    def labeled(self):
        m = len(self)
        n = len(self[0])
        return col([['',''] + [self.colkeys[j] for j in range(n)]] + [['',''] + ['-' for j in range(n)]] + [[self.rowkeys[i],':'] + [self[i][j] for j in range(n)] for i in range(m)])
                  
        
def form_table(L,rowkeys = [],colkeys = []):
    rowkeys = set(rowkeys)
    colkeys = set(colkeys)
    for e in L:
        rowkeys.add(e[0])
        colkeys.add(e[1])
    rowkeys = sorted(list(rowkeys))
    colkeys = sorted(list(colkeys))
    rowdict = {rowkeys[i]:i for i in range(len(rowkeys))}
    coldict = {colkeys[i]:i for i in range(len(colkeys))}
    ans = [[0 for j in range(len(colkeys))] for i in range(len(rowkeys))]
    for e in L:
        ans[rowdict[e[0]]][coldict[e[1]]] += 1
    res = form_table_result(ans)
    res.rowkeys = rowkeys
    res.colkeys = colkeys
    return res

    
class myhist_result(col):
    def pval(self):
        C = [e[1] for e in self]
        x,p = scipy.stats.chisquare(C)
        return p
    def bits(self):
        C = [e[1] for e in self]
        x,p = scipy.stats.chisquare(C)
        return -scipy.stats.chi2.logsf(x,len(C)-1)/math.log(2)

def myhist(x,sortit = True):
  C = Counter(x)
  KV = [[k,C[k]] for k in C.keys()]
  if sortit:
    KV = sorted(KV,key = lambda x:x[1])
  return myhist_result(KV)

def pick_save(file_name,x):
  open(file_name,'wb').write(pickle.dumps(x))

def pick_load(file_name):
  return pickle.loads(open(file_name,'rb').read())

def h_pad(L,padn = 5, padv = -1):
  n = L[0].shape[1]
  L0 = sum([[L[i],padv + np.zeros([padn,n])] for i in range(len(L)-1)],[]) + [L[-1]]
  return np.concatenate(L0)

def v_pad(L,padn = 5, padv = -1):
  return h_pad([e.T for e in L],padn = padn, padv = padv).T

def lmap(f,L):
    return list(map(f,L))

def make_patches(colors, labels):
    """USAGE: make_patches([[0,0,1], [0,1,0], [1,0,0]], ['10 trials', '100 trials', '1000 trials'])"""
    patch = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(colors))]
    plt.legend(handles=patch)

class Cat(np.ndarray):
    """Expose some pandas functionality to consistently index categorical data.
    
    Usage Example:
    >>> cat = Cat(['a','c','b','a'])
    >>> print(cat)
    [0, 2, 1, 0]

    >>> print("Number of categories:", cat.num_cat)
    Number of categories: 3

    >>> print("Categories:", cat.cat)
    Categories: Index(['a', 'b', 'c'], dtype='object')

    >>> new_cat = cat.index(['c','a'])
    >>> print("Pull in new data:", new_cat))
    >>> print("As a series:")
    >>> print(new_cat.series)
    >>> print("Lift to names:", new_cat.lift())
    Pull in new data: [2 0]
    As a series:
    0   c
    1   a
    dtype: category
    Categories (3, object): ['a', 'b', 'c']
    Lift to names: Index(['c', 'a'], dtype='object')
    """
    def __new__(self, data, categories = None):
        cat = pd.Categorical(data, categories = categories)
        series = pd.Series(cat)
        codes = series.cat.codes
        obj = np.asarray(codes.values).view(self)
        obj.series = series
        obj.cat = cat.categories
        obj.num_cat = obj.cat.nunique()
        return obj
    def index(self, data):
        """Apply Cat indexing to new data."""
        return Cat(data, categories=self.cat)
    def lift(self):
        return self.cat[self]
    def __hash__(self):
        return hash((tuple(self), tuple(self.cat)))

def compose_advice(list_of_arrays, array_names = None):
    """
    >>> compose_advice([[[[1,2],[4,5]]], [0,1,1,1,0,1,0,0,0,0]], ['count','meow'])
    count[:meow:]
    count[::meow]
    """
    if array_names is None:
        array_names = [f'x_{i}' for i in range(len(list_of_arrays))]
    L = [np.array(e) for e in list_of_arrays]
    SIZE = [e.max() + 1 for e in L]
    SHAPE = [e.shape for e in L]
    for array_num_1 in range(len(SHAPE)):
        for dim in range(len(SHAPE[array_num_1])):
            for array_num_2 in range(len(SIZE)):
                if SHAPE[array_num_1][dim] == SIZE[array_num_2]:
                    a = ':'.join([f'{array_names[array_num_2]}' if i==dim else ''
                                  for i in range(len(SHAPE[array_num_1]))])
                    print(f"{array_names[array_num_1]}[{a}]")
