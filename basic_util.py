from collections import Counter
import math
from cols import transpose, col
import pickle
import random
import numpy as np
import scipy.stats
import tempfile
import sys
import os
import pickle

def dim(x):
    if isinstance(x,str):
        return 0
    try:
        return dim(x[0]) + 1
    except:
        return 0

def blockup(x,n):
    return [x[i:i+n] for i in xrange(0,len(x)-n+1,n)]


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
    return [x[i+d] - x[i] for i in xrange(len(x) - d)]

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
  def labeled(self):
    m = len(self)
    n = len(self[0])
    return col([['',''] + [self.colkeys[j] for j in range(n)]] + [['',''] + ['-' for j in range(n)]] + [[self.rowkeys[i],':'] + [self[i][j] for j in range(n)] for i in range(m)])

def form_table(L,rowkeys = None,colkeys = None,counts_given = False):
    if rowkeys == None:
      SL = set()
      for l in L:
        if l[0] not in SL:
          SL.add(l[0])
      rowkeys = sorted(list(SL))
    VL = {rowkeys[i]:i for i in range(len(rowkeys))}
    if colkeys == None:
      SR = set()
      for l in L:
        if l[1] not in SR:
          SR.add(l[1])
      colkeys = sorted(list(SR))
    VR = {colkeys[i]:i for i in range(len(colkeys))}
    T = [[0 for j in range(len(colkeys))]  for i in range(len(rowkeys))]
    if counts_given:
      for l in L:
          T[VL[l[0]]][VR[l[1]]] += int(l[2])
    else:
      for l in L:
          T[VL[l[0]]][VR[l[1]]] += 1
    TAB = form_table_result(T)
    TAB.rowkeys = rowkeys
    TAB.colkeys = colkeys
    return TAB

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

def emacs(x):
  f = tempfile.NamedTemporaryFile()
  f.write(str(x))
  f.flush()
  os.system('emacs %s &' % f.name)
  sys.stdin.readline()

def pick_save(file_name,x):
  open(file_name,'w').write(pickle.dumps(x))

def pick_load(file_name):
  return pickle.loads(open(file_name).read())

def h_pad(L,padn = 5, padv = -1):
  n = L[0].shape[1]
  L0 = sum([[L[i],padv + np.zeros([padn,n])] for i in range(len(L)-1)],[]) + [L[-1]]
  return np.concatenate(L0)

def v_pad(L,padn = 5, padv = -1):
  return h_pad([e.T for e in L],padn = padn, padv = padv).T
