import math
from cols import transpose, col
import pickle
import random
import numpy as np

def dim(x):
    if isinstance(x,str):
        return 0
    try:
        return dim(x[0]) + 1
    except:
        return 0

def blockup(x,n):
    return [x[i:i+n] for i in xrange(0,len(x)-n+1,n)]

def entropy(x):
    c = dict()
    T = 0
    for e in x:
        c[e] = c.get(e,0) + 1
        T+=1
    P = [float(c[k])/T for k in c]
    return -sum([math.log(p,2) * p for p in P])

def rotate(x,n):
    return x[n:] + x[:n]

def reverse_it(x):
    if isinstance(x,str):
        return ''.join(reverse_it(list(x)))
    else:
        L = x[:]
        L.reverse()
        return L

def form_table(L,rowkeys = None,colkeys = None):
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
    for l in L:
        T[VL[l[0]]][VR[l[1]]] += 1
    return T

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
    return transpose(sorted(transpose([L,x])))[1]

def shuffle_it(x):
    y = x[:]
    random.shuffle(y)
    return y

def dither(x,e = 0.5):
  a = np.array(x)
  n = np.random.random(a.shape) * e
  return a + n
