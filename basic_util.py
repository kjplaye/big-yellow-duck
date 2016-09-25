import math
from cols import transpose
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

def form_table(L):
    VL = dict()
    VR = dict()
    l_cnt = 0
    r_cnt = 0
    for l in L:
        if l[0] not in VL:
            VL[l[0]] = l_cnt
            l_cnt+=1
        if l[1] not in VR:
            VR[l[1]] = r_cnt
            r_cnt+=1
    T = [[0 for j in range(r_cnt)]  for i in range(l_cnt)]
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
