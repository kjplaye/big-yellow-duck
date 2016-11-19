import numpy as np
import ctypes
import multiprocessing as mp
import os

my_path = os.path.dirname(os.path.abspath(__file__))
_wview = ctypes.cdll.LoadLibrary(my_path + '/_wview.so')

def _do_wview(*x):
    a = np.require(x,dtype = np.float64, requirements = 'C')
    ap = a.ctypes.get_as_parameter()
    s = a.shape
    print(s)
    if len(s) == 1:
        frames = 1
    elif len(s) == 2:
        frames = s[0]
    else:
        raise ValueError('Dimension must be 1 or 2')
    print(ap)
    _wview.wview(ap,ctypes.c_int64(s[-1]),ctypes.c_int32(frames))

def wview(x):
    p = mp.Process(target = _do_wview, args = x)
    p.start()

