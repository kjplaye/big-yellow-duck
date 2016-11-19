import numpy as np
import ctypes
import multiprocessing as mp
import os

my_path = os.path.dirname(os.path.abspath(__file__))
_wview = ctypes.cdll.LoadLibrary(my_path + '/_wview.so')

def _do_wview(*x):
    a = np.require(x,dtype = np.float64, requirements = 'C')
    b = (a - a.min()) / (a.max() - a.min()) - 0.5
    bp = b.ctypes.get_as_parameter()
    _wview.wview(bp,ctypes.c_int64(len(b)))

def wview(x):
    p = mp.Process(target = _do_wview, args = x)
    p.start()

