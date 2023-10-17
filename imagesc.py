import numpy as np
from ctypes import *
from basic_util import dim
import os
import multiprocessing as mp


my_path = os.path.dirname(os.path.abspath(__file__))
_imagesc = cdll.LoadLibrary(my_path + '/_imagesc.so')

def _reshape_nan(M, new_shape):
	new_M = np.full(new_shape, np.nan)
	new_M[:M.shape[0],:M.shape[1]] = M
	return new_M

def _unjag(M):
    try:
    	X = np.require(e, dtype = np.float64, requirements = 'C')
    	return X
    except:
    	max_len = max([len(e) for e in M])
    	Mp = [np.concatenate([e, [np.nan] * (max_len - len(e))]) for e in M]
    	X = np.require(Mp, dtype = np.float64, requirements = 'C')
    	return X

def _do_imagesc(*args):
    L = args[0]
    window_name_bytes = args[1].encode() if len(args) >= 2 else 'imagesc'.encode()
    window_name = c_char_p(window_name_bytes)
    _imagesc.imagesc.argtypes = (c_void_p, c_int64, c_int64, c_int64, c_char_p)
    if dim(L) == 2:
        x = _unjag(L)
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x.size),c_int64(x.shape[1]),c_int64(1),window_name)
    elif dim(L) == 3:
        X = [_unjag(x) for x in L]
        max_shape = np.max([x.shape for x in X], 0)
        x = np.require([_reshape_nan(xx, max_shape) for xx in X], dtype = np.float64, requirements = 'C')
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x[0].size),c_int64(x.shape[2]),c_int64(x.shape[0]),window_name)
    else:        
        raise ValueError('Need a 2-d object data[y,x] or 3-d object data[frame,y,x]')

def imagesc(*L):
    """
    Example:
    >>> import math
    >>> M = [[[math.sin(x*y*t/100) for x in range(100)] for y in range(100)] for t in range(10,100)]
    >>> imagesc(M)

    KEYS:
       Q             : Quit
       S             : Zoom Standard
       C             : Change Color Mode
       PgUp or PgDn  : Horizontal Zoom
       Home or End   : Vertical Zoom
       + or -        : Change Color Intensity
       [ or ]        : Change Width
       , or .        : Change Frames
       Arrow Keys    : Pan
       Return        : Print data index and value
       Mouse Buttons : Zoom     
    """
    p = mp.Process(target = _do_imagesc, args = L)
    p.start()
