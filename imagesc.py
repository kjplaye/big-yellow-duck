import numpy as np
from ctypes import *
from basic_util import dim
import os
import multiprocessing as mp

my_path = os.path.dirname(os.path.abspath(__file__))
_imagesc = cdll.LoadLibrary(my_path + '/_imagesc.so')

def _do_imagesc(*L):
    if dim(L) == 2:
        x = np.require(L,dtype = np.float64, requirements = 'C')
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x.size),c_int64(x.shape[1]),c_int64(1))
    elif dim(L) == 3:
        x = np.require(L,dtype = np.float64, requirements = 'C')
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x[0].size),c_int64(x.shape[2]),c_int64(x.shape[0]))
    else:        
        raise ValueError('Need a 2-d object data[y,x] or 3-d object data[frame,y,x]')

def imagesc(L):
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
