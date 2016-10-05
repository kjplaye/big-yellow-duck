import numpy as np
from ctypes import *
from basic_util import dim
import os

my_path = os.path.dirname(os.path.abspath(__file__))
_imagesc = cdll.LoadLibrary(my_path + '/_imagesc.so')

def imagesc(L):
    """
    Example:
    >>> import math
    >>> M = [[math.sin(x*y*0.1) for x in xrange(100)] for y in xrange(100)]
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
    if dim(L) == 2:
        x = np.require(L,dtype = np.float64, requirements = 'C')
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x.size),c_int64(x.shape[1]),c_int64(1))
    elif dim(L) == 3:
        x = np.require(L,dtype = np.float64, requirements = 'C')
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x[0].size),c_int64(x.shape[2]),c_int64(x.shape[0]))
    else:        
        raise ValueError('Need a 2-d object data[y,x] or 3-d object data[frame,y,x]')

