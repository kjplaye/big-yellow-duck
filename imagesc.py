import numpy
from ctypes import *
from basic_util import dim
from my_path import my_path

_imagesc = cdll.LoadLibrary(my_path + '_imagesc.so')
def imagesc(L):
    """
    Example:
    >>> import math
    >>> M = [[math.sin(x*y*0.1) for x in xrange(100)] for y in xrange(100)]
    >>> imagesc(M)
    """
    if dim(L) == 2:
        x = numpy.array(L,dtype = numpy.float64)
        _imagesc.imagesc(x.ctypes.get_as_parameter(),c_int64(x.size),c_int64(x.shape[1]))
    else:
        raise ValueError('Need a 2-d object')

