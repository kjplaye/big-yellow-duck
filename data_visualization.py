import os
import tempfile
import ctypes
import numpy as np

from imagesc import imagesc
from ggobi import *

my_path = os.path.dirname(os.path.abspath(__file__))
_wview = ctypes.cdll.LoadLibrary(my_path + '/_wview.so')

def plotn(x):
    """
    Hopefully we get a better version soon
    """
    if type(x) == str:
        os.system('plotn '+x)
    else:
        f = tempfile.NamedTemporaryFile()
        for L in x:
            f.write(' '.join(map(str,L)) + '\n')
        f.flush()
        os.system('plotn '+f.name)
        f.close()

def wview(x):
    a = np.require(x,dtype = np.float64, requirements = 'C')
    b = (a - a.min()) / (a.max() - a.min()) - 0.5
    bp = b.ctypes.get_as_parameter()
    _wview.wview(bp,ctypes.c_int64(len(b)))
    
