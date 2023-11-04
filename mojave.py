from ctypes import *
from subprocess import Popen, PIPE
import numpy as np
import os


my_path = os.path.dirname(os.path.abspath(__file__))
_mojave = cdll.LoadLibrary(my_path + '/_mojave.so')

def mojave(X, cl = None, name = 'Mojave'):
    """
    Parameters
    ----------
    X : array_like
        2-d array shape (data_size,dimension) usually data_size >> dimension.
    cl : array_like, optional
        Cluster labels (or colors), we make up colors and glyphs.

    USAGE EXAMPLE:
    >>> import numpy as np
    >>> cl_in = np.random.randint(4,size = (2000))
    >>> bit = np.random.randint(2,size = (2000))
    >>> V = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[1,1,1,0]])
    >>> ANGLE = np.random.random(2000) * 2 * np.pi
    >>> X = np.cos(ANGLE)*0.2
    >>> Y = np.sin(ANGLE)*0.2
    >>> N = np.random.normal(size = (2000,4))*0.03
    >>> D = V[cl_in] + np.array([X,Y,np.zeros(2000),bit]).T + N
    >>> cl_out = mojave(D,cl_in)
    """
    window_name_bytes = name.encode()
    X0 = np.array(X)
    delta = np.max(X0,0) - np.min(X0,0)
    delta[delta == 0] = 1
    X1 = 2.0 * ((X0 - np.min(X0,0)) / delta) - 1.0
    X1[:,delta == 0] = 0
    Xa = np.require(X1, dtype = 'float64')
    Xp = Xa.ctypes._as_parameter_
    if cl is None:
        cl = np.zeros(len(X))
    cl_a = np.require(cl, dtype = 'int32').copy()
    cl_p = cl_a.ctypes._as_parameter_
    _mojave.mojave(Xp, cl_p, Xa.shape[0], Xa.shape[1], window_name_bytes)
    return cl_a


