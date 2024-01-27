import os
import tempfile
import ctypes
import numpy as np

from imagesc import imagesc
from wview import wview
from ggobi import *

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
    
