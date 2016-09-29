import os
import tempfile
import speech

from imagesc import imagesc
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

def wview(x):
    if type(x) == str:
        os.system('wview '+x)
    else:
        speech.pcm(x).show();
    
