import os
import tempfile
import speech

from imagesc import imagesc
#def imagesc(x):
#    if type(x) == str:
#        os.system('imagesc '+x)
#    else:
#        f = tempfile.NamedTemporaryFile()
#        for L in x:
#            f.write(' '.join(map(str,L)) + '\n')
#        f.flush()
#        os.system('imagesc '+f.name)
#        f.close()

def plotn(x):
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
    
