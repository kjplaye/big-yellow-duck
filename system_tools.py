import sys,os
import tempfile
from subprocess import Popen,PIPE

def emacs(x):
    f = tempfile.NamedTemporaryFile()
    f.write(x)
    f.flush()
    os.system('emacs ' + f.name)
    sys.stdin.readline()
    f.close()
