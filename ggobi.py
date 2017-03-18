import tempfile
import os
from subprocess import Popen, PIPE
import xml.etree.cElementTree
from collections import defaultdict
import numpy as np
import six

if six.PY3:
    py3 = True
else:
    py3 = False

def ggobi(X, cl = None, smallest_glyph = 'auto'):
    """
    Simple python wrapper for GGobi.

    Parameters
    ----------
    X : array_like
        2-d array shape (data_size,dimension) usually data_size >> dimension.
    cl : array_like, optional
        Cluster labels (or colors), we make up colors and glyphs.
    smallest_glyph : string, optional
        'fast' - use smallest glyph (a pixil) so that drawing is fast
        'visible' - use visible glyphs, but it slows down for larger data sizes.
        'auto' - use 'fast' if we have 1000 data points or more.

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
    >>> cl_out = ggobi(D,cl_in)
    """
    if cl == None:
        cl = [0] * len(X)
    # Use colors, shapes, and sizes in a hopefully diverse and useful order
    # We choose "." as default shape at smallest size for speed on large data
    if smallest_glyph == 'auto' and len(X) >= 1000 or smallest_glyph == 'fast':
        use_pixil_glyph = True
    else:
        use_pixil_glyph = False
    color = ['0','1','2','3','4','5','6','7','8']
    glyph_shape = ['fc','plus','or','x','oc','fr']
    glyph_size = ['1','4','7','2','0','3','5','6']
    glyph = ['%s %s' % (shape,size) for size in glyph_size for shape in glyph_shape]
    if use_pixil_glyph:
        glyph = ['. 1'] + glyph
    label = ['color="%s" glyph="%s"' % (c,g)  for g in glyph for c in color]
    CL = [label[c] for c in cl]
    m = len(X)
    n = len(X[0])
    x = ''
    x += '<?xml version="1.0"?>\n'
    x += '<!DOCTYPE ggobidata SYSTEM "ggobi.dtd">\n'
    x += '\n'
    x += '<ggobidata count="1">\n'
    x += '<data name="fake_name.csv">\n'
    x += '<description>\n'
    x += 'This is XML created by GGobi\n'
    x += '</description>\n'
    x += '<variables count="%d">\n' % n
    for i in range(n):
        x += '    <realvariable name="x_%d" nickname="x_"/>\n' % i
    x += '</variables>\n'
    x += '<records count="%d" glyph="fc 1" color="0">\n' % m
    for i in range(m):
        x += '<record label="%d" %s>\n' % (i,CL[i])
        x += ' '.join(['<real>%f</real>' % X[i][j] for j in range(n)]) + '\n'
        x += '</record>\n'
    x += '</records>\n'
    x += '</data>\n'
    x += '</ggobidata>\n'
    f = tempfile.NamedTemporaryFile(suffix = '.xml')
    if py3:
        x = bytes(x,encoding='ascii')
    f.write(x)
    f.flush()
    print('Save brushed output as %s.xml' % f.name[:-4])
    print('Turn "cycle" on and off to unlock menu')
    p = Popen(['ggobi',f.name])
    p.communicate()
    xml_filename = f.name[:-4] + '.xml'
    if os.path.isfile(xml_filename):
        xml_tree = xml.etree.cElementTree.parse(xml_filename)
        root = xml_tree.getroot()
        R = root.getchildren()[0].getchildren()[2].getchildren()
        S = [str([e for e in r.items() if e[0] == 'color' or e[0] == 'glyph']) for r in R]
        label = dict()
        num_labels = 0
        for s in S:
            if s not in label:
                label[s] = num_labels
                num_labels+=1
        return np.array([label[s] for s in S])
