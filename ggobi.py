import tempfile
import os
from subprocess import Popen, PIPE
import xml.etree.cElementTree
from collections import defaultdict

def ggobi(X):
    head = ','.join(['x_%d' % i for i in range(len(X[0]))])                    
    body = head + '\n' + '\n'.join([','.join([str(ee) for ee in e]) for e in X])    
    f = tempfile.NamedTemporaryFile(suffix = '.csv')
    f.write(body)
    f.flush()
    print('Save brushed output as %s.xml' % f.name[:-4])
    p = Popen(['ggobi',f.name])
    p.communicate()
    xml_filename = f.name[:-4] + '.xml'
    if os.path.isfile(xml_filename):
        xml_tree = xml.etree.cElementTree.parse(xml_filename)
        os.unlink(xml_filename)
        root = xml_tree.getroot()
        R = root.getchildren()[0].getchildren()[2].getchildren()
        S = [str([e for e in r.items() if e[0] == 'color' or e[0] == 'glyph']) for r in R]
        label = dict()
        num_labels = 0
        for s in S:
            if s not in label:
                label[s] = num_labels
                num_labels+=1
        return [label[s] for s in S]
