class col(list):
    def __init__(self,data,the_type = None):
        if isinstance(data,str):
            if the_type:
                list.__init__(self,[[the_type(e) for e in l.split()] for l in open(data).readlines()])
            else:
                list.__init__(self,[l.split() for l in open(data).readlines()])
        else:
            list.__init__(self,data)
    def transpose(self):
        m = len(self)
        n = len(self[0])
        return col([[self[i][j] for i in range(m)] for j in range(n)])
    def __repr__(self):
        L = [map(str,l) for l in self]
        FL = [1] * len(self[0])
        for l in L:
            for i in range(len(FL)):                
                if len(l[i]) > FL[i]:
                    FL[i] = len(l[i])
        return '\n'.join([' '.join([("%" + str(FL[i]) + "s") % l[i] for i in range(len(FL))]) for l in L])

def transpose(L):
    return col(L).transpose()
