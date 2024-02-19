import multiprocessing as mp
from tqdm import tqdm, trange

def tmap(p, f, L):
    return list(tqdm(p.imap(f, L), total=len(L)))

def _mp_mapper(pair):
    """Multiprocessing helper mapper, mostly used in mp_dr."""
    f, x = pair
    return f(*x)    

def mp_dr(f, inputs):
    """Multiprocessing dereference.  Mainly useful for functions with multiple inputs.

    USAGE_EXAMPLE:
    >>> def f(x,y):
    >>> 	return x+y
    >>> inputs = [(1,2),(4,5),(41,25)]
    >>> out = map(*mp_dr(f, inputs))
    >>> print(list(out))
    """
    function_input_pairs = [[f, inp] for inp in inputs]
    return (_mp_mapper, function_input_pairs)

