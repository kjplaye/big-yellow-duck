import multiprocessing as mp
from tqdm import tqdm, trange

def tmap(p, f, L):
    return list(tqdm(p.imap(f, L), total=len(L)))

def _mp_mapper(pair):
    """Multiprocessing helper mapper, mostly used in mp_dr."""
    f0, x = pair
    f = dill.loads(f0) if isinstance(f0, str) else f0
    return f(*x)    


def mp_dr(f, inputs):
    """Multiprocessing dereference.  Mainly useful for functions with multiple inputs.
    Entries in inputs that are tuples will derefernce to f arguments, other types will
    act as single inputs.  
    
    USAGE_EXAMPLE:
    >>> def f(x,y):
    >>> 	return x+y
    >>> inputs = [(1,2),(4,5),(41,25)]
    >>> out = map(*mp_dr(f, inputs))
    >>> print(list(out))
    """
    function_input_pairs = [[f, inp if isinstance(inp, tuple) else (inp,)] for inp in inputs]
    return (_mp_mapper, function_input_pairs)


def vmap(f, inputs, num_cores = 1, verbose = True):
    """Maps f on inputs with various numbers of cores and verbosity.
    
    USAGE_EXAMPLES:
    >>> import time
    >>> def f(x): time.sleep(x); return x
    >>> print(sum(vmap(f, [0.1] * 20, num_cores = 1, verbose = False)))
    >>> print(sum(vmap(f, [0.1] * 20, num_cores = 1, verbose = True)))
    >>> print(sum(vmap(f, [0.1] * 20, num_cores = 2, verbose = False)))
    >>> print(sum(vmap(f, [0.1] * 20, num_cores = 2, verbose = True)))
    """
    MPDR = mp_dr(f, inputs)
    if num_cores == 1:
    if verbose:
    	return list(tqdm(map(*MPDR), total = len(inputs)))
    else:
    	return list(map(*MPDR))
    else:
    p = mp.Pool(num_cores)
    if verbose:
    	return list(tmap(p,*MPDR))
    else:
    	return p.map(*MPDR)
