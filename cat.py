import hashlib
import numpy as np
import pandas as pd


class Cat(np.ndarray):
    """Expose some pandas functionality to consistently index categorical data.
    Basically an augmented factorize.

    Usage Example:
    >>> from cat import Cat, KeyCat
    >>> key_cat = KeyCat(['c','b','a'])
    >>> cat = Cat(['a','c','b','a'])
    >>> new_cat = cat(['c','a'])
    >>> print("Codes:",cat)
    >>> print("Number of categories:", cat.ncat)
    >>> print("Uniques:", cat.cat)
    >>> print("Pull in new data:", new_cat)
    >>> print("Lift to names:", new_cat.lift())
    >>> print("Hash:", hash(cat))
    >>> print("Key use 1:", key_cat('a'))
    >>> print("Key use 2:", key_cat(['a','b','c']))
    Codes: [0 2 1 0]
    Number of categories: 3
    Uniques: Index(['a', 'b', 'c'], dtype='object')
    Pull in new data: [2 0]
    Lift to names: Index(['c', 'a'], dtype='object')
    Hash: 971233737139621537
    Key use 1: 2
    Key use 2: [2 1 0]
    """
    def __new__(self, data, categories=None):
    	categorical = pd.Categorical(data, categories=categories)
    	series = pd.Series(categorical)
    	codes = series.cat.codes
    	obj = np.asarray(codes.values).view(self)
    	obj.series = series
    	obj.cat = categorical.categories
    	obj.ncat = obj.cat.nunique()
    	return obj
    
    def __call__(self, data):
    	"""Apply Cat indexing to new data."""
    	if np.ndim(data)!=0:
        	return Cat(data, categories=self.cat)
    	else:
        	return Cat([data], categories=self.cat)[0]

    def lift(self):
    	return self.cat[self]

    def __hash__(self):
    	return int(hashlib.sha256(str((tuple(self), tuple(self.cat)))
                              	.encode()).hexdigest(), 16)

    def __repr__(self):
    	return str(self) + '\n' + str(self.cat)

def KeyCat(x):
    return Cat([],x)
