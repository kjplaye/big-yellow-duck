import hashlib
import numpy as np
import pandas as pd
from typing import Callable, List, Optional

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
    >>> print("Stripped to KeyCat:", list(Cat(['a','a','b'], strip = True)))
    >>> print("Using sort key:", Cat(['a','b','A','B','a','a'], key = str.casefold).cat)
    >>> print("Pickling works", pickle.loads(pickle.dumps(cat)))
    Codes: [0 2 1 0]
    Number of categories: 3
    Uniques: Index(['a', 'b', 'c'], dtype='object')
    Pull in new data: [2 0]
    Lift to names: Index(['c', 'a'], dtype='object')
    Hash: 971233737139621537
    Key use 1: 2
    Key use 2: [2 1 0]
    Stripped to KeyCat: []
    Using sort key: Index(['a', 'A', 'B', 'b'], dtype='object')
    Pickling works [0 2 1 0]
    """
    series: pd.Series
    cat: pd.Series
    ncat: int

    def __new__(self, data: List, categories: Optional[List] = None, strip = False,
            	key: Optional[Callable] = None) -> "Cat":
    	""" Create a new Cat.

    	Args:
            data: List of things to index.
            categories: Indexing set provided.
            strip: When set it strips the data creating essentially a KeyCat.
            key: When set, it uses sort_key as a key to sort the categories.

    	Returns:
        	A new Cat.
    	"""
    	if key is not None:
            sorted_categories = sorted(set(data), key = key)
            return Cat(data, sorted_categories, strip, key = None)
    	categorical: pd.Categorical = pd.Categorical(data, categories=categories)
    	series: pd.Series = pd.Series(categorical)
    	codes: pd.Series[int] = series.cat.codes[:0] if strip else series.cat.codes
    	obj: Cat = np.asarray(codes.values).view(self)
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

    def __reduce__(self):
    	# Get the parent's __reduce__ tuple
    	pickled_state = super(Cat, self).__reduce__()
    	# Create our own tuple to pass to __setstate__
    	new_state = pickled_state[2] + (self.cat, self.ncat)
    	# Return a tuple that replaces the parent's __setstate__ tuple with our own
    	return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
    	self.cat, self.ncat = state[-2:]  # Set the info attribute
    	# Call the parent's __setstate__ with the other tuple elements.
    	super(Cat, self).__setstate__(state[0:-2])
    	

def KeyCat(x):
    return Cat([],x)
