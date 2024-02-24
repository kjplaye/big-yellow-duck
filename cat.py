import numpy as np
import pandas as pd

class Cat(np.ndarray):
    """Expose some pandas functionality to consistently index categorical data.
    
    Usage Example:
    >>> cat = Cat(['a','c','b','a'])
    >>> print(cat)
    [0, 2, 1, 0]

    >>> print("Number of categories:", cat.num_cat)
    Number of categories: 3

    >>> print("Categories:", cat.cat)
    Categories: Index(['a', 'b', 'c'], dtype='object')

    >>> new_cat = cat.index(['c','a'])
    >>> print("Pull in new data:", new_cat))
    >>> print("As a series:")
    >>> print(new_cat.series)
    >>> print("Lift to names:", new_cat.lift())
    Pull in new data: [2 0]
    As a series:
    0   c
    1   a
    dtype: category
    Categories (3, object): ['a', 'b', 'c']
    Lift to names: Index(['c', 'a'], dtype='object')
    """
    def __new__(self, data, categories = None):
        cat = pd.Categorical(data, categories = categories)
        series = pd.Series(cat)
        codes = series.cat.codes
        obj = np.asarray(codes.values).view(self)
        obj.series = series
        obj.cat = cat.categories
        obj.num_cat = obj.cat.nunique()
        return obj
    def index(self, data):
        """Apply Cat indexing to new data."""
        return Cat(data, categories=self.cat)
    def lift(self):
        return self.cat[self]
    def __hash__(self):
        return hash((tuple(self), tuple(self.cat)))
