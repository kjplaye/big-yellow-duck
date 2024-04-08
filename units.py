import pint
pt = pint.UnitRegistry()

def bu(value):
    """Gets the named constant.
    
    Usage Example:
        >>> from units import pt, bu
        >>> x = 123 * pt.m / pt.s
        >>> print(bu(x))
    """
    return (1e0 * value).to_base_units()
