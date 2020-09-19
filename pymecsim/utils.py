"""Change logs

Version 2:
----------
1. Add depreceated function as a utility
"""

import re
import pdb

def pysed(oldstr, newstr, infile, outfile):
    """
    A utility function to modify a skeleton configuration file iteratively.
    See examples for a sample usage.
    
    Inputs:
    --------
        oldstr  : string uses as a variable in .sk file (usually starts with a $ sign)
        newstr  : string used to replace the `oldstr` variable (although its  float one should pass a float)
        infile  : file with variables defined
        outfile : file with variables replaced by the floats we want.
    
    
    --------------------------------------------------------------------
    Sed-like Replace function..
    Usage: pysed.replace(<Old string>, <Replacement String>, <Text File>)
    Example: pysed.replace('xyz', 'XYZ', '/path/to/file.txt')
    Example 'DRYRUN': pysed.replace('xyz', 'XYZ', '/path/to/file.txt', dryrun=True) 
    #This will dump the output to STDOUT instead of changing the input file.

    A Python equivalent of sed from bash:
    source : https://github.com/mahmoudadel2/pysed/blob/master/pysed.py
"""

    linelist = []
    with open(infile) as f:
        for line in f:
            items = line.split('!')[0].split(',')
            if oldstr in [x.replace(" ","") for x in items]:
                newline = line.replace(oldstr, newstr)
                linelist.append(newline)
            else:
                linelist.append(line)
                    
    with open(outfile, "w") as f:
        f.truncate()
        for line in linelist: f.writelines(line)
        f.writelines('\n')
        
        
import warnings
import functools
# from https://stackoverflow.com/a/30253848
def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning)  # turn off filter
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning,
                      stacklevel=2)
        warnings.simplefilter('default', DeprecationWarning)  # reset filter
        return func(*args, **kwargs)
    return new_func        
        
        