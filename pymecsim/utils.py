import re

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
            if oldstr in line:
                newline = line.replace(oldstr, newstr)
                linelist.append(newline)
            else:
                linelist.append(line)
                    
    with open(outfile, "w") as f:
        f.truncate()
        for line in linelist: f.writelines(line)
        f.writelines('\n')