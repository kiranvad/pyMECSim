"""
A Python equivalent of sed from bash:
source : https://github.com/mahmoudadel2/pysed/blob/master/pysed.py
"""

import re
import pdb

def replace(oldstr, newstr, infile, dryrun=False):
    '''
    Sed-like Replace function..
    Usage: pysed.replace(<Old string>, <Replacement String>, <Text File>)
    Example: pysed.replace('xyz', 'XYZ', '/path/to/file.txt')
    Example 'DRYRUN': pysed.replace('xyz', 'XYZ', '/path/to/file.txt', dryrun=True) #This will dump the output to STDOUT instead of changing the input file.
    '''
    outfile = './Master.inp'
    linelist = []
    with open(infile) as f:
        for item in f:
            newitem = item.replace(oldstr, newstr)
            linelist.append(newitem)
    if dryrun == False:
        with open(outfile, "w") as f:
            f.truncate()
            for line in linelist: f.writelines(line)
            f.writelines('\n')
    elif dryrun == True:
        for line in linelist: print(line, end='')
    else:
        exit("Unknown option specified to 'dryrun' argument, Usage: dryrun=<True|False>.")