import subprocess, shlex
import pdb
import os
import sys
if '../' not in sys.path:
    sys.path.insert(0,'../')
    
import os
dirname = os.path.dirname(__file__)

from shutil import copyfile

def MECSIM(configfile='Master.sk'):
    """
    This is a python wrapper for the original MECSIM software re-distrubuted along with this packages as a .exe file
    It takes an MECSim configuration file in the configfile as input and returns a dictonarry.
    
    INPUT
    -----
        configfile  :  an MECSim configuration file (see ../mechanisms/cvexamples.sk for an example)
                       For advanced usage, refer to the original MECSim user manual
    
    OUTPUT
    ------
        output  : a dictonary with the 'current', 'voltage', 'time' and 'info':output from MECSim software
        
    
    NOTES:
    ------
        A log.txt file will be created in the `src` folder with detailed output from MECSim
        (in future this willbe moved to be as an output file)
    
    """
    main_configfile = os.path.join(dirname, 'Master.inp')
    copyfile(configfile, main_configfile)
    
    args = shlex.split('chmod u+x ./MECSim')
    process = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=dirname)
    args = shlex.split('./MECSim')
    process = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=dirname)
    stdout, stderr = process.communicate()
    
    T,V,I = ReadMECSimOut(os.path.join(dirname, './MECSimOutput_Pot.txt'))
    
    output = {'current':I,'voltage':V,'time':T,'info':[stdout,stderr]}
    
    return output


import re

def pysed(oldstr, newstr, infile, outfile):
    """
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
        for item in f:
            newitem = item.replace(oldstr, newstr)
            linelist.append(newitem)
    with open(outfile, "w") as f:
        f.truncate()
        for line in linelist: f.writelines(line)
        f.writelines('\n')
        
        
import numpy as np

def ReadMECSimOut(filename):
    """
    Internal function used to read output from MECSim to simple numpy arrays
    """
    f = open(filename, 'r')
    # search for last line of header that is made by MECSim (always this line)
    for line in f:
        if line.strip() == "Post(ms):       0.000000E+00": break
    time = []
    eapp = []
    current = []
    for line in f:
        columns = line.split()
        time.append(float(columns[2]))
        eapp.append(float(columns[0]))
        current.append(float(columns[1]))
        
    return np.asfarray(time), np.asfarray(eapp), np.asfarray(current)

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

def plotcv(I,V,ax = None):
    """
    A utility function to plot CV curve. I is output['current'] and V is output['voltage']
    """
    if ax is None:
        fig = plt.figure(figsize = (4,4))
        ax = fig.add_subplot(111)
    ax.plot(V,I)

    ax.set_xlabel('Voltage(V)')
    ax.set_ylabel('Current(mA)')
    
    plt.tight_layout()
    sns.despine()
    
    return ax
   