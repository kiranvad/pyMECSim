# import required python packages
import numpy as np

def ReadMECSimOut(filename):
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
    if ax is None:
        fig = plt.figure(figsize = (4,4))
        ax = fig.add_subplot(111)
    ax.plot(V,I)

    ax.set_xlabel('Voltage(V)')
    ax.set_ylabel('Current(mA)')
    
    plt.tight_layout()
    sns.despine()
    
    return ax