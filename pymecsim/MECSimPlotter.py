
# coding: utf-8

# MECSim plotter for ``MECSimOutput_Pot.txt``
# -----
# 
# Output the current against time and against voltage as 2 sets of images. Each is output as PNG, PS and PDF in a clear format for inclusion in scientific journal articles.

# In[1]:

# import required python packages
import numpy as np
import matplotlib.pyplot as plt


# In[2]:

# load MECSim output file
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
        time.append(float(columns[0]))
        eapp.append(float(columns[1]))
        current.append(float(columns[2]))
    return np.asfarray(time), np.asfarray(eapp), np.asfarray(current)


# In[3]:

eapp, current, time = ReadMECSimOut('MECSimOutput_Pot.txt')


# In[4]:

plt.figure(figsize=(8,6),dpi=100)
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.width'] = 2
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.plot(time, eapp, 'k')#, label='A', linewidth=2)
plt.xlim([0,time[-1]])
plt.xlabel('Time (s)', fontsize=20)
plt.ylabel('Voltage (V)', fontsize=20)
plt.savefig('MECSimOutputVT.png')



# In[5]:

plt.figure(figsize=(8,6),dpi=100)
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.plot(eapp, current, 'k')#, label='A', linewidth=2)
plt.xlabel('Voltage (V)', fontsize=20)
plt.ylabel('Current (A)', fontsize=20)
plt.savefig('MECSimOutputCV.png')


