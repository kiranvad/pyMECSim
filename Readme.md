Python Cyclic Voltammetry Simulation Software
=======================================
This is a python wrapper for [MECSim](http://www.garethkennedy.net/MECSim.html) software.
It works completely in python in a Linux environment. I wrote this while working on [GPCV](https://github.com/kiranvad/gpcv) related work. 

If you use this software in your work please cite the original MECSim software along with this repository:
```
@misc{pymecsim,
    author       = {Kiran Vaddi},
    title        = {{pyMECSim: A Python wrapper for MECSim}},
    month        = April,
    year         = 2020,
    version      = {1.0},
    publisher    = {github},
    url          = {https://github.com/kiranvad/pyMECSim}
    }
```

This repository is arranged as follows:
1. src  : Contains the original MECSim software distributed under the same license as MECSim
2. notebooks :  An example usage of pyMECSim is shown.
3. mechanisms : Folder where you can host a original mechanism as an input file 

A sample usage is as follows:

Import `pymecsim` using the following: 
```python
from src.pymecsim import MECSIM, pysed, plotcv
```
We can perform a simulation on a one electron transfer mechanism and visualize the effect of changing the formal potential using the following code:

```python
import matplotlib.pyplot as plt
import numpy as np
import os

configfile  = '../mechanisms/cvexamples.sk'
E0 = [-0.25,0.0,0.25]
fig = plt.figure(figsize = (4,4))
ax = fig.add_subplot(111)
dirname = os.getcwd()
for i,e0 in enumerate(E0):
    outfile = dirname + '/outfile.sk'
    pysed('$E0', str(e0), configfile, outfile)
    out = MECSIM(outfile)
    ax = plotcv(out['current'],out['voltage'], ax = ax)
    ax.set_label("E0 = "+str(e0))
plt.legend([r'$E_0=0.5$',r'$E_0=0.1$',r'$E_0=1e-2$'],loc='lower right')
#plt.savefig('cvexample.png',dpi=500,bbox_inches='tight')
plt.show()
```

This will plot the following:
<img src="notebooks/cvexample.png" width="600">


Naturally, you would want to be able to run simulations for different mechanisms and confgurations. You can do that by just defining a mechanism that MECSim can model(see examples for some possible reaction mechanisms [here](http://www.garethkennedy.net/MECSimScripts.html)).

Once you have the mechanism file in say `/path/to/folder/mechanism.sk` format, turn it in as an input to pyMECSim using the following:
```python
out = MECSIM('/path/to/folder/mechanism.sk')
plotcv(out['current'],out['voltage'])
```









