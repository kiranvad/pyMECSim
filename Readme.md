Python Cyclic Voltammetry Simulation Software
=======================================
This is a python wrapper for [MECSim](http://www.garethkennedy.net/MECSim.html) software.
It works completely in python in a Linux environment. I wrote while working on [GPCV](https://github.com/kiranvad/gpcv) related work and software. 

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
1. pymecsim :  Contains python wrapper files for MECSim
2. src  : Contains the original MECSim software distributed under the same license as MECSim
3. notebooks :  An example usage of pyMECSim is shown.
4. mechanisms : Folder where you can host a original mechanism as an input file 

A sample usage is as follows:

Import two main modules of `pyMECSim` using the following: 
```python
from pymecsim.core import MECSIM
from pymecsim import pysed
```
We can perform a simulation on a one electron transfer mechanism and visualize the effect of changing the formal potential using the following code:

```python
import pdb
from pymecsim.utils import plotcv
import matplotlib.pyplot as plt
import numpy as np

mechanism  = '../mechanisms/cvexamples.sk'
E0 = [-0.25,0.0,0.25]
fig = plt.figure(figsize = (4,4))
ax = fig.add_subplot(111)
for i,e0 in enumerate(E0):
    pysed.replace('$E0', str(e0), mechanism)
    out = MECSIM()
    T = out['time']
    forward_sweep = np.arange(len(T)/2,len(T)).astype(int)
    ax = plotcv(out['current'],out['voltage'], ax = ax)
    ax.set_label("E0 = "+str(e0))
plt.legend([r'$E_0=0.5$',r'$E_0=0.1$',r'$E_0=1e-2$'],loc='lower right')
plt.show()
```

This will plot the following:
<img src="noteboos/cvexample.png" width="600">
