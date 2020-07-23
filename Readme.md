Python Cyclic Voltammetry Simulation Software
=======================================
This is a python wrapper for [MECSim](http://www.garethkennedy.net/MECSim.html) software.
It works completely in python in a Linux environment. I wrote this while working on [GPCV](https://github.com/kiranvad/gpcv) related work. 

If you use this software in your work please cite the original MECSim software along with this repository:
```bibtex
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
To install as a package, run 
```bash
pip install git+https://github.com/kiranvad/pyMECSim#egg=pyMECSIM.` 
```
Dependencies will be checked and installed from the setup.py file.

A sample usage is as follows:

Import `pymecsim` using the following: 
```python
from pymecsim import MECSIM, pysed
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
    model = MECSIM(outfile)
    ax = model.plot(ax = ax)
    ax.set_label("E0 = "+str(e0))
plt.legend([r'$E_0=0.5$',r'$E_0=0.1$',r'$E_0=1e-2$'],loc='lower right')
#plt.savefig('cvexample.png',dpi=500,bbox_inches='tight')
plt.show()
```

This will plot the following:
<img src="notebooks/cvexample.png" width="400">


Naturally, you would want to be able to run simulations for different mechanisms and confgurations. You can do that by just defining a mechanism that MECSim can model(see examples for some possible reaction mechanisms [here](http://www.garethkennedy.net/MECSimScripts.html)).

Once you have the mechanism file in say `/path/to/folder/mechanism.sk` format, turn it in as an input to pyMECSim using the following:
```python
model = MECSIM('/path/to/folder/mechanism.sk')
model.plot()
plt.show()
```

One can also get concentration profiles by first indicating `MECSIM` to return concentration profiles in the configuration file by setting `1		! show debug output files as well as MECSimOutput.txt (1=yes; 0=no)`. `pymecsim` will then be able to return concentration profiles as numpy arrays. see `notebooks/Cyclic Voltammetry Simulation Example for Single Electron Transfer Mechanism.ipynb` for an example use case.


## Notes
Please free to contribute to this repository both interms of code and documetation or simple example use cases in jupyter notebook. Submit a pull request and I would be happy to integrate into this repository.








