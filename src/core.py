import subprocess, shlex
import pdb, os, sys
import numpy as np
import pandas as pd
import os
dirname = os.path.dirname(__file__)

from shutil import copyfile
import warnings

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

class MECSIM:
    def __init__(self, configfile):
        """
        Simulate a voltammetry response for a given configuration.
        
        Input:
        -----
            configfile : An MECSIM configuration file usually in the .inp form
        
        Attributes:
        -----------
            T  : Time profile
            I  : Current profile
            V  : Voltage profile
                    
        Note
        -----
            This class is a python wrapper for the MECSIM software by Gareth F Kennerdy (http://www.garethkennedy.net/MECSim.html).
            It simulates a mechanism with the configurations specified and throws errors or warnings based on the results of simulations.
            While the python wrapper provides access to many functions of MECSIM, it sometimes fails to understand the errors from MECSIM.
            Best practice is to always look at the log.txt file generated in your working directory and figure out what has gone wrong.
        """
        self.configfile = configfile
        self.dirname = os.path.dirname(__file__)
        self.solve()
        self._get_errors()
        self.T,self.V,self.I = self._read_mecsim_out(self._get_filename('./'+self.outfile))
        
    def solve(self):
        main_configfile = os.path.join(self.dirname, 'Master.inp')
        copyfile(self.configfile, main_configfile)

        args = shlex.split('chmod u+x ./MECSim')
        process = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=dirname)
        args = shlex.split('./MECSim')
        process = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=dirname)
        self.stdout, self.stderr = process.communicate()

        logfile = os.path.join(os.getcwd(),'log.txt')
        copyfile(os.path.join(dirname, 'log.txt'), logfile)
        file_path = self._get_filename('log.txt')
        f = open(file_path, 'r')
        self.logs = list(f)
        self.flag_concs = False
        for line in reversed(self.logs):
            if 'Output file' in line:
                self.outfile = line.split(' ')[-1].split('\n')[0]

            elif 'Additional files' in line:
                self.flag_concs= True
                
    def _get_errors(self):
        for line in reversed(self.logs):
            if "Error" in line:
                error = line.replace('Error: ','')
                break

        if not len(self.stderr)==0:
            raise RuntimeError('MECSIM threw the following error that pyMECSIM could not understand \n'+
                               self.stderr.decode("utf-8"))

        if self.stderr.decode("utf-8") == 'STOP 4\n':
            raise RuntimeError('pyMECSIM could not understand the error message from MECSIM.')

        if len(self.logs)<30:
            if 'error' in locals():
                raise RuntimeError(error)
            else:
                raise RuntimeError('Could not identify a clear error from MECSIM but here is what is avaiable: \n'
                                   + logs[-1])
        else:
            if 'error' in locals():
                warnings.warn(error)
                
    def get_current_profile(self):
        """ return current values over simulated time scale """
        return self.I
    def get_voltage_profile(self):
        """ 
        return voltage applied over simulated time scale. 
        Read MECSIM paper for more details on what is applied voltage.
        https://doi.org/10.1016/j.coelec.2016.12.001
        
        """
        return self.V
    
    def get_time_profile(self):
        """
        return time profile of the simulation.
        """
        return self.T

    def _read_mecsim_out(self, filename):
        """ internal function to read the text output """
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
    
    def get_bulk_concentrations(self):
        """
        returns bulk concentrations C(x,t=t_end).
        
        output is a numpy array of shape (N,P) where N is number of grid points in x-direction and P is number of species
        """
        if not self.flag_concs:
            raise KeyError('MECSIM did not output any concentration profiles.'+
                           'Did you specifiy in the configuration file to output the concentration profiles?')
        
        fin_file = self._get_filename('EC_Model.fin')
        f = open(fin_file, 'r')
        concs = []
        for i, line in enumerate(f):
            columns = line.split()
            concs.append(columns)
        concs = np.asarray(concs) 
        concs = self._get_float_matrix(concs)

        return concs
    
    def get_surface_concentrations(self):
        """
        returns surface concentrations profiles over time C(x=0,t).
        
        output is a numpy array of shape (N,P) where N is number of time points and P is number of species
        """
        if not self.flag_concs:
            raise KeyError('MECSIM did not output any concentration profiles.'+
                       'Did you specifiy in the configuration file to output the concentration profiles?')
        tcv_file = self._get_filename('EC_Model.tvc')    
        f = open(tcv_file, 'r')
        concs = []
        for i, line in enumerate(f):
            if i>0:
                columns = line.split()
                concs.append(columns[3:])
        concs = np.asarray(concs) 
        concs = self._get_float_matrix(concs)   

        return concs
    
    def plot(self,ax = None):
        """ simple visualization tool for CV curve """
        if ax is None:
            fig = plt.figure(figsize = (4,4))
            ax = fig.add_subplot(111)
        ax.plot(self.V,self.I)

        ax.set_xlabel('Voltage(V)')
        ax.set_ylabel('Current(A)')

        plt.tight_layout()
        sns.despine()

        return ax
    
    def _get_filename(self,filename):
        file_path = os.path.join(dirname, './', filename)

        if os.path.exists(file_path):
            return file_path
        else:
            return False

    def _get_float_matrix(self, matrix):
        """ 
        Given a numpy matrix of floats, returns a float matrix 
        This is useful because sometimes you get a string that can not be made to a float
        """
        def to_float(string):
            try:
                f = float(string)
            except ValueError:
                f = string.replace("-", 'E-')

            return f

        df = pd.DataFrame(matrix)
        df = df.applymap(to_float)

        return df.to_numpy()