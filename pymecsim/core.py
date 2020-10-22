import subprocess, shlex
import pdb, os, sys
import numpy as np
import pandas as pd
import os
dirname = os.path.dirname(__file__)

from shutil import copyfile
import warnings

import matplotlib.pyplot as plt

from .utils import pysed

class MECSIM:
    def __init__(self, configfile=None, exp=None):
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
        self.dirname = os.path.dirname(__file__)
        self.exp = exp
        self.configfile = configfile
        self.outfile = None

        if self.configfile is None:
            inpfile = exp.get_inpfile_lines()
            self.configfile = os.path.join(self.dirname, 'from_expt.inp')
            with open(self.configfile, 'w') as f:
                for item in inpfile:
                    f.write(item)
        elif self.exp is None:            
            self.configfile = configfile
        else:
            raise RuntimeError('At least one input required. Use either a pymecsim exp class or a MECSIM .inp file')
        

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
        self._get_errors()
        self.T,self.V,self.I = self._read_mecsim_out(self._get_filename('./'+self.outfile))  
        
        if self.flag_concs:
            _file = self._get_filename('EC_Model.fin')
            self.fin_file = list(open(_file, 'r'))
            _file = self._get_filename('EC_Model.tvc')
            self.tvc_file = list(open(_file, 'r'))
        
        return self.T, self.V , self.I
    
    def _get_errors(self):
        
        converged = False
        # try if MECSim has thrown any fortran errors
        if not len(self.stderr)==0:
            raise RuntimeError('Unknown error but could be that some ' 
                                   'of the simulation parameters are not properly set')

        
        # try to get meaninful errors from log.txt
        errors = []
        for line in reversed(self.logs):
            if "Error" in line:
                errors.append(line.replace('Error: ',''))
            if "unphysical" in line:
                errors.append('Simulation diverged: '+ line.strip('with \n'))
            if "Done 100%" in line:
                converged = True
                
        error = '\n'.join(e for e in errors)        
        if not converged:
            warnings.warn('Simulation diverged\nLook at log.txt file in the working directory')
            if len(errors)>0:
                raise RuntimeError(error)
            else:
                raise RuntimeError('Could not identify a clear error from MECSIM but here is what is avaiable: \n'
                                   + logs[-1])

        
        # if nothing from log files makes sense:
        if self.outfile is None:
            raise RuntimeError('Simulation exited without any errors. Please look at log.txt file')
                
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
        
        output is dictonaty with keys as species names. Access the normalized length using key 'x'
        """
        if not self.flag_concs:
            raise Exception('MECSIM did not output any concentration profiles.'+
                           'Did you specifiy in the configuration file to output the concentration profiles?')
            
        if len(self.fin_file)==0:
            warnings.warn('No bulk concentration profiles found')
            
        concs = []
        for i, line in enumerate(self.fin_file):
            columns = line.split()
            concs.append(columns)
        concs = np.asarray(concs) 
        concs = self._get_float_matrix(concs)
        out = {}
        out['x'] = concs[1:,0]
        for i,s in enumerate(self.exp.mechanism.species):
            out[s.name] = concs[1:,i+1]
            
        return out
    
    def get_surface_concentrations(self):
        """
        returns surface concentrations profiles over time C(x=0,t).
        
        output is dictonaty with keys as species names. Access the time array over t using key 'T'
        
        Notes:
        ------
                EC_Model.tvc is structured such that first line is number of species
                next few lines can be trerated as a table with columns ordered as time, voltage,unkown, concentrations
        """
        if not self.flag_concs:
            raise Exception('MECSIM did not output any concentration profiles.'+
                       'Did you specifiy in the configuration file to output the concentration profiles?')
        
        if len(self.tvc_file)==0:
            warnings.warn('No surface concentration profiles found')
        
        table = []
        for i, line in enumerate(self.tvc_file):
            if i>0:
                columns = line.split()
                table.append(columns)
        table = np.asarray(table) 
        table = self._get_float_matrix(table)
        
        # only collect rows with a positive time stamp
        positive_time = table[:,0]>0
        table = table[positive_time,:]
        
        out = {}
        out['T'] = table[:,0]
        concs = table[:,3:]
        for i, s in enumerate(self.exp.mechanism.species):
            out[s.name] = concs[:,i]
            
        return out
    
    def plot(self,ax = None):
        """ simple visualization tool for CV curve """
        if ax is None:
            fig = plt.figure(figsize = (4,4))
            ax = fig.add_subplot(111)
        ax.plot(self.V,self.I)

        ax.set_xlabel('Voltage(V)')
        ax.set_ylabel('Current(A)')

        plt.tight_layout()

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
    
    def set_parameter(self,string, value):
        """
        Utility function to set a parameter value when using Experiment class as input
        """
        value = '{:.2E}'.format(value).replace('E','e')
        pysed(string,value, self.configfile, self.configfile)
        
        
    def plot_concentrations(self, **kwargs):
        """Utility function to plot concentrations
        
        You can pass any kwargs of the class `matplotlib.pyplot.subplots`
        Returns axis for two subplots [0] Bulk concentrations, [1] Sufface concentrations
        """
        
        species = self.exp.mechanism.species
        C_bulk = self.get_bulk_concentrations()

        fig, axs = plt.subplots(1,2, **kwargs)
        fig.subplots_adjust(wspace=0.3)

        for s in species:
            axs[0].plot(C_bulk['x'], C_bulk[s.name], label=s.name)
            
        axs[0].set_xlabel('Normalized length')
        axs[0].set_ylabel('Normalized concentration')
        axs[0].legend()
        axs[0].set_title(r'Species concentration in the bulk (Final t=$\infty$)')

        C_surface = self.get_surface_concentrations()
        time = C_surface['T']

        for s in species:
            axs[1].plot(time, C_surface[s.name], label=s.name)

        axs[1].set_xlabel('time (s)')
        axs[1].set_ylabel('Normalized concentration')
        axs[1].legend()
        axs[1].set_title('Species concentration at the surface')
        
        return axs
        
        
        
        
        