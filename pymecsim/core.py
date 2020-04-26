import subprocess, shlex
from pymecsim.utils import ReadMECSimOut, plotcv

def MECSIM():
    args = shlex.split('chmod u+x ../src/MECSim')
    process = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    args = shlex.split('../src/MECSim')
    process = subprocess.Popen(args,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    T,V,I = ReadMECSimOut('./MECSimOutput_Pot.txt')
    
    output = {'current':I,'voltage':V,'time':T,'info':[stdout,stderr]}
    
    return output