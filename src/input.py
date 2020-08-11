import pdb

def process_parameter(param):
    if isinstance(param, str):
        out = '{}'.format(param)
    elif isinstance(param, float):
        out = '{:.2E}'.format(param).replace('E','e')
    elif isinstance(param, int):
        out = '{:d}'.format(param)
    else:
        raise ValueError('{} type is not recognised for {}'.format(type(param)))
    
    return out

class Specie:
    """
    Usage:
    >>> A = Specie('A')
    >>> print(A)
    Solution specie with Diffusion coefficent 1.00E-05 cm2/s , intial concentration 0.00E+00 mol/cm3
    """
    
    def __init__(self, name,D=1e-5,C0 = 0.0, surface_confined=False):
        self.name = name
        self._set_specie_type(surface_confined)
        self.diffusion(D)
        self.concentration(C0)
        
    def diffusion(self, D):
        self.D = D
        self.D_units = 'cm2/s'

    def concentration(self, C0):
        self.C0 = C0  
        if self._specie_type==1:
            self.C_units = 'mol/cm2'
        else:
            self.C_units = 'mol/cm3'
            
    def _set_specie_type(self, surface_confined):
        if surface_confined:
            self.specie_type = 'Surface confined'
            self._specie_type = 1
        else:
            self.specie_type = 'Solution'
            self._specie_type = 0
    
    def __repr__(self):
        line = '{} specie with Diffusion coefficent {:.2E} {} , intial concentration {:.2E} {}'
        return line.format(self.specie_type, self.D, self.D_units, self.C0, self.C_units)
        
    def get_input(self):
        _D = process_parameter(self.D)    
        _C0 = process_parameter(self.C0)
        line = '{}, {}, {} ! {} \n'.format(_C0, _D, self._specie_type, self.name)
        
        return line
        
class Capacitance:
    """
    Usage:
    import numpy as np
    >>> cap = Capacitance(0.0, np.zeros(5)) 
    """
    def __init__(self,Epzc, coeffs):
        self.epzc = Epzc
        self.coeffs = coeffs
    def get_input(self):
        lines = ['{:d} ! Number of terms\n'.format(len(self.coeffs)-1)]
        lines.append('{} ! E_pzc (V)\n'.format(process_parameter(self.epzc)))
        for i, a in enumerate(self.coeffs):
            lines.append('{} ! a_{:d}\n'.format(process_parameter(a), i).replace('E','e')) 
            
        return lines
    
class Reaction:
    def __init__(self, reactants, products):
        self.reactants = reactants
        self.products = products
        self.num_species = len(self.reactants) + len(self.products)
    
    def __repr__(self):
        lhs = ['{} {}'.format(value, key) for key, value in self.reactants.items()]
        rhs = ['{} {}'.format(value, key) for key, value in self.products.items()]
        
        reaction =self.mode[0] +' : '
        for i, s in enumerate(lhs):
            if not i==len(lhs)-1:
                reaction += s + ' + '
            else:
                reaction += s 
                
        reaction += ' <=> '
        
        for i, s in enumerate(rhs):
            if not i==len(rhs)-1:
                reaction += s + ' + '
            else:
                reaction += s 
                
        reaction += self.params
        
        return reaction
    
class ChargeTransfer(Reaction):
    """
    Usage:
    >>> R1 = ChargeTransfer({'A':1,'e':2},{'B':1},0.0)
    >>> print(R1) 
    Charge Transfer : 1 A + 2 e <=> 1 B  ks= 1.00E+04, E0 = 0.00E+00, alpha = 0.50
    """
    def __init__(self, reactants, products, E0, ks=1e4, alpha=0.5):
        self.ks = ks
        self.E0 = E0
        self.alpha = alpha
        self.mode = ['Charge Transfer', 0]
        self.params = self.get_params_as_string()
        super().__init__(reactants, products)

     
    def get_params_as_string(self):
        params = '  ks= {:.2E}, '.format(self.ks)
        params += 'E0 = {:.2E}, '.format(self.E0)
        params += 'alpha = {:.2f}'.format(self.alpha)
        
        return params

class ChemicalReaction(Reaction):
    """
    usage:
    >>> R2 = ChemicalReaction({'A':1,'B':1},{'C':1,'D':1})
    >>> print(R2)
    Chemical Reaction : 1 A + 1 B <=> 1 C + 1 D  kf= 1.00E+04, kb= 1.00E+04  
    """
    def __init__(self, reactants, products, kf=1e4, kb=1e4):
        self.kf = kf
        self.kb = kb
        self.mode = ['Chemical Reaction', 2]
        self.params = self.get_params_as_string()
        super().__init__(reactants, products)

     
    def get_params_as_string(self):
        params = '  kf= {:.2E}, '.format(self.kf)
        params += 'kb= {:.2E} '.format(self.kb)
        
        return params    

class CatalyticReaction(Reaction):
    """
    Usage:
    >>> R3 = CatalyticReaction({'B':1},{'C':1})
    >>> print(R3)
    Catalytic Reaction : 1 B <=> 1 C  kf= 1.00E+04, kb= 1.00E+04 
    """
    def __init__(self, reactants, products, kf=1e4, kb=1e4):
        self.kf = kf
        self.kb = kb
        self.mode = ['Catalytic Reaction', 1]
        self.params = self.get_params_as_string()
        super().__init__(reactants, products)

     
    def get_params_as_string(self):
        params = '  kf= {:.2E}, '.format(self.kf)
        params += 'kb= {:.2E} '.format(self.kb)
        
        return params 
           
class Mechanism:
    """
    Usage:
    >>> A = Specie('A', C0=1e-6)
    >>> B = Specie('B')
    >>> C = Specie('C')
    >>> D = Specie('D')
    >>> species = [A, B, C, D]

    >>> R1 = ChargeTransfer({'A':1,'e':1},{'B':1},0.0)
    >>> R2 = ChemicalReaction({'B':1,'C':1},{'A':1,'D':1})
    >>> rxn = [R1, R2]
    >>> mech = Mechanism(species, rxn)
    >>> print(mech)  

    Charge Transfer : 1 A + 1 e <=> 1 B  ks= 1.00E+04, E0 = 0.00E+00, alpha = 0.50
    Chemical Reaction : 1 B + 1 C <=> 1 A + 1 D  kf= 1.00E+04, kb= 1.00E+04
    """
    def __init__(self, species, reactions):
        self.species = species
        self.reactions = reactions
        self.num_species = len(species)
        self.specie_list = [s.name for s in species]   
        self.input = self.get_input()

    def get_reaction_dict(self):
        reaction_input_dict = {}
        for s in self.specie_list:
            reaction_input_dict[s]=0
        
        return reaction_input_dict
        
    def process_reaction(self, reaction):
        if reaction.mode[1]==0:
            line = '{}, '.format(process_parameter(reaction.mode[1]))
            reaction_input_dict = self.get_reaction_dict()
            for key, coeff in reaction.reactants.items():
                reaction_input_dict[key]= '-{}'.format(process_parameter(coeff))
            for key, coeff in reaction.products.items():
                reaction_input_dict[key]='{}'.format(process_parameter(reaction.reactants['e'])) 
            del reaction_input_dict['e']
            
            line += ', '.join(process_parameter(value) for key, value in reaction_input_dict.items() )
            line += ', 0.0e0, 0.0e0, '
            line += '{}, '.format(process_parameter(reaction.E0))
            line += '{}, '.format(process_parameter(reaction.ks))
            line += '{}\t!'.format(process_parameter(reaction.alpha))
            line += reaction.__repr__()+'\n'
            return line
        
        else:
            line = '{}, '.format(process_parameter(reaction.mode[1]))
            reaction_input_dict = self.get_reaction_dict()
            for key, coeff in reaction.reactants.items():
                reaction_input_dict[key]= '-{}'.format(process_parameter(coeff))
            for key, coeff in reaction.products.items():
                reaction_input_dict[key]= '{}'.format(process_parameter(coeff))
                
            line += ', '.join(process_parameter(value) for key, value in reaction_input_dict.items() )
            line += ', '
            line += '{}, '.format(process_parameter(reaction.kf))
            line += '{}, '.format(process_parameter(reaction.kb))
            line += '0.0e0, 0.0e0, 0.50\t!'
            line += reaction.__repr__()+'\n'
            return line
        
    def __repr__(self):
        line = ''
        for r in self.reactions:
            line += r.__repr__()
            line += '\n'
            
        return line
    
    def get_input(self):
        lines = []
        for reaction in self.reactions:
            lines.append(self.process_reaction(reaction))
        
        return lines
    
from shutil import copy
import os

class Voltammetry:
    """
    Usage:
    To use just a DC cyclic voltammetry
    >>> cv = DCVoltammetry(E_min = 1.0, E_max=-1.0, nu=75e-3)
    >>> volt = Voltammetry(objs=[cv])
    
    To add a AC voltammetry load
    >>> cv = DCVoltammetry(E_min = 1.0, E_max=-1.0, nu=75e-3)
    >>> ac_cv = ACVoltammetry(1,1.0e0,18.0e0)
    >>> volt = Voltammetry(objs=[cv, ac_cv])

    To load a voltage input from a text file
    >>> volt = Voltammetry(objs=['voltage.txt'])
    (This function has not been tested)
    
    """
    def __init__(self, objs=None):
        self.set_defaults()
        self.objs = objs
        self.ramp = 0
        if self.objs is not None:
            for obj in self.objs:
                if isinstance(obj, str):
                    self.from_file(obj)
                    self.ramp = 2
                elif obj.type=='DC':
                    self.dcvolt_params.update(obj.params)
                    self.ramp = 0
                elif obj.type=='AC':
                    self.acvolt_params.update(obj.params)
                    self.ramp = 0
                else:
                    raise KeyError('Did not understand type {}'.type(obj))
    
    def from_file(self):
        dirname = os.path.dirname(__file__)
        fname = os.path.join(dirname, 'Einput.txt')
        copy(self.obj, fname)

    def set_defaults(self):
        self.dcvolt_params = {'T':298.2, 'Rh':0.0e0, 'E_start':0.50, 'E_rev':-0.50, 'num_cycles':1, 'scan_rate':1.0e0}
        self.acvolt_params = {'num_sources': 1, 'amplitude':0.0e0, 'frequency':18.0e0}
      
    def get_input(self):
        lines = []
        for key, value in self.dcvolt_params.items():
            lines.append('{}\t!{}\n'.format(process_parameter(value), key))
            
        lines.append('12\t!points in time across n cycle \n') 
        lines.append('0\t! correct vscan and freq for DigiPot/FFT \n')
        lines.append('1\t! output type: 0=E,i,t; 1=DigiPot compatible \n')
        lines.append('0\t! EC type: 0 = Butler-Volmer, 1 = Marcus theory \n')
        lines.append('0\t! Pre-equilibrium switch: 0=stay with user entered, 1 = apply Pre-eqm operation \n')
        lines.append('0\t! fix number of timesteps (1 = yes; 0 = no) \n')
        lines.append('4000\t! Use a fixed number of timesteps rather than 2^N \n')
        lines.append('0.10\t! beta \n')
        lines.append('10.0\t! Dstar_min \n')
        lines.append('0.005e0\t! max voltage step \n')
        lines.append('25.6e0\t! time resolution experimentally to correct vscan/f (us) \n')
        lines.append('1\t! show debug output files as well as MECSimOutput.txt (1=yes; 0=no) \n')
        lines.append('{}\t! 0 = E_start=E_end, 1 = use advanced ramp below, 2=From file "EInput.txt"\n'.format(self.ramp))
        lines.append('2\t! Not used \n')
        lines.append('0.0\t!Not used E_start (V) \n')
        lines.append('0.0\t! Not use E_end (V) \n')
        lines.append('0.50\t! Not used E_rev - REPEAT for more complicated ramps \n')
        lines.append('-0.50\t! Not used E_rev - REPEAT for more complicated ramps \n')
        
        _line = '{}\t! number of AC sources to add (keep 1 with zero amplitude if want DC only) \n'
        lines.append(_line.format(process_parameter(self.acvolt_params['num_sources'])))
        _line = '{}, {}\t! AC sin wave: amp (mV), freq(Hz) (REPEAT) \n'
        lines.append(_line.format(process_parameter(self.acvolt_params['amplitude']),
                                  process_parameter(self.acvolt_params['frequency'] )))
        
        return lines
                    
class DCVoltammetry:
    """
    Usage:
    >>> cv = DCVoltammetry(E_start = 1.0, E_rev=-1.0, nu=75e-3)
    """
    def __init__(self, E_start=0.50, E_rev=-0.50, N = 1, nu=1.0e0, T=298.2, Rh=0.0e0):
        self.params = {'T':T, 'E_start':E_start, 'E_rev':E_rev, 'num_cycles':N, 'scan_rate':nu, 'Rh':Rh}
        self.type = 'DC'
    
class ACVoltammetry:
    """
    Usage:
    >>> cv = DCVoltammetry(E_min = 1.0, E_max=-1.0, nu=75e-3)
    >>> ac_cv = ACVoltammetry(1,1.0e0,18.0e0)
    """
    def __init__(self, num_sources = 1, amplitude=0.0e0, frequency=18.0e0):
        self.params= {'num_sources': num_sources, 'amplitude':amplitude, 'frequency':frequency}
        self.type = 'AC'
        
class Electrode:
    def __init__(self, params):
        self.params = params
        self.set_defaults()
        self.defaults.update(params)

    def set_defaults(self):
        keys = ['planar_surface_area','num_spheres','r_sphere','num_cylinders','r_cylinder','l_cylinder',\
                'cylinder_resolution','r_RDE','v_RDE','nu_RDE']
        values = [1.0e0, 1.0e0, 1.0e-4, 0.5e0, 0.001e0, 0.10e0, 100, 1.0e-1, 1.0e2,1.0e-5]
        self.defaults = dict(zip(keys, values))
     
    def get_input(self):
        lines = ['{}\t! Geometry type\n'.format(process_parameter(self.geometry))]
        for key, value in self.defaults.items():
            lines.append('{}\t! {} \n'.format(process_parameter(value), key))
            
        return lines

    def __repr__(self):    
        strings = [0,0,0,0]
        strings[0] = 'Planar electrode of area :{:.2E} cm2'.format(self.defaults['planar_surface_area'])
        strings[1] = 'Spherical electrode of radius :{:.2E} cm, {:.2f} spheres'
        strings[1] = strings[1].format(self.defaults['r_sphere'], self.defaults['num_spheres'])
        strings[2] = 'Cylindrical electrode with {} cylinders of radius :{:.2E} cm, length : {:.2E} cm'
        strings[2] = strings[2].format(self.defaults['num_cylinders'],self.defaults['r_cylinder'], self.defaults['l_cylinder'])
        strings[3] = 'RDE electrode of radius :{:.2E} cm, speed : {:.2E} rad/s, kinematic viscosity {:.2E} cm2/s'
        strings[3] = strings[3].format(self.defaults['r_RDE'],self.defaults['v_RDE'], self.defaults['nu_RDE'])

        return strings[int(self.geometry-1)]
    
class PlanarElectrode(Electrode):
    """
    Usage:
    >>> planar = PlanarElectrode()
    >>> print(planar)
    Planar electrode of area :1.00E+00 cm2

    """
    def __init__(self,area=1.0e0):
        self.geometry = 1
        kwargs = {'planar_surface_area':area}
        super().__init__(kwargs)
            
class SphericalElectrode(Electrode):
    """
    Usage:
    >>> spherical = SphericalElectrode()
    >>> print(spherical)
    Spherical electrode of radius :1.00E-04 cm, 0.50 spheres
    
    """
    def __init__(self,N=0.5, R=1.0e-4):
        self.geometry = 2
        kwargs = {'num_spheres':N, 'r_sphere':R}
        super().__init__(kwargs)

    
class CylindricalElectrode(Electrode):
    """
    >>> cylinder = CylindricalElectrode()
    >>> print(cylinder)
    Cylindrical electrode with 0.5 cylinders of radius :1.00E-04 cm, length : 1.00E-01 cm
    """
    
    def __init__(self,N=0.5, R=1.0e-4, l=0.1e0):
        self.geometry = 3
        kwargs = {'num_cylinders':N, 'r_cylinder':R,'l_cylinder':l }
        super().__init__(kwargs)
  
    
class RDEElectrode(Electrode):
    """
    >>> rde = RDEElectrode()
    >>> print(rde)
    RDE electrode of radius :1.00E-01 cm, speed : 1.00E+01 rad/s, kinematic viscosity 1.00E-05 cm2/s

    """
    
    def __init__(self, R=1.0e-1, v=0.1e2, nu=1.0e-5):
        self.geometry = 4
        kwargs = {'r_RDE':R,'v_RDE':v,'nu_RDE':nu }
        super().__init__(kwargs)
        
        
class Experiment:
    """
    Main eperiment class that is used to set up a voltammetry experiment.
    
    Usage:
    # first create species involved in the reaction
    >>> A = Specie('A', C0=1e-6)
    >>> B = Specie('B')
    >>> C = Specie('C')
    >>> D = Specie('D')
    >>> species = [A, B, C, D]
    
    # Define what is a reaction mechanism we are interested in simulating
    >>> R1 = ChargeTransfer({'A':1,'e':1},{'B':1},0.0)
    >>> R2 = ChemicalReaction({'B':1,'C':1},{'A':1,'D':1})
    >>> rxn = [R1, R2]
    
    # make a reaction mechanism from the species and reactions
    >>> mech = Mechanism(species, rxn)
    
    # Additionally, we can specify the electrode type (seel Electrodes above)
    >>> electrode = PlanarElectrode()
    
    # and a specific voltammetry experiment
    >>> cv = DCVoltammetry(E_min = 1.0, E_max=-1.0, nu=75e-3)
    >>> volt = Voltammetry(objs=[cv])  
    
    # create an exeriment
    >>> exp = Experiment(mech, electrode=electrode, voltammetry=volt)
    
    """
    
    def __init__(self, mechanism, **kwargs):
        self.mechanism = mechanism
        self.electrode = kwargs.get('electrode',PlanarElectrode())
        self.voltammetry = kwargs.get('voltammetry', Voltammetry())
        self.capacitance = kwargs.get('capacitance', Capacitance(0, [0,0,0,0]))
    
    def get_inpfile_lines(self):
        lines = {}
        lines['mechanism'] = self.mechanism.get_input()
        lines['voltammetry'] = self.voltammetry.get_input()
        lines['capacitance'] = self.capacitance.get_input()
        lines['species'] = [s.get_input() for s in self.mechanism.species]
        lines['electrode'] = self.electrode.get_input()
        
        inpfile = []
        inpfile.append(lines['voltammetry'][:-2])
        inpfile.append(lines['electrode'])
        inpfile.append(lines['voltammetry'][-2:])
        _line = '{}\t! Number of species \n'.format(process_parameter(len(lines['species'])))
        inpfile.append([_line])
        inpfile.append(lines['species'])
        inpfile.append(lines['capacitance'])
        inpfile.append(lines['mechanism'])
        
        inpfile_flatten = [item for sublist in inpfile for item in sublist]
        
        return inpfile_flatten
    
   
    
    
    
    
    
    
    
    
    
    