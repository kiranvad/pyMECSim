"""
Change logs
-----------
Version 2 :
    1. reactions now requires input to be list of tuples with (species, coefficent)
    2. mechanism requires only list of inputs

"""

import pdb
import warnings
import itertools
from shutil import copy
import os
import numpy as np

from .utils import deprecated

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
        self.surface_confined = surface_confined
        self._set_specie_type()
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
            
    def _set_specie_type(self):
        if self.surface_confined:
            self.specie_type = 'Surface confined'
            self._specie_type = 1
        else:
            self.specie_type = 'Solution'
            self._specie_type = 0
    
    def __repr__(self):
        line = '{} specie {} with Diffusion coefficent {} {} , intial concentration {} {}'
        return line.format(self.specie_type, self.name, 
                           process_parameter(self.D), self.D_units, process_parameter(self.C0), self.C_units)
        
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
        """
        Main reaction base class that takes list of tuples in reactants and products and produces a pymecsim reaction.
        Each tuple is (specie, coefficient) in the reaction. where species needs to be `pymecsim::Specie` type.
        
        _get_formula : produces the reaction as a formula and returns string
        _to_dict : converts the reactants/products from tuple to dict 
        """
        self._check_all_species(reactants)
        self.reactants = reactants
        self._check_all_species(products)
        self.products = products
        
        self.reactants_dict = self._to_dict(self.reactants)
        self.products_dict = self._to_dict(self.products)
        self._check_surface_reaction()
        self.num_species = len(self.reactants_dict) + len(self.products_dict)
    
    def __repr__(self):
        reaction = self._get_formula()  
        reaction += self.params
        
        return reaction
    
    def _get_formula(self):
        lhs = ['{} {}'.format(value, key) for key, value in self.reactants_dict.items()]
        rhs = ['{} {}'.format(value, key) for key, value in self.products_dict.items()]
        
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
        
        return reaction
    

    def _to_dict(self, x):
        d = {}
        for specie, coeff in x:
            if specie=='e':
                d[specie] = coeff
            else:
                d[specie.name] = coeff
            
        return d
    
    def _check_surface_reaction(self):
        
        surface_reactants = 0
        for (r,_) in self.reactants:
            if r!='e':
                if r._specie_type==1:
                    surface_reactants += 1
                
        surface_products = 0
        for (p,_) in self.products:
            if p._specie_type==1:
                surface_products += 1
        if surface_reactants!=surface_products:
            message = 'Number of surface confined species on both sides of'\
            'the reactions needs to be same for \n' 
            message += self._get_formula()
            message += '\nGiven {} as reactants, {} as products'.format(surface_reactants, surface_products)

            raise Exception(message) 
    
    def _check_all_species(self, specielist):
        from pymecsim import Specie
        
        for i,_ in specielist:
            if isinstance(i, str) and i=='e':
                continue
            elif isinstance(i, Specie):
                continue
            else:
                raise RuntimeError('Species {} is not recognized `pymecsim::Specie`')
    
class ChargeTransfer(Reaction):
    """
    A subtype of `pymecsim::Reaction` class for charge transfer reactions.
    In this rections, one needs to use electron as a specie with a number of electrons transfered as coefficient.
    To represent electron as a tuple to be used in `Reaction` use the tuple: ('e',n) where n is number of electron transfered.
    
    Usage:
    >>> R1 = ChargeTransfer([(A,1),('e',2)],[(B, 1)],0.0)
    >>> print(R1) 
    Charge Transfer : 1 A + 2 e <=> 1 B  ks= 1.00E+04, E0 = 0.00E+00, alpha = 0.50
    """
    def __init__(self, reactants, products, E0=0.0, ks=1e4, alpha=0.5):
        self.ks = ks
        if isinstance(self.ks, float):
            if self.ks>1e14:
                warnings.warn('reaction rate {:.2E} may be too high'.format(self.ks))
                
        self.E0 = E0
        self.alpha = alpha
        self.mode = ['Charge Transfer', 0]
        self.params = self.get_params_as_string()
        super().__init__(reactants, products)

    def get_params_as_string(self):
        params = '  ks= {}, '.format(process_parameter(self.ks))
        params += 'E0 = {}, '.format(process_parameter(self.E0))
        params += 'alpha = {}'.format(process_parameter(self.alpha))
        
        return params

class ChemicalReaction(Reaction):
    """
    A subtype of `pymecsim::Reaction` class for chemical reactions.
    usage:
    >>> R2 = ChemicalReaction({'A':1,'B':1},{'C':1,'D':1}) (depreceated)
    >>> R2 = ChemicalReaction([(A,1),(B,1)],[(C,1),(D,1)])
    >>> print(R2)
    Chemical Reaction : 1 A + 1 B <=> 1 C + 1 D  kf= 1.00E+04, kb= 1.00E+04  
    """
    def __init__(self, reactants, products, kf=1e4, kb=1e4):
        self.kf = kf
        if isinstance(self.kf, float):
            if self.kf>1e14:
                warnings.warn('reaction rate {:.2E} may be too high'.format(self.kf))
                
        self.kb = kb
        if isinstance(self.kb, float):
            if self.kb>1e14:
                warnings.warn('reaction rate {:.2E} may be too high'.format(self.kb))
        self.mode = ['Chemical Reaction', 2]
        self.params = self.get_params_as_string()
        super().__init__(reactants, products)

     
    def get_params_as_string(self):
        params = '  kf= {}, '.format(process_parameter(self.kf))
        params += 'kb= {} '.format(process_parameter(self.kb))
        
        return params    

class CatalyticReaction(Reaction):
    """
    A subtype of `pymecsim::Reaction` class for chemical reactions.
    Usage:
    >>> R3 = CatalyticReaction({'B':1},{'C':1}) (depreceated)
    >>> R3 = CatalyticReaction([(B,1)],[(C,1)])
    >>> print(R3)
    Catalytic Reaction : 1 B <=> 1 C  kf= 1.00E+04, kb= 1.00E+04 
    """
    def __init__(self, reactants, products, kf=1e4, kb=1e4):
        self.kf = kf
        if isinstance(self.kf, float):
            if self.kf>1e14:
                warnings.warn('reaction rate {:.2E} may be too high'.format(self.kf))
                
        self.kb = kb
        if isinstance(self.kb, float):
            if self.kb>1e14:
                warnings.warn('reaction rate {:.2E} may be too high'.format(self.kb))
                
        self.mode = ['Catalytic Reaction', 1]
        self.params = self.get_params_as_string()
        super().__init__(reactants, products)

     
    def get_params_as_string(self):
        params = '  kf= {}, '.format(process_parameter(self.kf))
        params += 'kb= {} '.format(process_parameter(self.kb))
        
        return params 
           
class Mechanism:
    """
    A mechanism class that can be created by passing a list of `pymecsim::Reaction` class objects
    Usage:
    >>> A = Specie('A', C0=1e-6)
    >>> B = Specie('B')
    >>> C = Specie('C')
    >>> D = Specie('D')
    >>> species = [A, B, C, D]

    >>> R1 = ChargeTransfer({'A':1,'e':1},{'B':1},0.0) (depreceated)
    >>> R1 = ChargeTransfer([(A, 1), ('e', 1)],[(B, 1)],0.0) 
    
    >>> R2 = ChemicalReaction({'B':1,'C':1},{'A':1,'D':1}) (depreceated)
    >>> R2 = ChemicalReaction([(B, 1),(C,1)],[(A, 1),(D, 1)])    
    >>> rxn = [R1, R2]
    >>> mech = Mechanism(rxn)
    >>> print(mech)  

    Charge Transfer : 1 A + 1 e <=> 1 B  ks= 1.00E+04, E0 = 0.00E+00, alpha = 0.50
    Chemical Reaction : 1 B + 1 C <=> 1 A + 1 D  kf= 1.00E+04, kb= 1.00E+04
    """
    def __init__(self, reactions):
        self.reactions = reactions
        self.species = self._get_species_from_reactions()    
        self.num_species = len(self.species)
        self.input = self.get_input()

    def _get_reaction_dict(self):
        """ Creates a dummy dictornay with species names as keys
        
        makes sure that species are labelled in the same order in reactions and 
        species parts of  input to MECSim
        
        """
        reaction_input_dict = {}
        for s in self.species:
            reaction_input_dict[s.name]=0
        
        return reaction_input_dict
        
    def process_reaction(self, reaction):
        """Main function to process the reaction into a MECSim readable reaction line
        
        """
        if reaction.mode[1]==0:
            line = '{}, '.format(process_parameter(reaction.mode[1]))
            reaction_input_dict = self._get_reaction_dict()
            for key, coeff in reaction.reactants_dict.items():
                reaction_input_dict[key]= '-{}'.format(process_parameter(coeff))
            for key, coeff in reaction.products_dict.items():
                reaction_input_dict[key]='{}'.format(process_parameter(coeff)) 
            del reaction_input_dict['e']
            
            line += ', '.join(process_parameter(value) for key, value in reaction_input_dict.items() )
            line += ', 0.0e0, 0.0e0, '
            line += '{}, '.format(process_parameter(reaction.E0))
            line += '{}, '.format(process_parameter(reaction.ks))
            line += '{}\t!'.format(process_parameter(reaction.alpha))
            line += reaction._get_formula()+'\n'
            return line
        
        else:
            line = '{}, '.format(process_parameter(reaction.mode[1]))
            reaction_input_dict = self._get_reaction_dict()
            for key, coeff in reaction.reactants_dict.items():
                reaction_input_dict[key]= '-{}'.format(process_parameter(coeff))
            for key, coeff in reaction.products_dict.items():
                reaction_input_dict[key]= '{}'.format(process_parameter(coeff))
                
            line += ', '.join(process_parameter(value) for key, value in reaction_input_dict.items() )
            line += ', '
            line += '{}, '.format(process_parameter(reaction.kf))
            line += '{}, '.format(process_parameter(reaction.kb))
            line += '0.0e0, 0.0e0, 0.50\t!'
            line += reaction._get_formula()+'\n'
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
    
    @deprecated
    def _check_species_reactions(self):
        species_names = set([s.name for s in self.species])
        species_in_reactions = []
        for rxn in self.reactions:
            for r in rxn.reactants:
                species_in_reactions.append(r)
            for p in rxn.products:
                species_in_reactions.append(p) 
        pdb.set_trace()
        species_in_reactions = set(species_in_reactions)
        species_in_reactions.remove('e')
        if not species_names==species_in_reactions:
            raise Exception('species in reactions and the list of Species provided did not match')
         
    def _get_species_from_reactions(self):
        solution_species = []
        surface_species = []
        for reaction in self.reactions:
            for side in [reaction.reactants, reaction.products]:
                for specie,_ in side:
                    if not specie=='e':
                        if not specie.surface_confined:
                            if specie not in solution_species:
                                solution_species.append(specie)
                        else:
                            if specie not in surface_species:
                                surface_species.append(specie)
                        
        species = list(itertools.chain(solution_species, surface_species))
        
        return species
            
        
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
    def __init__(self, objs=[], N=12):
        """
        N : Number of spatial points as a power of two
        """
        self.N = N
        self.set_defaults()
        self.objs = objs
        self.ramp = 0
        self._check_objs()
        self._check_errors()

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
            
        lines.append('{}\t!points in time across n cycle \n'.format(process_parameter(self.N))) 
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
    
    def _check_objs(self):
        self.obj_mode = []
        if len(self.objs)==0:
            msg = 'using the default DC Voltammetry loading: \n'
            msg += '{}'.format(self.dcvolt_params)
            warnings.warn(msg)
        else:
            for obj in self.objs:
                if isinstance(obj, str):
                    self.from_file(obj)
                    self.ramp = 2
                elif obj.type=='DC':
                    self.dcvolt_params.update(obj.params)
                    self.ramp = 0
                    self.obj_mode.append('DC')
                elif obj.type=='AC':
                    self.acvolt_params.update(obj.params)
                    self.ramp = 0
                    if 'DC' not in self.obj_mode:
                        msg = 'DC voltammetry is not specied using default loadings in : \n'
                        msg += '{}'.format(self.dcvolt_params)
                        warnings.warn(msg)
                        
                    self.obj_mode.append('AC')
                    
                else:
                    raise KeyError('Did not understand type {}'.type(obj))
    
    def _check_errors(self):
        # check if voltage sweep is None
        E1 = self.dcvolt_params['E_start']
        E2 = self.dcvolt_params['E_rev']
        if np.isclose(E1, E2):
            raise RuntimeError('Voltage sweep has identical start and reverse voltages.')
     
    
    
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
    
   
    
    
    
    
    
    
    
    
    
    