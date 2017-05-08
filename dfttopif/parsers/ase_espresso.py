from pypif.obj.common.property import Property

from .base import DFTParser, Value_if_true
from .pwscf import PwscfParser
import os
from pypif.obj.common.value import Value
from ase import Atoms, io
import json
from ase.io.jsonio import encode
from collections import OrderedDict
from ase.constraints import dict2constraint

class AseEspressoParser(PwscfParser):
    '''
    Parser for PWSCF calculations
    '''
    
    def get_name(self): return "ASE-ESPRESSO"

    def test_if_from(self, directory):
        '''Look for PWSCF input and output files'''
        self.outputf = ''
        files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
        subdirs = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]
        for sd in subdirs:
            files += [os.path.join(sd,f) for f in os.listdir(os.path.join(directory,sd)) if os.path.isfile(os.path.join(directory,sd, f))]
        for f in files:
            try:
                if self._get_line('Program PWSCF', f, basedir=directory, return_string=False):
                    self.outputf = f
                if self.outputf: 
                    return True
            except UnicodeDecodeError:
                pass
        return False
        
    def get_setting_functions(self):
        '''Get a dictionary containing the names of methods
        that return settings of the calculation
        
        Returns:
            dict, where the key is the name of the setting,
                and the value is function name of this parser
        '''
        return {
            'XC Functional':'get_xc_functional',
            'Relaxed':'is_relaxed',
            'Cutoff Energy':'get_cutoff_energy',
            'k-Points per Reciprocal Atom':'get_KPPRA',
            'Spin-Orbit Coupling':'uses_SOC',
            'DFT+U':'get_U_settings',
            'vdW Interactions':'get_vdW_settings',
            'Psuedopotentials':'get_pp_name',
            'INCAR':'get_incar',
            'ASE atoms':'get_atoms',
            'POSCAR':'get_poscar',
        }

    def get_KPPRA(self):
        '''Determine the no. of k-points in the BZ (from the input) times the
        no. of atoms (from the output)'''
        return None

    def get_vdW_settings(self):
        '''Determine the vdW type if using vdW xc functional or correction
        scheme from the input otherwise'''
        return None
        
    def get_total_energy(self):
        '''Determine the total energy from the output'''
        hist = self.get_total_energy_histogram()
        with open(os.path.join(self._directory, self.outputf)) as fp:
            # reading file backwards in case relaxation run
            for line in reversed(fp.readlines()):
                if "!" in line and "total energy" in line:
                    energy = line.split()[4:]
                    return Property(scalars=float(energy[0]), units=energy[1], histogram=hist)
            raise Exception('%s not found in %s'%('! & total energy',os.path.join(self._directory, self.outputf)))

    def get_total_energy_histogram(self):
        '''Return the 2000 value BEEF ensemble'''
        with open(os.path.join(self._directory, self.outputf)) as fp:
            txt = fp.read()
            _,E_total = txt.rsplit('total energy              =',1)
            E_total,_ = E_total.split('Ry',1)
            E_total = float(E_total.strip())
            E_total *= 13.605698
            _, ens = txt.rsplit('BEEFens 2000 ensemble energies',1)
            ens,_ = ens.split('BEEF-vdW xc energy contributions',1)
            ens.strip()
            ens_ryd = []
            for Ei in ens.split('\n'):
                if Ei.strip():
                    Ei = float(Ei.strip()) + E_total
                    ens_ryd.append(Ei)
            ens_ryd = ens_ryd
            return ens_ryd

            # reading file backwards in case relaxation run
            log_lines = reversed(fp.readlines())
            
            for line_number,line in enumerate(log_lines):
                if "BEEFens" in line and "ensemble energies" in line:
                    ensemble_location = range(line_number-2001,line_number-1)
                    ens = []
                    for ens_line in ensemble_location[::-1]:
                        
                        ens.append(float(log_lines[ens_line].strip()))
                    return ens
            raise Exception('%s not found in %s'%('BEEFens & ensemble energies',os.path.join(self._directory, self.outputf)))

    def atoms_to_dict(self):
        """
        converts an atoms object into a dictionary of the properties. Mostly 
        copied from the Kitchin group
        """
        atoms = io.read(os.path.join(self._directory, 'converged_slab.traj'),)
#        atoms = io.read(os.path.join(self._directory, self.outputf),-1,format='espresso-out')
        
        d = OrderedDict(atoms=[{'symbol': atom.symbol,
                            'position': json.loads(encode(atom.position)),
#                            'position': atom.position,
                            'tag': atom.tag,
                            'index': atom.index,
                            'charge': atom.charge,
                            'momentum': json.loads(encode(atom.momentum)),
#                            'momentum': atom.momentum,
                            'magmom': atom.magmom}
                           for atom in atoms],
                    cell=atoms.cell,
                    pbc=atoms.pbc,
                    info=atoms.info,
                    constraints=[c.todict() for c in atoms.constraints])
                        # redundant information for search convenience.
        d['natoms'] = len(atoms)
        cell = atoms.get_cell()
        if cell is not None and np.linalg.det(cell) > 0:
            d['volume'] = atoms.get_volume()
    
        d['mass'] = sum(atoms.get_masses())
    
        syms = atoms.get_chemical_symbols()
        d['chemical_symbols'] = list(set(syms))
        d['symbol_counts'] = {sym: syms.count(sym) for sym in syms}
        d['spacegroup'] = spglib.get_spacegroup(atoms)
        return Property(scalars=float(energy[0]), units=energy[1], histogram=hist)
        
    def dict_to_atoms(doc):
        """
        Takes in a PIF dictionary and creates an atoms object. Mostly copied 
        from Kitchin group.
        """
        atoms = Atoms([Atom(atom['symbol'],
                                atom['position'],
                                tag=atom['tag'],
                                momentum=atom['momentum'],
                                magmom=atom['magmom'],
                                charge=atom['charge'])
                           for atom in doc['atoms']['atoms']],
                          cell=doc['atoms']['cell'],
                          pbc=doc['atoms']['pbc'],
                          info=doc['atoms']['info'],
                          constraint=[dict2constraint(c) for c in doc['atoms']['constraints']])
    
#        from ase.calculators.singlepoint import SinglePointCalculator
#        results = doc['results']
#        calc = SinglePointCalculator(energy=results.get('energy', None),
#                                     forces=results.get('forces', None),
#                                     stress=results.get('stress', None),
#                                     atoms=atoms)
#        atoms.set_calculator(calc)
        return atoms
