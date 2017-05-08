#This is a heavy calculation that requires the use of a supercomputer. It must be submitted using the run.sh script.

#Import all necessary modules
from ase import io
from ase.parallel import paropen as open #ensures that open works in parallel environment
from ase.optimize import QuasiNewton #geometry optimization algorithm; QuasiNewton links to BFGS line search, which is the best general-purpose optimizer, but other options are available: https://wiki.fysik.dtu.dk/ase/ase/optimize.html
from espresso import espresso
import numpy as np

#setup calculator

calcargs = dict(xc='BEEF-vdW',
        kpts=(4, 4, 1), #only need 1 kpt in z-direction
        pw=400.,
        dw=4000.,
        spinpol=True,
        beefensemble=True,
        printensemble=True,
        convergence={'energy':1e-6,
                    'mixing':0.05,
                    'maxsteps':1000,
                    'diag':'david'},
        startingwfc='atomic',
        smearing='fd', #fermi-dirac electron smearing
        sigma=0.1, #smearing width
        dipole={'status':True}, #dipole corrections True turns them on
        #parflags='-nk 2',
        outdir ='esp.log')

calc = espresso(**calcargs)

atoms = io.read('POSCAR') #Read in the structure built by the other script


atoms.set_calculator(calc)
init_mag = np.zeros(len(atoms))
init_mag[-1] = 0.1
atoms.set_initial_magnetic_moments(magmoms=init_mag)

relax = QuasiNewton(atoms,logfile='opt.log',trajectory='opt.json',restart='opt.pckl')
#set up the optimization algorithm. It has a logfile, a "trajectory" file that tracks progress, and a restart file in case the algorithm has to be restarted.
relax.run(fmax=0.05) #execute the relaxation algorithm. It will run until the maximum force on any atom is <0.05 eV/Angstrom.

energy = atoms.get_potential_energy() #this is the potential energy of the electrons as computed by DFT. It will be closely related to the enthalpy.

atoms.write('converged_slab.traj')

f = open('converged.log','w')
f.write(str(energy))
f.close()
#This is not technically necessary, but it is often helpful to have a file that confirms whether or not a simulation has converged.
#Now we know that if 'converged.log' exists then the calculation has finished without having to check the queue.
