#This is a heavy calculation that requires the use of a supercomputer. It must be submitted using the run.sh script.

#Import all necessary modules
from gpaw import GPAW
from gpaw import Mixer
from gpaw import FermiDirac
from gpaw.poisson import PoissonSolver
from gpaw.dipole_correction import DipoleCorrection
from ase import io
from ase.parallel import paropen as open #ensures that open works in parallel environment
from ase.optimize import QuasiNewton #geometry optimization algorithm; QuasiNewton links to BFGS line search, which is the best general-purpose optimizer, but other options are available: https://wiki.fysik.dtu.dk/ase/ase/optimize.html

#setup calculator

mixer = Mixer(beta=0.07, nmaxold = 3, weight = 120.0)
fermidirac = FermiDirac(0.1)
convergence = {'energy': 0.0005, 
        'density': 1e-4,
        'eigenstates': 4e-8,
        'bands': 'occupied'}

p = PoissonSolver(nn = 3, relax = 'J', eps=1e-12)
poisson = DipoleCorrection(p,2) #add dipole correction

calcargs = dict(xc='BEEF-vdW',
        poissonsolver=poisson,
        kpts=(4, 4, 1), #only need 1 kpt in z-direction
        h=0.16,
        mixer=mixer,
        eigensolver='rmm-diis',
        mode = 'fd',
        spinpol=False,
        charge=0,
        setups = 'paw',
        maxiter=500,
        occupations=fermidirac,
        convergence = convergence,
        txt ='gpaw.log')
#Many of these arguments are technically not necessary, or set to default values, but it is a good idea to understand what they mean.
#In general I think it is better to *explicitly* define all arguments to a calculator definition rather than falling back on defaults
#that may not be the right choice or may change in different versions of the underlying calculator.
#Explictly stating the arguments forces you to think about all the assumptions you are really making. Note that this is still not a fully
#complete set of arguments.

calc = GPAW(**calcargs)

atoms = io.read('Pt_111+CO.json') #Read in the structure built by the other script

atoms.set_calculator(calc)

relax = QuasiNewton(atoms,logfile='opt.log',trajectory='opt.json',restart='opt.pckl')
#set up the optimization algorithm. It has a logfile, a "trajectory" file that tracks progress, and a restart file in case the algorithm has to be restarted.
relax.run(fmax=0.05) #execute the relaxation algorithm. It will run until the maximum force on any atom is <0.05 eV/Angstrom.

energy = atoms.get_potential_energy() #this is the potential energy of the electrons as computed by DFT. It will be closely related to the enthalpy.

atoms.write('Pt_111+CO-converged.json')

f = open('converged.log','w')
f.write(str(energy))
f.close()
#This is not technically necessary, but it is often helpful to have a file that confirms whether or not a simulation has converged.
#Now we know that if 'converged.log' exists then the calculation has finished without having to check the queue.
