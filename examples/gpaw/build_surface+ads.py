from ase.build import fcc111, add_adsorbate
from ase.visualize import view #this allows us to quickly view atomic structures from within the script using ASE's lightweight GUI.
from ase.constraints import FixAtoms

#This script does not require a supercomputer, and should run from the login node (the node you first "land" on in a supercomputer) or from a desktop/laptop assuming you have python and the ASE python package installed.

#create the (111) surface slab using one of ASE's handy utility functions. Note that not all surfaces are so easy to create.
slab = fcc111('Pt', 
              size=(2,2,3), 
              vacuum=6.0,
              a=4.02,
              orthogonal=True)

#now we need to constrain some layers to simulate bulk metal. Let's start by constraining the lower 2/3 of the slab:

z_positions = [a.z for a in slab] #slab is an ASE Atoms object that acts like a list of atoms, allowing us to use Python's list comprehension syntax. 
#Each atom has some useful position attributes: a.position is the x,y,z coordinates, a.x is the x coordinate, a.y is the y coordinate, and a.z is the z coordinate.

z_positions.sort() #sort the z-positions from lowest to highest
cutoff_index = int(2*(len(slab)/3.)) #choose a cutoff index for the atom corresponding to 2/3 of the total atoms
z_cutoff = z_positions[cutoff_index] #find the z-coordinate of the lowest unconstrained atom

constraint = FixAtoms(mask=[a.z < z_cutoff for a in slab])  #create an ASE constraint object to fix the atoms. The "mask" argument allows using a list comprehension filter to define the atoms that are or aren't constrained. Read more here: https://wiki.fysik.dtu.dk/ase/ase/constraints.html

slab.constraints += [constraint] #add the constraints to the surface slab

add_adsorbate(slab,adsorbate='C',height=1.5,position='fcc') #Adding adsorbates isn't always this easy, but fcc(111) is a common surface then we can use this convenience function. See https://wiki.fysik.dtu.dk/ase/ase/build/surface.html for more.
add_adsorbate(slab,adsorbate='O',height=3.0,position='fcc') #Add the O 1.5 Anstrom above the C. We will let the optimizer find the correct bond length and position, this is just an initial guess.

slab.set_pbc((1,1,0)) #set periodic boundary conditions in x,y directions (x, y, z) - 0=not periodic, 1=periodic

#view(slab) # uncomment this line to view the slab as the script executes. You can move the "view" call up and down in the script to see how things change. Note that atoms with an "X" on them are constrained.

slab.write('Pt_111+CO.json') #we can write the atoms object to a .json file and then view it from the command line with "ase-gui" or load it into another script.

#It is generally better to think of the process of creating the atomic model as a separate step from actually calculating on it. This will make it easier to define arbitrarily complex atomic models, to transfer models between calculators, and to restart calculations.
