# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 14:59:17 2017

@author: benjamin
"""

#imports
from ase.build import surface, add_adsorbate, molecule
from ase import io
from ase.visualize import view
from ase.constraints import FixAtoms
from ase import Atoms
import numpy as np


y_num = 2
layers = 4
vacuum = 6
unfrozen = 3
displ = 2.0
rutile = io.read('Bulk_Rutile.json')
slb = surface(rutile,indices=(1,1,0),layers=layers,vacuum=vacuum)
n2o2 = Atoms('ON2O',[(0, 0, 0),
                      (0.1936, 1.144, 0),
                     (2.3, 1.144, 0),
                       (2.624, 0, 0)])
n2o2.rotate('x',3*np.pi/2, center=(0,0,0))
n2o2.rotate('z',np.pi/2, center=(0,0,0))

#move around some unruly atoms
pos = slb[len(slb)-1].get('position') +[0,0,-3.288*layers]
slb[len(slb)-1].set('position',pos)
pos = slb[len(slb)-2].get('position') +[0,0,-3.288*layers]
slb[len(slb)-2].set('position',pos)
pos = slb[len(slb)-4].get('position') +[0,0,-3.288*layers]
slb[len(slb)-4].set('position',pos)


slb = slb*(1,y_num,1)
#view(slb)

if (layers % 2 ==0):
    add_adsorbate(slb,molecule('N2'), displ, position=(3.2,3))
else:
    add_adsorbate(slb,molecule('N2'), displ, position=(0.1,3))
#z_cutoff = 8.29+3.289*(layers-4)+vacuum-3
z_cutoff = 3.288*(layers-unfrozen)+vacuum
d = FixAtoms(mask=[a.z < z_cutoff for a in slb])
slb.set_constraint(d)
#rutile.rattle()
#slb.write_potcar()
#slb.initialize
view(slb)
slb.write('POSCAR')
