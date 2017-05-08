from dfttopif import directory_to_pif

from dfttopif.parsers.ase_espresso import PIF_to_calculator


import json

data = directory_to_pif('examples/ase_espresso')
dd = data.as_dictionary()

json.dump(dd, open('TiO2_N2.pif','w'))

PIF_to_calculator(data)

a = json.dumps(data.as_dictionary())
f = open('PIF_view.json','w')

f.write(a)
#print(data.properties[-1].histogram)
