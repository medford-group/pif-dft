from dfttopif import directory_to_pif
import json

data = directory_to_pif('examples/ase_espresso')
dd = data.as_dictionary()

json.dump(dd, open('TiO2_N2.pif','w'))

print(data.as_dictionary())

#print(data.properties[-1].histogram)
