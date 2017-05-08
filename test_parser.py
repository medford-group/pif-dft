from dfttopif import directory_to_pif

data = directory_to_pif('examples/ase_espresso')

print(data.histogram)
