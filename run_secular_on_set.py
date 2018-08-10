from sim_secularmultiples import run_secular_multiples 
from amuse.lab import *
import sys	

file_path = sys.argv[1]

bodies = read_set_from_file(file_path, "hdf5")
print bodies
data = run_secular_multiples(bodies, 1.0 | units.yr, 50)
print data
