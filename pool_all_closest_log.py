from multiprocessing import Pool
import closest_all_test, sys
import numpy as np
from amuse.lab import units
from amuse.ic.brokenimf import MultiplePartIMF

def apply_tuple_to_run(tuple):
    param, mass, seed = tuple
    e_changes_to_check = [0.01, 0.001]
    initial = 0.009833817293648425
    for e in e_changes_to_check:
        closest_all_test.run(param+str(e), mass, seed, initial, e)

def get_mass_from_kroupa(n):
    return MultiplePartIMF(mass_boundaries=[0.01, 0.08, 0.5, 15.0] | units.MSun, \
        alphas=[-0.3, -1.3, -2.3], random=True).next_mass(n)

p = Pool(6)
start = int(sys.argv[1])
n = int(sys.argv[2])

masses = get_mass_from_kroupa(n)
labels = map("{0:03d}".format, range(start, start+n))
seeds = np.random.randint(100000, size=n)
masses = np.logspace(np.log2(0.01), np.log2(10.0), base=2.0, num=n, endpoint=True)
p.map(apply_tuple_to_run, zip(labels, masses, seeds))
