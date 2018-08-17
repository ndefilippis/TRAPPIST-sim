from multiprocessing import Pool
import closest_all_test, sys
import numpy as np

def apply_tuple_to_run(tuple):
    param, mass, seed = tuple
    initial = 0.009833817293648425
    percent_increase = 0.1
    closest_all_test.run(param, mass, seed, initial, percent_increase)

p = Pool(6)
start = int(sys.argv[1])
n = int(sys.argv[2])
mass = float(sys.argv[3])
end_mass = float(sys.argv[4])
samples = int(sys.argv[5])

dm = (end_mass - mass) / n
masses = np.logspace(start=np.log2(mass), stop=np.log2(end_mass), num=n, base=2.0)
masses = np.repeat(masses, samples)
labels = map("{0:03d}".format, range(start, start+n))
labels = np.repeat(labels, samples)
seeds = np.random.randint(100000, size=n*samples)
print masses
p.map(apply_tuple_to_run, zip(labels, masses, seeds))
