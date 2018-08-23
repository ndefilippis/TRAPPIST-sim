from multiprocessing import Pool
import angles_test, sys

p = Pool(6)
start = int(sys.argv[1])
n = int(sys.argv[2])
p.map(angles_test.run, map("{0:03d}".format, range(start, start+n)))

