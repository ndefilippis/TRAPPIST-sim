from multiprocessing import Pool
import angles_test

p = Pool(3)
p.map(angles_test.run, map(str, range(120, 180)))
