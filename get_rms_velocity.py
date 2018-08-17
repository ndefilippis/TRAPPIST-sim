from amuse.lab import units, constants
import os, pickle
import numpy as np

path = '/home/draco/mgiovinazzi/ce_databases'
for dir in os.listdir(path):
   for dd in os.listdir(os.path.join(path, dir)):
     f = open(os.path.join(path, dir, dd))
     q = pickle.load(f)
     vs = []
     if "W3" in dd:
       continue
     for key in q.keys():
       for en in q[key]:
         if len(en) == 2:
           rel_vel = en[0].velocity - en[1].velocity
	   rel_pos = en[0].position - en[1].position
	   KE = 0.5 * en[0].mass * rel_vel.norm() ** 2
	   PE = -constants.G * en[1].mass * en[0].mass / rel_pos.norm()
	   v_inf = np.sqrt(2 * (KE-PE) / en[0].mass)
	   vs.append(v_inf)
     aa = 0 | (units.kms)
     for v in vs:
       aa += v
     print aa / len(vs)
