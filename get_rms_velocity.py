'''
Find the root mean squared velocity given a particular close encounter
database.
Created by: Nick DeFilippis, Last Modified: 8/22/2018
'''

from amuse.lab import units, constants
import os, pickle, sys
import numpy as np

if len(sys.argv) < 2:
    path = '/home/draco/mgiovinazzi/ce_databases'
else:
    path = sys.argv[1]

for dir in os.listdir(path):
   for dd in os.listdir(os.path.join(path, dir)):
     f = open(os.path.join(path, dir, dd))
     q = pickle.load(f)
     vs = []
     if "W3" in dd: # Separate by density parameter
       continue
     for key in q.keys():
       for en in q[key]:
         if len(en) == 2: # Only look at two bodies interactng
           rel_vel = en[0].velocity - en[1].velocity

       # Account for energy effects
	   rel_pos = en[0].position - en[1].position
	   KE = 0.5 * en[0].mass * rel_vel.norm() ** 2
	   PE = -constants.G * en[1].mass * en[0].mass / rel_pos.norm()
	   v_inf = np.sqrt(2 * (KE-PE) / en[0].mass)
	   vs.append(v_inf)

     # Average of all velocities
     aa = 0 | (units.kms)
     for v in vs:
       aa += v
     print aa / len(vs)
