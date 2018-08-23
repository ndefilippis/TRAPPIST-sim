# Short script to run SmallN on the TRAPPIST-1 system
# Create by Nick DeFilippis, Last Modified: 8/21/2018

from trappist_system import gen_trappist_system
from sim_encounters import run_collision
from amuse.lab import *


bodies = gen_trappist_system(10)

converter = nbody_system.nbody_to_si(bodies.mass.sum(), 0.2 | units.RSun)

print run_collision(bodies, 2|units.yr, 1|units.hour, "001", "ce_directory/rebuild/TRAPPIST/000/1", converter=converter)
