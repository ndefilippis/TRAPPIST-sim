#------------------------------------------------------------------------------#
#------------------------------trappist_system.py------------------------------#
#------------------------------------------------------------------------------#
#--------------------------Created by Nick DeFilippis--------------------------#
#------------------------------------------------------------------------------#
#-This program was created for the purpose of generating the TRAPPIST-1 system-#
#---------------------------in a particle set format---------------------------#
#------------------------------------------------------------------------------#
#------------------------Date Last Modified: 08/16/2018------------------------#
#------------------------------------------------------------------------------#

from amuse.lab import *
from amuse.ext.orbital_elements import rel_posvel_arrays_from_orbital_elements
import numpy as np
import random


#List all the data we do know

#Mass and radius of the star
star_mass = 0.089 | units.MSun
star_radius = 0.121 | units.RSun

#Masses and radii of the individual planets, from b-h
masses = [1.017, 1.156, 0.297, 0.772, 0.934, 1.148, 0.331] | units.MEarth
radii = [1.121, 1.095, 0.784, 0.910, 1.046, 1.148, 0.773] | units.REarth

#Orbital parameters of the planets around their host star
semimajor_axes = [0.01154775, 0.01581512, 0.02228038, 0.02928285, 0.03853361, 0.04687692, 0.06193488] | units.AU
eccentricities = [0.00622, 0.00654, 0.00837, 0.00510, 0.01007, 0.00208, 0.00567] | units.rad
inclinations = [89.56, 89.70, 89.89, 89.736, 89.719, 89.721, 89.796] | units.deg
argument_of_periapses = [336.86, 282.45, -8.73, 108.37, 368.81, 191.34, 338.92] | units.deg
mean_anomalies = [203.12, 69.86, 173.92, 347.95, 113.61, 265.08, 268.72] | units.deg

# Since we don't know longitude of ascending node,
# we'll need to randomly generate it
def gen_random_data(n, min, max, unit):
    data = [] | unit
    for _ in range(n):
        data.append(random.uniform(min, max) | unit)
    return data

# Generate the trappist-1 system from a given seed
def gen_trappist_system(seed):
   random.seed(seed)
   bodies = Particles()
   
   # Put the star at the center of the system
   star = Particle()
   star.mass = star_mass
   star.radius = star_radius
   star.position = [0, 0, 0] | units.AU
   star.velocity = [0, 0, 0] | units.kms
   bodies.add_particle(star)

   # Randomly generate the LANs
   longitude_of_ascending_nodes = gen_random_data(len(masses), -np.pi/20, np.pi/20, units.rad)

   # Get position and velocities by calculating them from orbital parameters.
   positions, velocities = rel_posvel_arrays_from_orbital_elements(star_mass, masses, semimajor_axes, eccentricities, mean_anomalies, inclinations.number, \
       longitude_of_ascending_nodes, argument_of_periapses.number, G=constants.G)
   for i in range(len(masses)):
       planet = Particle()
       planet.mass = masses[i]
       planet.radius = radii[i]
       planet.position = positions[i]
       planet.velocity = velocities[i]
       bodies.add_particle(planet)

   return bodies

if __name__ == "__main__":
    gen_trappist_system(10)
       
