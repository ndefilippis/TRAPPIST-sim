#------------------------------------------------------------------------------#
#----------------------------closest_all_test.py------------------------------#
#------------------------------------------------------------------------------#
#--------------------------Created by Nick DeFilippis--------------------------#
#------------------------------------------------------------------------------#
#---This program was created for the purpose of finding the closest approach---#
#-----distance needed to account for the a certain percentage change in the----#
#-----------------maximum eccentricity of the outermost planet-----------------#
#------------------------------------------------------------------------------#
#------------------------Date Last Modified: 08/16/2018------------------------#
#------------------------------------------------------------------------------#

import math, sys, os
import numpy as np
from amuse.lab import *
from amuse.community.kepler.interface import kepler
from sim_secularmultiples import run_secular_multiples
from sim_encounters import run_collision
from trappist_system import gen_trappist_system

def preform_EulerRotation():
    ''' Preforms a randomly oriented Euler Transformation to a set of AMUSE Particles.
        particle_set: AMUSE particle set which it will preform the transform on.
        !! Based on James Arvo's 1996 "Fast Random Rotation Matrices"
        !! https://pdfs.semanticscholar.org/04f3/beeee1ce89b9adf17a6fabde1221a328dbad.pdf
    '''
# First: Generate the three Uniformly Distributed Numbers (Two Angles, One Decimal)
    n_1 = np.random.uniform(0.0, math.pi*2.0)
    n_2 = np.random.uniform(0.0, math.pi*2.0)
    n_3 = np.random.uniform(0.0, 1.0)
# Second: Calculate Matrix & Vector Values
    c1 = np.cos(n_1)
    c2 = np.cos(n_2)
    s1 = np.sin(n_1)
    s2 = np.sin(n_2)
    r3 = np.sqrt(n_3)
    R = [[  c1,  s1, 0.0],
         [ -s1,  c1, 0.0],
         [ 0.0, 0.0, 1.0]]
    V = [[c2*r3],
         [s2*r3],
         [np.sqrt(1-n_3)]]
# Third: Create the Rotation Matrix

    # This was the old rotation matrix calculation...
    #rotate = (np.outer(V, V) - np.dot(np.eye(3),(R)))

    # But here is the new one which more correctly implements the equations from the paper referenced above...
    rotate = (2 * np.dot(np.outer(V, V), R) - np.dot(np.eye(3), R))
    return rotate

def extract_EulerAngles(rotation_matrix):
    #psi = np.arcsin(X_prime[1] / np.sqrt(1 - X_prime[2] ** 2))
    if rotation_matrix[2, 0] < 1:
        if rotation_matrix[2, 0] > -1:
	    theta = np.arcsin(-rotation_matrix[2, 0])
	    psi = np.arctan2(rotation_matrix[1, 0], rotation_matrix[0, 0])
	    phi = np.arctan2(rotation_matrix[2, 1], rotation_matrix[2, 2])
	else:
	    theta = np.pi / 2
	    psi = -np.arctan2(-rotation_matrix[1, 2], rotation_matrix[1, 1])
	    phi = 0
    else:
        theta = -np.pi / 2
	psi = np.arctan2(-rotation_matrix[1, 2], rotation_matrix[1,1])
	phi = 0
    #phi = np.arcsin(Y_prime[2] / np.sqrt(1 - X_prime[2] ** 2))
    return psi, theta, phi

def rotation_matrix_from_angles(psi, theta, phi):
    cs = np.cos([psi, theta, phi])
    ss = np.sin([psi, theta, phi])

    c = cs[0]
    s = ss[0]
    r1 = np.matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    c = cs[1]
    s = ss[1]
    r2 = np.matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    c = cs[2]
    s = ss[2]
    r3 = np.matrix([[1, 0, 0], [0, c, -s], [0, s, c]])

    return r1 * r2 * r3

def get_largest_delta_e(bodies, perturber, timescale, converter, n, param):
    path = "../ce_directory/close"+param+"/TRAPPIST/1"
    if not os.path.exists(path): os.makedirs(path)
    new_bodies = run_collision(bodies, timescale, timescale / 100, str(n), path, converter=converter)

    #print new_bodies
    t, a, e, i, _ = run_secular_multiples(new_bodies-perturber, 10000.0|units.yr, 300)

    es = []
    for line in e:
        es.append(line)
    return es

# Generate perturber based on given parameters.
def generate_random_perturber_orientation(r_min, ecc, M, other_M, kep, psi, theta, phi):
    perturber = Particle()
    perturber.mass = M
    perturber.radius = np.cbrt(M.in_(units.MSun).number) | units.RSun
    total_mass = M + other_M

    mean_anomaly = np.deg2rad(-180)
    semi = r_min / (1 - ecc)
    kep.initialize_from_elements(mass=total_mass, semi=semi, ecc=ecc, mean_anomaly=mean_anomaly, time=0|units.yr, periastron=r_min)

    position = kep.get_separation_vector()
    velocity = kep.get_velocity_vector()


    kep.advance_to_periastron()
    timescale = kep.get_time()
    print timescale.in_(units.day)
    rotation_matrix = rotation_matrix_from_angles(psi, theta, phi)#preform_EulerRotation()
    position_matrix = np.matrix(([[position[0].number], [position[2].number], [position[1].number]]))
    velocity_matrix = np.matrix(([[velocity[0].number], [velocity[2].number], [velocity[1].number]]))
    position_new = np.dot(rotation_matrix, position_matrix) | position[0].unit
    velocity_new = np.dot(rotation_matrix, velocity_matrix) | velocity[0].unit

    perturber.position = [position_new[0][0], position_new[1][0], position_new[2][0]]
    perturber.velocity = [velocity_new[0][0], velocity_new[1][0], velocity_new[2][0]]
    #print np.rad2deg([psi, theta, phi])
    #print perturber.position.in_(units.AU)

    return perturber, (psi, theta, phi), 2 * timescale

def run(param, mass, name_int, initial_e, percent_increase):
    path = "data/distance/closer"+param+"/"
    if not os.path.exists(path): os.makedirs(path)

    # Do a binary search for closest approach distance
    max_distance = 100000.0 | units.AU
    min_distance = 0.05 | units.AU
    closest_distance = 0 | units.AU

    # Randomly sample about the most effecitve orientation
    psi = np.random.uniform(0, np.pi/20, 1)[0]
    theta = np.random.uniform(np.pi/2, np.pi/20, 1)[0]
    phi = np.random.uniform(0, np.pi/20, 1)[0]

    n = 0
    f = open(path+"distance%05d.txt" % name_int, "w")
    f.write(str(mass)+"\n")
    f.write(str(percent_increase)+"\n")

    # Binary search
    while max_distance - min_distance > (10000000 | units.km):
        bodies = gen_trappist_system(10)
        closest_distance = (max_distance - min_distance) / 2 + min_distance
        print "Checking", closest_distance.in_(units.AU)
        r_min = closest_distance
        v_inf = 3.5 | units.kms

        mu = constants.G * bodies.mass.sum()
        semimajor_axis = - mu / (v_inf * v_inf)
        ecc = 1 - r_min / semimajor_axis

        # Find the right eccentricity
        print ecc
        M = mass | units.MSun
        converter = nbody_system.nbody_to_si(bodies.mass.sum() + M, 2 * (M.number)**(1.0 / 3.0) | units.RSun)
        kep = Kepler(unit_converter=converter)

        perturber, angles, timescale = generate_random_perturber_orientation(r_min, ecc, M, bodies.mass.sum(), kep, psi, theta, phi)
	bodies.add_particle(perturber)
	kep.stop()
        es = get_largest_delta_e(bodies, perturber, timescale, converter, n, param)
	eject_stable = True
        check_for_e = percent_increase != -1
        for ee in es:
	    for e in ee:
	        if np.isnan(e):
		    eject_stable = False
		    break
	if eject_stable and (not check_for_e or initial_e * (1 + percent_increase) > max(es[-1])):
	    max_distance = closest_distance
	    print "is stable!"
	else:
	    print "is not stable!"
	    min_distance = closest_distance
	n += 1
    f.write(str(mass)+"\n")
    f.write(str(closest_distance)+"\n")
    f.close()

if __name__ == "__main__":
    initial = 0.009833817293648425
    percent_increase = 0.1
    run(sys.argv[1], float(sys.argv[2]), 45, initial, percent_increase)
