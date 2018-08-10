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
    X = np.array([1, 0, 0])
    Y = np.array([0, 1, 0])
    X_prime = np.dot(rotation_matrix, X)
    Y_prime = np.dot(rotation_matrix, Y)
    psi = np.arcsin(X_prime[1] / np.sqrt(1 - X_prime[2] ** 2))
    theta = np.arcsin(-X_prime[2])
    phi = np.arcsin(Y_prime[2] / np.sqrt(1 - X_prime[2] ** 2))
    return psi, theta, phi
   
def get_largest_delta_e(bodies, perturber, timescale, converter, n, param):
    path = "ce_directory/angles"+param+"/TRAPPIST/1"
    if not os.path.exists(path): os.makedirs(path)
    _, new_bodies = run_collision(bodies, timescale, 10|units.hour, str(n), path, converter=converter)
    
    #print new_bodies
    t, a, e, i, _ = run_secular_multiples(new_bodies-perturber, 10000|units.yr, 50)

    delta_es = []
    for line in e:
        delta_es.append(max(line) - min(line))
    return delta_es

def generate_random_perturber_orientation(r_min, ecc, M, other_M, kep):
    perturber = Particle()
    perturber.mass = M
    perturber.radius = np.cbrt(M.in_(units.MSun).number) | units.RSun
    total_mass = M + other_M

    mean_anomaly = np.deg2rad(-135)
    semi = r_min / (1 - ecc)
    kep.initialize_from_elements(mass=total_mass, semi=semi, ecc=ecc, mean_anomaly=mean_anomaly, time=0|units.yr, periastron=r_min)
    
    position = kep.get_separation_vector()
    velocity = kep.get_velocity_vector()


    kep.advance_to_periastron()
    timescale = kep.get_time()

    rotation_matrix = preform_EulerRotation()
    psi, theta, phi = extract_EulerAngles(rotation_matrix)

    position_matrix = np.matrix(([[position[0].number], [position[2].number], [position[1].number]]))
    velocity_matrix = np.matrix(([[velocity[0].number], [velocity[2].number], [velocity[1].number]]))
    
    position_new = np.dot(rotation_matrix, position_matrix) | position[0].unit
    velocity_new = np.dot(rotation_matrix, velocity_matrix) | velocity[0].unit

    perturber.position = [position_new[0][0], position_new[1][0], position_new[2][0]]
    perturber.velocity = [velocity_new[0][0], velocity_new[1][0], velocity_new[2][0]]
    return perturber, (psi, theta, phi), 2 * timescale

def run(param):
    path = "rebuild_angles/angles"+param+"/"
    if not os.path.exists(path): os.makedirs(path)

    for n in range(100):
        bodies = gen_trappist_system(10)
   
        r_min = 1.0 | units.AU
        ecc = 1.1
        M = 4.0 | units.MSun

        converter = nbody_system.nbody_to_si(bodies.mass.sum() + M, 2 * (M.number)**(1.0 / 3.0) | units.RSun)
        kep = Kepler(unit_converter=converter)

        perturber, angles, timescale = generate_random_perturber_orientation(r_min, ecc, M, bodies.mass.sum(), kep)
        bodies.add_particle(perturber)
        kep.stop() 
        f = open(path+"angle%03d.txt" % n, "w")
	f.write(str(angles) + "\n")
        f.write(str(get_largest_delta_e(bodies, perturber, timescale, converter, n, param)) + "\n")
	f.close()

if __name__ == "__main__":
    run(sys.argv[1])