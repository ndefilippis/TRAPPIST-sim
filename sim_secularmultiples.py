from amuse.lab import *
from amuse.community.secularmultiple.interface import SecularMultiple
from amuse.ext.orbital_elements import orbital_elements_for_rel_posvel_arrays
from amuse.community.kepler.interface import kepler
import numpy as np
import os

# Functions to separate stars from planets - 
# or use Mark's code in sim_encounters
def get_stars(bodies):
    limiting_mass_for_planets = 13 | units.MJupiter
    
    stars = bodies[bodies.mass > limiting_mass_for_planets]
    return stars

def get_planets(bodies):
    limiting_mass_for_planets = 13 | units.MJupiter
    
    planets = bodies[bodies.mass <= limiting_mass_for_planets]
    planets = planets[planets.mass > 0 | units.kg]
    return planets

# Returns the relative distance between to bodies
def rel_distance(body1, body2):
    return np.linalg.norm(body1.position - body2.position)


# Returns a list of binaries from orbiting bodies.
# Pulls out stars and categorizes planets that are
# orbiting each one.
def make_binaries(bodies, converter = None):
    
    # Use kepler to verify orbital parameters have an elliptical orbit
    # kep = Kepler(unit_converter = converter, redirection = 'none')
    # kep.initialize_code()
    
    stars, planets = get_stars(bodies), get_planets(bodies)
    num_stars, num_planets = len(stars), len(planets)
    
    root_binary = None
    
    # Each planetary system can be converted to a hierarchical binary system,
    # This dictionary is indexed by the star and contains the binary and a particle
    # set as the value
    
    star_system_binaries = {}
    
    for star in stars:
        star.is_binary = False
        star.child1 = None
        star.child2 = None
        star_system_binaries[star] = star
    
    # Associate each planet with its closest star
    # Note: Does not take into account planets orbiting binary stars
    converter = nbody_system.nbody_to_si(1|units.MEarth, 1|units.AU)
    #kep = Kepler(unit_converter=converter)
    for planet in planets:
        planet.is_binary = False
        planet.child1 = None
        planet.child2 = None
        # find planet with closest distance to a star
        # If this method becomes a bottleneck, look
        # into R*-trees for a better data structure
        # https://en.wikipedia.org/wiki/R-tree
        
        min_distance = np.inf | units.m
        min_star = None
        for star in stars:
            
            #TODO: check if relative velocity * position ~= 0,
            # or, check eccentricity of the orbit
	    #rel_pos = star.position - planet.position
	    #rel_vel = star.velocity - planet.velocity
	    #total_mass = star.mass + planet.mass

	    #kep.initialize_from_dyn(total_mass, rel_pos[0], rel_pos[1], rel_pos[2], rel_vel[0], rel_vel[1], rel_vel[2])
	    #a, e = kep.get_elements()
            d = rel_distance(planet, star)
	    #if e > 1:
	    #    "Planet of mass", planet.mass, " in hyperbolic orbit with star of mass", star.mass, "(ecc=",e,")"
            if d < min_distance:# and e < 1:
                min_distance = d
                min_star = star
	    
        if min_star != None:
            # Update binary system with next closest plant
            old_binary = star_system_binaries[min_star]
        
            root_binary = Particle()
            root_binary = bodies.add_particle(root_binary)
            root_binary.is_binary = True
            root_binary.child1 = old_binary
            root_binary.child2 = planet
             
            star_system_binaries[min_star] = root_binary
    #kep.stop()  
    # Put planetary systems in a heirarchical binary configuration
    
    # Easy if we only have one or two stars
    if len(stars) == 1:
        return bodies[bodies.is_binary], root_binary
    
    if len(stars) == 2:
        # Make pairing between the only two stars
        root_binary = Particle()
        root_binary = bodies.add_particle(root_binary)
        
        root_binary.is_binary = True
        star1, star2 = star_system_binaries.keys()
        root_binary.child1 = star_system_binaries[star1]
        root_binary.child2 = star_system_binaries[star2]
        return bodies[bodies.is_binary], root_binary
    
    # This part is a bit more complicated and may be unnecessary
    
    # Sort by closest distance between stars
    closest = stars.nearest_neighbour()
    for star, neighbour in zip(stars, closest):
        star.pairing_distance = rel_distance(star, neighbour)
        star.neighbour = neighbour

    stars_by_closest = stars.sorted_by_attribute("pairing_distance")
    
    # The hierarchical binary each star is associated with
    # We basically build the tree from the ground up
    star_binaries = {}
    for star in stars_by_closest:
        star_binaries[star] = star
        
    for star in stars_by_closest:
        root_binary = Particle()
        root_binary = bodies.add_particle(root_binary)
        
        root_binary.is_binary = True
        root_binary.child1 = star_binaries[star]
        root_binary.child2 = star_binaries[star.neighbour]
        
        #Update binaries
        star_binaries[star] = root_binary
        star_binaries[star.neighbour] = root_binary
    
    return bodies[bodies.is_binary], root_binary
   
def print_tree(binary): 
    if binary.is_binary:
        print_tree(binary.child1)
        print "---",
    elif not binary.is_binary:
        print binary.mass
        return
    print "---", 
    print_tree(binary.child2)
 


# evolve model using Secular Multiples. Assumes that hiearchy of binaries
# has already been set up
def evolve_secular_multiples(bodies, binaries, end_time, n_steps, root_binary):
    
    dt = end_time / n_steps
    time = 0 | units.yr
    
    # Initialize secular mutliples
    code = SecularMultiple()
    code.particles.add_particles(bodies)
    
    # Set up channels
    channel_from_bodies_to_code = bodies.new_channel_to(code.particles)
    channel_from_code_to_bodies = code.particles.new_channel_to(bodies)
    channel_from_bodies_to_code.copy()
    
    
    N_binaries = len(binaries)
    
    # Create arrays for output
    print_times_Myr = [] | units.Myr
    print_semimajor_axis_AU = [[] | units.AU for _ in range(N_binaries)]
    print_eccentricity = [[] for _ in range(N_binaries)]
    print_inclination = [[] for _ in range(N_binaries)]
    
    # Append initial conditons
    print_times_Myr.append(time.in_(units.Myr))
    for index in range(N_binaries):
        print_semimajor_axis_AU[index].append(binaries[index].semimajor_axis)
        print_eccentricity[index].append(binaries[index].eccentricity)
        print_inclination[index].append(binaries[index].inclination)
    
    # Do the actual simulation
    while time <= end_time:

        #if is_stable(bodies):
        #    break

        time += dt
        #print time.in_(units.Myr)
        code.evolve_model(time)
        
        channel_from_code_to_bodies.copy()
        
        # write to output arrays
        print_times_Myr.append(time.in_(units.Myr))
        for index in range(N_binaries):
            print_semimajor_axis_AU[index].append(binaries[index].semimajor_axis)
            print_eccentricity[index].append(binaries[index].eccentricity)
            print_inclination[index].append(binaries[index].inclination)
    code.stop()        
    return print_times_Myr, print_semimajor_axis_AU, print_eccentricity, print_inclination, bodies

# determines whether a planetary system is completely stable or not
# checks if any body is within a Hill Sphere radius of another's orbit

def get_hill_radius(binary):
    min_mass = min(binary.child1.mass, binary.child2.mass)
    max_mass = max(binary.child1.mass, binary.child2.mass)
    hill_radius = binary.semimajor_axis * (1 - binary.eccentricity) * np.cbrt(min_mass / (3 * max_mass))
    return hill_radius

def is_stable(bodies):
    binaries = bodies[bodies.is_binary]
    binaries = binaries.sorted_by_attribute("semimajor_axis")
    for index in range(len(binaries)-1):
        inner_binary = binaries[index]
        outer_binary = binaries[index+1]
        
        inner_hill = get_hill_radius(inner_binary)
        outer_hill = get_hill_radius(outer_binary)
        
        print outer_hill.in_(units.AU), outer_binary.mass
        
        a_inner = inner_binary.semimajor_axis
        a_outer = outer_binary.semimajor_axis
        e_inner = inner_binary.eccentricity
        e_outer = inner_binary.eccentricity
        
        max_inner_distance = a_inner * (1 + e_inner)
        min_outer_distance = a_outer * (1 - e_inner)
        
        if max_inner_distance + inner_hill > min_outer_distance - outer_hill:
            return False
    return True


# parse planetary systems after running Mark's code
def get_planetary_systems(bodies):
    # Make a dictionary of all star ids, populate list with star
    stars = bodies[bodies.host_star == 0]
    
    # link all associated bodies into one solar system set
    planetary_systems = []
    for star in stars:
        mask = np.logical_or(bodies.host_star == star.id, bodies.id == star.id)
        planetary_systems.append(bodies[mask])
    
    # filter out lone stars
    return [system for system in planetary_systems if len(system) > 1]
    
    
def set_orbital_parameters(binary, G, binary_dict={}):
    
    # Get the mass, position, and velocity from a particle or binary
    def get_properties_from_body(body):
        if body.is_binary: # assumes that body is in the binary map
            mass = binary_dict[body.key]["mass"]
            position = binary_dict[body.key]["position"]
            velocity = binary_dict[body.key]["velocity"]
        else:
            mass = body.mass
            position = body.position
            velocity = body.velocity
        return mass, position, velocity
            
    body1 = binary.child1
    body2 = binary.child2
    
    # Make sure children get their orbital properties set first
    if body1.is_binary and not body1.key in binary_dict:
        set_orbital_parameters(body1, G)
    if body2.is_binary and not body2.key in binary_dict:
        set_orbital_parameters(body2, G)
    
    mass_1, position_1, velocity_1 = get_properties_from_body(body1)
    mass_2, position_2, velocity_2 = get_properties_from_body(body2)
    total_mass = mass_1 + mass_2
    center_of_mass = (position_1 * mass_1 + position_2 * mass_2) / total_mass
    center_of_mass_velocity = (velocity_1 * mass_1 + velocity_2 * mass_2) / total_mass
    
    relative_position = position_2 - position_1
    relative_velocity = velocity_2 - velocity_1
    
    orbital_parameters = orbital_elements_for_rel_posvel_arrays(relative_position, relative_velocity, total_mass, G)
    a, e, _, i, Omega, omega = [i[0] for i in orbital_parameters] 
    #a, e, i, o, u = get_orbital_parameters_from_state_vectors(relative_position, relative_velocity, G * total_mass)
    binary.semimajor_axis = a
    binary.eccentricity = e
    binary.inclination = np.deg2rad(i)
    binary.mass = total_mass
    binary.argument_of_pericenter = np.deg2rad(omega)
    binary.longitude_of_ascending_node = np.deg2rad(Omega)
    
    binary_dict[binary.key] = {}
    binary_dict[binary.key]["mass"] = total_mass
    binary_dict[binary.key]["position"] = center_of_mass
    binary_dict[binary.key]["velocity"] = center_of_mass_velocity

# extract all necessary information out of a list of particles
def convert_solar_system_to_particle_set(bodies):
    num_particles = len(bodies)
    particles = Particles(num_particles)
    stars = bodies[bodies.host_star == 0]
    for index in range(num_particles):
        particle = particles[index]
        particle.mass = bodies[index].mass
        particle.position = bodies[index].position
        particle.velocity = bodies[index].velocity
        
        # convert relative position to absolute position
        if bodies[index].type == "planet":
            for star in stars:
                if bodies[index].host_star == star.id:
                    particle.position += star.position 
    return particles

def run_secular_multiples(bodies, t_max, n_steps):
    bodies = bodies.copy()
    binaries, root_binary = make_binaries(bodies)
    # print_tree(root_binary)
    set_orbital_parameters(root_binary, constants.G)
    #print bodies
    return evolve_secular_multiples(bodies, binaries, t_max, n_steps, root_binary)


# Main runner function. Pass in bodies and output
# data from secular multiples
def run_secular_multiples_from_smalln(data_path, bodies, t_max, n_steps):
    print "Running Secular Multiples"
    bodies = bodies.copy()
    systems = get_planetary_systems(bodies)
    results = []
    for system in systems:
        planetary_bodies = convert_solar_system_to_particle_set(system)
        print planetary_bodies
        results.append(run_secular_multiples(planetary_bodies, t_max, n_steps))
    filename = os.path.join(data_path, "secular_run")
    if not os.path.exists(data_path): os.mkdir(data_path)

    f = open(filename, "w")
    f.write(str(results))
    f.close()  
    return bodies
