### TRAPPIST-sim 

This package was developed to investigate the effects of close stellar encounters on the TRAPPIST-1 system, particularly because of its relatively small size.

## Pre-requisites

This package uses Python 2.7, the [AMUSE](https://github.com/amusecode/amuse) framework, and the [Tycho](https://github.com/JPGlaser/Tycho) community code.

## Package Overview

# Packages 

**sim_secularmultiples.py**

This package contains several important functions to run the SecularMultiples integrator on a bound system. 

`run_secular_multiples(bodies, t_max, n_steps)`

`bodies` - Particle set of bodies to be simulated (Particles)

`t_max` - the amount of time to simulate for (in amuse.units.time)

`n_step` - the number of steps to get data for (int)

Returns a tuple of (timesteps, semimajor axes, eccentricities, inclinations, bodies)

`timesteps` is a list of timesteps based on the number of steps and the maximum time

`semimajor_axes` is a list of lists, each representing the semimajor axis of each binary for each timestep

`eccentricities` is a list of lists, each representing the eccentricities of each binary for each timestep

`inclinations` is a list of lists, each representing the inclinations of each binary for each timestep

`bodies` is a Particle Set of the system after the simulation is over

**trappist_system.py**

This package is responsible for generating an accurate representation of the TRAPPIST-1 system, based on available information.

`gen_trappist_system(seed)`

`seed` - a random seed for generating the randomly based data (at this point, only the longitude of ascending node) (int)

Returns a ParticleSet of the TRAPPIST-1 system.

**angles_test.py**

Randomly samples an orientation for a perturber and finds the maximum eccentricity of the outermost planet

**closest_all_test.py**

Simulates an encounter at a particular mass and finds the closest appraoch distance to cause a certain percent change in 
eccentricity

Usage: `python closest_all_test.py folder_id solar_mass percent_change`

Pass in -1 for `percent_change` to find the ejection limit, i.e. the point at which the outermost planet gets ejected from the system

# Multithreading

This package also contains several scripts to pool running several runs at once. All these scripts begin with `pool_`. To
run these, the syntax is usually

python pool_*.py start_folder num_samples

# Generating Graphics

The plot subfolder has several scripts for generating matplotlib graphics based on simulated data
