### TRAPPIST-sim 

## This package was developed to investigate the effects of close stellar encounters on the TRAPPIST-1 system, particularly because of its relatively small size.

## Pre-requisites

This package uses Python 2.7, the [AMUSE](https://github.com/amusecode/amuse) framework, and the [Tycho](https://github.com/JPGlaser/Tycho) community code.

## Package Overview ##

* sim_secularmultiples.py *

This package contains several important functions to run the SecularMultiples integrator on a bound system. 
The main function in this library is `run_secular_multiples(bodies, t_max, n_steps)`
