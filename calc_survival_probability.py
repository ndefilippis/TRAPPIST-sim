from plot_closest import generate_mass_map
import numpy as np
from amuse.lab import units, constants

def calc_survival_rate(distance_list):
    distance_list.sort(key=lambda tuple: tuple[0])
    
    masses = [x for x, _ in distance_list]
    distances = [y for _, y in distance_list]
    bins = np.logspace(np.log2(0.01), np.log2(15.0), endpoint=True, base=2.0)
    
    distance_hist = [(0 | units.AU,0)]
    masses_hist = [(0 | units.MSun,0)]
    for m, d in zip(masses, distances):
        if m > (bins[len(distance_hist)] | units.MSun):
	    distance_hist.append((d, 1))
	    masses_hist.append((m, 1))
	else:
	    x_curr, n_curr = distance_hist[-1]
	    new_n = n_curr + 1.0
	    distance_hist[-1] = ((x_curr * n_curr + d) / new_n, new_n)
	    x_curr, n_curr = masses_hist[-1]
	    new_n = n_curr + 1.0
	    masses_hist[-1] = ((x_curr * n_curr + m) / new_n, new_n)

    distance_hist = [d for d, _ in distance_hist]
    masses_hist = [m for m, _ in masses_hist]
    hist = np.histogram([m.number for m in masses], bins=bins, density=True)

    survival_prob = 0 | (units.s**-1)
    for i in range(len(hist[0])-1):
	closest_approach = distance_hist[i]
	average_mass = masses_hist[i]
	v_rms = 0.5 | units.kms
	TRAPPIST_mass = (0.089 | units.MSun) + (4 | units.MSun)#average_mass
	#print closest_approach
	factor = (constants.G * TRAPPIST_mass) / (closest_approach * v_rms**2)
	#print factor
	impact_parameter_2 = closest_approach ** 2 * (1 + factor)
	dynamical_time = 100 | units.Myr
        #print dynamical_time.in_(units.Myr)
	dm = hist[1][i+1] - hist[1][i]
        number_density = hist[0][i] * dm * 10.5 / (1 | units.lightyear)**3
	survival_prob += number_density * v_rms * impact_parameter_2
    return dynamical_time * survival_prob

if __name__ == "__main__":
    data_path = "/home/draco/ndefilippis/distance_ecc"
    data = generate_mass_map(data_path, contains_ecc=True)
    print calc_survival_rate(data[0.01])
