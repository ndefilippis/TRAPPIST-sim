from amuse.lab import *
import os
import matplotlib; matplotlib.use('agg')
import numpy as np, matplotlib.pyplot as plt
import matplotlib.patches as mpatches

path = "distances"
path2 = "distance_ecc_list"
def parse_file(filename, contains_ecc=False):
    f = open(filename, "r")
    mass = f.readline()
    if mass == "":
        return None
    if contains_ecc:
        ecc = float(f.readline())
	f.readline() # get rid of extra mass line
    distance = f.readline()
    mass = float(mass) | units.MSun
    distance = float(distance.split(" ")[0]) | units.AU
    if contains_ecc:
        return mass, distance, ecc
    else:
        return mass, distance

def generate_mass_map(path, contains_ecc=False):
    ret = {}
    ret[-1] = []
    mass_dirs = os.listdir(path)
    for mass_dir in mass_dirs:
        mass_path = os.path.join(path, mass_dir)
	file_paths = os.listdir(mass_path)
	for file_path in file_paths:
	    filename = os.path.join(path, mass_dir, file_path)
	    data = parse_file(filename, contains_ecc=contains_ecc)
	    if data != None:
	        if contains_ecc:
		    mass, distance, ecc = data
		    if ecc in ret:
		        ret[ecc].append((mass, distance))
		    else:
		        ret[ecc] = [(mass, distance)]
		else:
                    ret[-1].append(data)
    return ret

def running_average(x_old, n_old, x_new):
    return (x_old * n_old + x_new) / (n_old + 1.0)

def standard_deviation(data, average):
    variance = sum([(d - average) * (d - average) for d in data]) / float(len(data)-1)
    return np.sqrt(variance)

def plot_distances(mass_list):
    mass_list.sort(key=lambda x: x[0]) # sort by mass
    masses = []
    distances = []
    max_distances = []
    min_distances = []

    mass_list = [(m.number, d.number) for m, d in mass_list]
    
    masses.append([mass_list[0][0]])
    distances.append([mass_list[0][1]])

    for mass, distance in mass_list[1:]:
        old_mass = masses[-1][0]
	if (mass - old_mass) / old_mass < 0.4:
	    #d_old, n_old = distances[-1]
	    masses[-1].append(mass)
	    distances[-1].append(distance)# = (running_average(d_old, n_old, distance), n_old+1)
	else:
	    masses.append([mass])
	    distances.append([distance])

    average_masses = [sum(m) / len(m) for m in masses]
    min_masses = [min(m) for m in masses]
    max_masses = [max(m) for m in masses]

    distances_histogram = [len(d) for d in distances]
    print distances_histogram
    print average_masses
    min_distances = [min(d) for d in distances]
    max_distances = [max(d) for d in distances]
    average_distances = [sum(d) / len(d) for d in distances]
    upper_std = [av + 2 * standard_deviation(d, av) for d, av in zip(distances, average_distances)]
    lower_std = [av - 2 * standard_deviation(d, av) for d, av in zip(distances, average_distances)]
    #----------------------
    #|   Begin Drawing    |
    #----------------------
    
    plt.figure(figsize=(20, 10))
    
    ax = plt.subplot(313)
    #plt.title("Mass Histogram") 
    ax.invert_yaxis()

    plt.hist(average_masses, bins=min_masses+[max_masses[-1]], weights = distances_histogram, color='#07294d')
    plt.xscale('log')
    plt.xlim(0.01, 100)
    ax.xaxis.tick_top()
    ax.xaxis.set_ticklabels([])
    plt.ylabel("Number of Samples", fontsize=15)
    plt.subplots_adjust(hspace=0.25)

    ax = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
        
    #mss = [m for m, d in mass_list]
    #dss = [d for m, d in mass_list]
    #plt.scatter(mss, dss, zorder=1)
    line1 = plt.semilogx(average_masses, average_distances, linewidth=4, color='#07294d')
    line2 = plt.semilogx(average_masses, upper_std, linestyle='--', color='#07294d')
    line3 = plt.semilogx(average_masses, lower_std, linestyle='--', color='#07294d')
    plt.xlim(0.01, 100)
    xs = average_masses + average_masses[::-1]
    ys = upper_std + lower_std[::-1]
    plt.fill(xs, ys, '#ffc600')

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)

    plt.xlabel("Mass [MSun]", fontsize=15)
    plt.ylabel("Closest Approach [AU]", fontsize=15)

    interval = mpatches.Patch(facecolor='#ffc600', edgecolor='#07294d', linestyle='--')
  
    leg = plt.legend((line1[0], interval), ("Closest Approach", "95% confidence interval"), loc='upper left', bbox_to_anchor=(0,1), fontsize=15, handlelength=2)#loc = (1.5, 1))#, loc = 'upper right')
    leg.get_frame().set_edgecolor('white')



    plt.suptitle('Mass versus Closest Approach Distance', fontsize=22)
    plt.savefig("mass_v_distance.png", dpi=1000)
    
if __name__ == "__main__":
    data_eject = generate_mass_map(path)
    data = generate_mass_map(path2, contains_ecc=True) 
    data[-1] = data_eject[-1]
    print data.keys()
    print [len(data[k]) for k in data.keys()]
    plot_distances(data)
