#------------------------------------------------------------------------------#
#--------------------------------plot_closest.py-------------------------------#
#------------------------------------------------------------------------------#
#--------------------------Created by Nick DeFilippis--------------------------#
#------------------------------------------------------------------------------#
#---This program was created for the purpose of plotting the closest approach--#
#-distance that a stellar perturber needs to take in order to disrupt a system-#
#------------------------------------------------------------------------------#
#---------------------------Date Created: 08/10/2018---------------------------#
#------------------------Date Last Modified: 08/28/2018------------------------#
#------------------------------------------------------------------------------#
from amuse.lab import *
import os
import matplotlib; matplotlib.use('agg')
import numpy as np, matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import interpolate

path = "../../distance_ecc" # path to data files
def parse_file(filename, contains_ecc=False):
    '''
    Parse data file. Data is in the format
    mass
    <required eccentricity change>
    <mass>
    closest approach distance
    <> - means that these fields exist if contains_ecc is true
    '''
    f = open(filename, "r")
    mass = f.readline()
    if mass == "":
        return None
    if contains_ecc:
        ecc = float(f.readline())
	f.readline() # get rid of extra mass line
    distance = f.readline()
    if distance == "":
        return None
    mass = float(mass.split(" ")[0]) | units.MSun
    distance = float(distance.split(" ")[0]) | units.AU
    if contains_ecc:
        return mass, distance, ecc
    else:
        return mass, distance

def generate_mass_map(path, contains_ecc=False):
    '''
    Take in all data files in `path` and generate a dictionary
    of eccentricity changes mapped to lists of (mass, distance)
    pairs.
    '''
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
    '''
    Calculate running average based on old average and count and new data
    '''
    return (x_old * n_old + x_new) / (n_old + 1.0)

def standard_deviation(data, average):
    '''
    Calculate standard deviation based on data and average
    '''
    if len(data) == 1:
       return 0
    variance = sum([(d - average) * (d - average) for d in data]) / float(len(data)-1)
    return np.sqrt(variance)

def plot_distances(mass_dict, hist_data):

    #----------------------
    #|   Begin Drawing    |
    #----------------------

    fig = plt.figure(figsize=(18, 9))
    hist_ax = plt.subplot(313)
    hist_ax.invert_yaxis()
    plt.xscale('log') # Plot masses on a log-scale
    plt.xlim(0.01, 10)
    hist_ax.xaxis.tick_top()
    hist_ax.xaxis.set_ticklabels([])
    plt.ylabel("Mass Distribution", fontsize=22)
    plt.subplots_adjust(hspace=0.25)

    line_ax = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    plt.ylim(0, 2.0)
    plt.xlim(0.01, 10)
    plt.xticks([0.1, 1.0, 10.0])
    for tick in line_ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)

    hist_ax.tick_params(labelsize=20)
    line_ax.tick_params(labelsize=20, pad=10)
    
    masses = []
    distances = []
    mass_list = [(m.number, d.number) for m, d in hist_data[-1]]
    mass_list.sort(key= lambda x: x[0])
    
    bins = np.logspace(np.log10(0.01), np.log10(100.0), num=25, endpoint=True, base=10.0)
    # Create a histogram based on distances
    # Bins are spaced logarithmically
    for mass, distance in mass_list[1:]:
	if mass < bins[len(masses)]:
    	    #d_old, n_old = distances[-1]
    	    masses[-1].append(mass)
    	    distances[-1].append(distance)# = (running_average(d_old, n_old, distance), n_old+1)
    	else:
    	    masses.append([mass])
    	    distances.append([distance])
    
    average_masses = [sum(m) / len(m) for m in masses] # keep track of various statistics for each bin
    min_masses = [min(m) for m in masses]
    max_masses = [max(m) for m in masses]
    distances_histogram = [len(d) for d in distances]
    distances_histogram = [float(d) / sum(distances_histogram) for d in distances_histogram]
    
    xnew = np.logspace(np.log2(min(average_masses)+0.001), np.log2(max(average_masses)-0.01), base=2.0, num=1000)
    f = interpolate.interp1d(average_masses, distances_histogram, kind='cubic')
    ynew = f(xnew) # smoothing process

    hist_ax.hist(xnew, bins=xnew, weights = ynew, color='#07294d')
    plt.xlabel("Mass [MSun]", fontsize=22, labelpad=-20)
    plt.ylabel("Closest Approach [AU]", fontsize=22)

    interval = mpatches.Patch(facecolor=(151/255.0, 27/255.0, 141/255.0, 0.2), edgecolor='#971b2f', linestyle='--', label='95% confidence interval')

    colors = ['#b7bf10', '#de7c00', '#971b2f']
    lines = []
    labels = ['10% increase', '100% increase', 'Ejection']
    c_i = 0
    for eccentricity_change in ([0.1, 1.0, -1]):
        mass_list = mass_dict[eccentricity_change]
	mass_list.sort(key=lambda x: x[0]) # sort by mass
	masses = []
        distances = []
        max_distances = []
        min_distances = []

        mass_list = [(m.number, d.number) for m, d in mass_list]

        bins = np.logspace(np.log10(0.01), np.log10(100.0), num=20, endpoint=True, base=10.0)

        for mass, distance in mass_list[1:]:
    	    if mass < bins[len(masses)]:
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
        print eccentricity_change, distances_histogram
        min_distances = [min(d) for d in distances]
        max_distances = [max(d) for d in distances]
        average_distances = [sum(d) / len(d) for d in distances]

	# Chebyshev Inequality
        upper_std = [av + 3 * standard_deviation(d, av) for d, av in zip(distances, average_distances)]
        lower_std = [av - 3 * standard_deviation(d, av) for d, av in zip(distances, average_distances)]

        xnew = np.logspace(np.log2(min(average_masses)+0.001), np.log2(max(average_masses)-0.01), base=2.0, num=250)
	f = interpolate.interp1d(average_masses, average_distances, kind='cubic')	
	ynew = f(xnew)
        #mss = [m for m, d in mass_list]
        #dss = [d for m, d in mass_list]
        #plt.scatter(mss, dss, zorder=1)
        line1 = line_ax.semilogx(xnew, ynew, linewidth=4, color=colors[c_i], label=labels[c_i])
        lines.append(line1[0])
	f = interpolate.interp1d(average_masses, upper_std, kind='cubic')
	upper_new = f(xnew)
        line2 = plt.semilogx(xnew, upper_new, linestyle='--', color=colors[c_i])
	f = interpolate.interp1d(average_masses, lower_std, kind='cubic')
	lower_new = f(xnew)
        line3 = plt.semilogx(xnew, lower_new, linestyle='--', color=colors[c_i])

        xs = np.concatenate((xnew, xnew[::-1]))
        ys = np.concatenate((upper_new, lower_new[::-1]))
	plt.fill(xs, ys, colors[c_i], alpha=0.2)
        c_i += 1

    fig.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.05, hspace=0.25)
    lines.append(interval)
    leg = plt.legend(handles=lines, loc='upper left', bbox_to_anchor=(0,1), fontsize=20, handlelength=2)#loc = (1.5, 1))#, loc = 'upper right')
    leg.get_frame().set_edgecolor('white')

        
    plt.suptitle('Mass versus Closest Approach Distance', fontsize=28)
    plt.savefig("mass_v_distance.png", dpi=1000)

if __name__ == "__main__":
    data = generate_mass_map(path, contains_ecc=True)
    hist_data = generate_mass_map(path, contains_ecc=True)
    #data = generate_mass_map('/home/draco/ndefilippis/distance_ecc_lowe', contains_ecc=True)
    data2 = generate_mass_map("/home/draco/ndefilippis/distance_ecc_list", contains_ecc=True)
    data1 = generate_mass_map("/home/draco/ndefilippis/distances")
    data3 = generate_mass_map("/home/draco/ndefilippis/distance_ecc_low", contains_ecc=True)
    for key in data2:
         if key in data:
            data[key] = data[key] + data2[key]
         else:
    	    data[key] = data2[key]
    data[0.001] = data3[0.001]
    data[-1] = data[-1] + data1[-1]
    print data.keys()
    print [len(data[k]) for k in data.keys()]
    plot_distances(data, hist_data)
