#------------------------------------------------------------------------------#
#------------------------------plot_encounters.py------------------------------#
#------------------------------------------------------------------------------#
#--------------------------Created by Mark Giovinazzi--------------------------#
#------------------------------------------------------------------------------#
#-This program was created for the purpose of loading in saved encounter files-#
#-to plot our encounters in a more general sense and will also generate movies-#
#------------------------------------------------------------------------------#
#---------------------------Date Created: 12/??/2017---------------------------#
#------------------------Date Last Modified: 05/02/2018------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#-------------------------------Import Libraries-------------------------------#
#------------------------------------------------------------------------------#

import matplotlib; matplotlib.use('agg')
import numpy as np, matplotlib.pyplot as plt
import os, sys, glob, time, pickle
from amuse.lab import *
from amuse import datamodel
from amuse.community.kepler.interface import kepler
#from sim_encounters import get_stars, get_planets
import matplotlib.patches as mpatches

def parse_data_from_file(filename):
    f = open(filename)
    angles = f.readline()
    eccs = f.readline()
    psi, theta, phi = map(np.rad2deg, map(float, angles[1:-2].split(", ")))
    eccs = map(float, eccs[1:-2].split(", "))
    return psi, theta, phi, eccs

def get_delta_e_metric(eccs):
    return sum(eccs) / len(eccs)

def running_average(x_old, n_old, x_new, weight):
    n_new = n_old + weight
    return (x_old * n_old + weight * x_new) / n_new, n_new

def gen_heatmap(xs, ys, ws, min, max, bins):
    arr = np.empty((bins, bins, 2))
    dt = float(max - min) / bins
    
    for x, y, w in zip(xs, ys, ws):
        x_b = (x - min) / dt
	y_b = (y - min) / dt
        
	j = int(x_b)
	i = int(y_b)
	xf = x_b - j
	yf = y_b - i

	if j+1 >= bins:
	    xf = 1
	if i+1 >= bins:
	    yf = 1
	
	x_old, n_old = arr[i, j]
	arr[i, j] = running_average(x_old, n_old, w, xf*yf)

	if j+1 < bins and i+1 < bins:
	    x_old, n_old = arr[i+1, j+1]
	    arr[i+1,j+1] = running_average(x_old, n_old, w, (1-xf) * (1-yf))
	if j+1 < bins:
	    x_old, n_old = arr[i, j+1]
	    arr[i, j+1] = running_average(x_old, n_old, w, (1-xf) * yf)
	if i+1 < bins:
	    x_old, n_old = arr[i+1, j]
	    arr[i+1, j] = running_average(x_old, n_old, w, xf * (1-yf))
    
    # deal with zero values
    for i in range(len(arr)):
        for j in range(len(arr[i])):
	    if arr[i, j, 1] == 0:
	        arr[i, j, 0] = 0.01
		arr[i, j, 1] = 1
    return arr[:,:,0]

#------------------------------------------------------------------------------#
#-The following function will handle all of this project's visualization tools-#
#------------------------------------------------------------------------------#
def generate_visuals(psis, thetas, phis, eccs, path, name, make_movie = False, make_map = True):
    bins = 35
    
    ecc_min = min(eccs)
    ecc_max = max(eccs)
    
    plot_max = 90.0
    plot_min = -90.0

    # If the user wants to make trail maps in their run, the following code will be executed
    if make_map == True:

        plt.figure(figsize = (5.5, 5))
        #print xs[0], zs[0], cs[0], ss[0]
        # x by z plot (top left corner)
        plt.subplot(221, aspect = 'equal')
        ax = plt.gca()
        ax.set_facecolor('black')
        ax.set_xticklabels([])
        #print range(len(xs))
        #print [(xs[i], zs[i], cs[i], ss[i]) for i in range(len(xs))]
        #[plt.scatter(float(xs[i]), float(zs[i]), c = cs[i], s = int(ss[i])) for i in range(len(xs))]
        #print '---------------', xs, zs
        #print xs[:12]
	heatmap = gen_heatmap(psis, phis, eccs, plot_min, plot_max, bins)
	extent = [-90.0, 90.0, -90.0, 90.0]

        plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max)
        plt.ylabel('$\phi [deg]$')
        plt.xlim(plot_min, plot_max)
        plt.ylim(plot_min, plot_max)

        # x by y plot (lower left corner)
        plt.subplot(223, aspect = 'equal')
        ax = plt.gca()
        ax.set_facecolor('black')
	heatmap = gen_heatmap(psis, thetas, eccs, plot_min, plot_max, bins)
        plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max)
        plt.xlabel('$\psi [deg]$')
        plt.ylabel(r"$\theta [deg]$")
        plt.xlim(plot_min, plot_max)
        plt.ylim(plot_min, plot_max)        
        #ax.set_xticklabels([min, max])
        #ax.set_yticklabels([min, max])
        #ax.set_xicks([min, max])
        #ax.set_yicks([min, max])
        # z by y plot (lower right corner)
        plt.subplot(224, aspect = 'equal')
        ax = plt.gca()
        ax.set_facecolor('black')
        ax.set_yticklabels([])
        #[plt.scatter(float(zs[i]), float(ys[i]), c = cs[i], s = int(ss[i])) for i in range(len(xs))]

        heatmap = gen_heatmap(phis, thetas, eccs, plot_min, plot_max, bins)
	plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max)
        plt.xlabel('$\phi [deg]$')
        plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1475, right=0.95, hspace=0.25, wspace=0.01)
        #plt.subplots_adjust(wspace = 0.001) # shorten the width between left and right since there aren't tick marks
        plt.xlim(plot_min, plot_max)
        plt.ylim(plot_min, plot_max)
	plt.colorbar()
        #leg = plt.colorbar(loc='upper right', bbox_to_anchor=(0.9, 2))#loc = (1.5, 1))#, loc = 'upper right')
        #leg.get_frame().set_edgecolor('white')

        # make some plot_radius function to determine value for s for stars and planets...
        # stars and planets should look noticeably different in size, but stars/stars and planets/planets should be subtle
        plt.suptitle("$\Delta\epsilon$ versus Orientation")
        print 'saving///' # just for the user's sake..
        plt.savefig('angle_maps/' + str(name) + '.pdf', dpi = 1000)
        plt.clf()


#--------------------------------------
# Here is where we begin the program...
#--------------------------------------

# set up some conditionals to determine what we're going to use this porgram for
gen_vis = True   # if `gen_vis` is True, you will get a trail map for each close encounter 
gen_dict = False # if `gen_dict` is True, you will get a dictionary of orbital elements for each close encounter
gen_movie = True # if `gen_movie` is True, you will get snapshots of each close encounter (which you can later stitch together to make a movie)

# Define the path for all of our hdf5 files
par_str = sys.argv[1]
data_path = par_str #+ par_str # where the data is going to be pulled from
data_dirs = os.listdir(data_path)

nan_counter = 0


psis = []
thetas = []
phis = []
eccs = []
# Each data file has its own subdirectory; let's chalk through them in this loop 
for data_dir in data_dirs:

    sample_paths = os.listdir(os.path.join(data_path, data_dir))

    for sample_path in sample_paths:
        fin_data_path = os.path.join(data_path, data_dir, sample_path)
	# If we compute how many stars and planets there are up top, we can ignore doing so in the next loop       
	psi, theta, phi, ecc_arr = parse_data_from_file(fin_data_path)
	psis.append(psi)
	thetas.append(theta)
	phis.append(phi)
	eccs.append(get_delta_e_metric(ecc_arr))
generate_visuals(psis, thetas, phis, eccs, path = None, name = str(fin_data_path).replace('/', '_').replace('__', '_').replace('_.', '.'), make_movie = gen_movie, make_map = True)

print 'There are currently', nan_counter, 'files that have turned into nans.'
