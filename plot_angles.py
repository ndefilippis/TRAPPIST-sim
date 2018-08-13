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
    if angles == "":
        print filename, "is empty"
	return None
    eccs_all = f.readline()
    psi, theta, phi = map(np.rad2deg, map(float, angles[1:-2].split(", ")))
    eccs_all = eccs_all[2:-3].split("], [")
    eccs = []
    for ecc_planet in eccs_all:
        eccs.append(map(float, ecc_planet.split(", ")))
    return psi, theta, phi, eccs

def get_delta_e_metric(eccs):
    delta_eccs = []
    for ecc in eccs:
        delta_eccs.append(max(ecc))
    return delta_eccs[-1]

def running_average(x_old, n_old, x_new, weight):
    n_new = n_old + weight
    return (x_old * n_old + weight * x_new) / n_new, n_new

def gen_heatmap(xs, ys, ws, x_min, x_max, y_min, y_max, bins):
    arr = np.empty((bins, bins, 2))
    dx = float(x_max - x_min) / bins
    dy = float(y_max - y_min) / bins
    
    for x, y, w in zip(xs, ys, ws):
        x_b = (x - x_min) / dx
	y_b = (y - y_min) / dy
        
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
    count = 0
    for i in range(len(arr)):
        for j in range(len(arr[i])):
	    if arr[i, j, 1] == 0:
	        count += 1
	        arr[i, j, 0] = 0.00
		arr[i, j, 1] = 1
    print "%d bins are 0." % count
    return arr[:,:,0]

#------------------------------------------------------------------------------#
#-The following function will handle all of this project's visualization tools-#
#------------------------------------------------------------------------------#
def generate_visuals(psis, thetas, phis, eccs, path, name, make_movie = False, make_map = True):
    bins = 35
    
    ecc_min = 0.006 
    ecc_max = 0.0125
   
    print ecc_min, ecc_max

    psi_max = 180.0
    psi_min = -180.0
    theta_max = 90.0
    theta_min = -90.0
    phi_max = 180.0
    phi_min = -180.0

    t_bins = 31
    l = np.empty((t_bins, 2))
    dt = (theta_max - theta_min) / t_bins
    for theta, e in zip(thetas, eccs):
        i = int((theta - theta_min) / dt)
	l[i] = max(l[i, 0], e)#running_average(l[i, 0], l[i, 1], e, 1)

    print l

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
	heatmap = gen_heatmap(psis, phis, eccs, psi_min, psi_max, phi_min, phi_max, bins)
	extent = [psi_min, psi_max, phi_min, phi_max]

        plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max)
        plt.ylabel('$\phi [deg]$')
        plt.xlim(psi_min, psi_max)
        plt.ylim(phi_min, phi_max)
        plt.colorbar()
        # x by y plot (lower left corner)
        plt.subplot(223, aspect = 'equal')
        ax = plt.gca()
        ax.set_facecolor('black')
	heatmap = gen_heatmap(psis, thetas, eccs, psi_min, psi_max, theta_min, theta_max, bins)
	extent = [psi_min, psi_max, theta_min, theta_max]
        plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max, aspect=2)
        plt.xlabel('$\psi [deg]$')
        plt.ylabel(r"$\theta [deg]$")
        plt.xlim(psi_min, psi_max)
        plt.ylim(theta_min, theta_max)        
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

        heatmap = gen_heatmap(phis, thetas, eccs, phi_min, phi_max, theta_min, theta_max, bins)
	extent = [phi_min, phi_max, theta_min, theta_max]
	plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max, aspect=2)
        plt.xlabel('$\phi [deg]$')
        plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1475, right=0.95, hspace=0.25, wspace=0.01)
        #plt.subplots_adjust(wspace = 0.001) # shorten the width between left and right since there aren't tick marks
        plt.xlim(phi_min, phi_max)
        plt.ylim(theta_min, theta_max)
	#plt.colorbar()
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
aae = 0
aap = None

#Each data file has its own subdirectory; let's chalk through them in this loop 
for data_dir in data_dirs:

    sample_paths = os.listdir(os.path.join(data_path, data_dir))
    for sample_path in sample_paths:
	fin_data_path = os.path.join(data_path, data_dir, sample_path)
	# If we compute how many stars and planets there are up top, we can ignore doing so in the next loop 
	results = parse_data_from_file(fin_data_path)
	if results == None:
	    continue
	psi, theta, phi, ecc_arr = results
	if theta == 0:
	    print theta
	ecc = get_delta_e_metric(ecc_arr)
	psis.append(psi)
	thetas.append(theta)
	phis.append(phi)
        eccs.append(ecc)
generate_visuals(psis, thetas, phis, eccs, path = None, name = str(fin_data_path).replace('/', '_').replace('__', '_').replace('_.', '.'), make_movie = gen_movie, make_map = True)

print 'There are currently', nan_counter, 'files that have turned into nans.'
