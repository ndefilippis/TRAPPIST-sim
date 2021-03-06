#------------------------------------------------------------------------------#
#--------------------------------plot_angles.py--------------------------------#
#------------------------------------------------------------------------------#
#--------------------------Created by Nick DeFilippis--------------------------#
#--------------------------Inspired by Mark Giovinazzi-------------------------#
#------------------------------------------------------------------------------#
#-This program was created for the purpose of creating heatmaps for encounters-#
#-of differently oriented perturbers. The value of the heatmap depends on how--#
#------much of an affect the perturber had on the eccentricity of the outer----#
#-----------------------------------planet.------------------------------------#
#------------------------------------------------------------------------------#
#---------------------------Date Created: 07/15/2018---------------------------#
#------------------------Date Last Modified: 08/28/2018------------------------#
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
    '''
    Parse data file. Files are in the format
    psi, theta, phi
    <Time series of eccentricity>
    '''
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
    '''
    Values used for the heatmap.
    Using the maximum eccentricity of the outermost planet
    '''
    delta_eccs = []
    for ecc in eccs:
        delta_eccs.append(max(ecc))
    return delta_eccs[-1]

def running_average(x_old, n_old, x_new, weight):
    '''
    Calculate running average based on old average, old count, new data, and its respective weight
    '''
    n_new = n_old + weight
    return (x_old * n_old + weight * x_new) / n_new, n_new

def gen_heatmap(xs, ys, ws, x_min, x_max, y_min, y_max, bins):
    '''
    Generate the heatmap that we are going to plot
    '''
    arr = np.empty((bins, bins, 2)) # Heatmap array
    dx = float(x_max - x_min) / bins
    dy = float(y_max - y_min) / bins
    
    for x, y, w in zip(xs, ys, ws):
        '''
        For each data point, distribute it among the four closest bins, with a 
	weight respective to how close it is to each bin.
	'''
        x_b = (x - x_min) / dx
	y_b = (y - y_min) / dy
        
	j = int(x_b) # bin index
	i = int(y_b) # bin index
	xf = x_b - j # how much of the point is in the bin
	yf = y_b - i

        # Deal with edge case
	if j+1 >= bins:
	    xf = 1
	if i+1 >= bins:
	    yf = 1
	
	# Keep track of running average for each bin

	# lower-left
	x_old, n_old = arr[i, j]
	arr[i, j] = running_average(x_old, n_old, w, xf*yf)
        
        # upper-right
	if j+1 < bins and i+1 < bins:
	    x_old, n_old = arr[i+1, j+1]
	    arr[i+1,j+1] = running_average(x_old, n_old, w, (1-xf) * (1-yf))
        # lower right
	if j+1 < bins:
	    x_old, n_old = arr[i, j+1]
	    arr[i, j+1] = running_average(x_old, n_old, w, (1-xf) * yf)
	# upper left
	if i+1 < bins:
	    x_old, n_old = arr[i+1, j]
	    arr[i+1, j] = running_average(x_old, n_old, w, xf * (1-yf))
	
    # deal with zero values
    count = 0
    for i in range(len(arr)):
        for j in range(len(arr[i])):
	    if arr[i, j, 1] == 0:
	        count += 1
	        arr[i, j, 0] = get_average(arr, i, j) # get average of surrounding values
		arr[i, j, 1] = 0
    print "%d bins are 0." % count
    return arr[:,:,0]

def get_average(arr, i, j):
    '''
    Get average values of surrounding squares, weighted by how close they are to this square
    '''
    c = 0
    t = 0
    for ii in range(-3, 4): # take a 7x7 grid around this square, 3 squares in each direction
        for jj in range(-3, 4):
            i_n = i + ii
	    j_n = j + jj
	    if i_n >= 0 and i_n < len(arr) and j_n >= 0 and j_n < len(arr):
	        x, n = arr[i_n, j_n]
		if n != 0: # if this square isnt already empty
                    c += 1 / float(ii*ii + jj*jj) # inverse square weighting based on distance
		    t += x / float(ii*ii + jj*jj)
    if c != 0:
        return t / c
    else:
        return 0


#------------------------------------------------------------------------------#
#-The following function will handle all of this project's visualization tools-#
#------------------------------------------------------------------------------#
def generate_visuals(psis, thetas, phis, eccs, path, name, make_movie = False, make_map = True):
    '''
    Generate three separate heatmaps, one for each pair of Euler angles
    psis - list of psis
    thetas - list of thetas
    phis - list of phis
    eccs - list of eccentricities
    path - where to save image
    name - name of image file
    '''
    bins = 35 # number of bins in each dimension
    
    ecc_min = 0.006 # upper and lower limit for heatmap values
    ecc_max = 0.0125
   
    print ecc_min, ecc_max

    # boundary for Euler angles
    psi_max = 180.0
    psi_min = -180.0
    theta_max = 90.0
    theta_min = -90.0
    phi_max = 180.0
    phi_min = -180.0

    # Print out histogram for theta distribution
    # Not necessary, just something to look at
    t_bins = 31
    l = np.empty((t_bins, 2))
    dt = (theta_max - theta_min) / t_bins
    for theta, e in zip(thetas, eccs):
        i = int((theta - theta_min) / dt)
	l[i] = max(l[i, 0], e)#running_average(l[i, 0], l[i, 1], e, 1)

    print l

    # If the user wants to make trail maps in their run, the following code will be executed
    if make_map == True:
        fig = plt.figure(figsize = (11.5, 10))
        #print xs[0], zs[0], cs[0], ss[0]
        # psi v phi plot (top left corner)
        plt.subplot(221, aspect = 'equal')
        ax = plt.gca()
        ax.tick_params(labelsize=20)
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
        plt.yticks(np.arange(-180.0, 181.0, 90.0))
	plt.ylabel('$\phi$ [deg]', fontsize=22)
        plt.xlim(psi_min, psi_max)
        plt.ylim(phi_min, phi_max)

        # psi v theta plot (lower left corner)
        plt.subplot(223, aspect = 'equal')
        ax = plt.gca()
	ax.tick_params(labelsize=20)
        ax.set_facecolor('black')
	heatmap = gen_heatmap(psis, thetas, eccs, psi_min, psi_max, theta_min, theta_max, bins)
	extent = [psi_min, psi_max, theta_min, theta_max]
        plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max, aspect=2)
        plt.yticks(np.arange(-90.0, 91.0, 90.0))
	plt.xticks(np.arange(-180.0, 181.0, 90.0))
	plt.xlabel('$\psi$ [deg]', fontsize=22)
        plt.ylabel(r"$\theta$ [deg]",fontsize=22)
        plt.xlim(psi_min, psi_max)
        plt.ylim(theta_min, theta_max)        
        #ax.set_xticklabels([min, max])
        #ax.set_yticklabels([min, max])
        #ax.set_xicks([min, max])
        #ax.set_yicks([min, max])
        # phi by theta plot (lower right corner)
        plt.subplot(224, aspect = 'equal')
        ax = plt.gca()
	ax.tick_params(labelsize=20)
        ax.set_facecolor('black')
        ax.set_yticklabels([])
        #[plt.scatter(float(zs[i]), float(ys[i]), c = cs[i], s = int(ss[i])) for i in range(len(xs))]

        heatmap = gen_heatmap(phis, thetas, eccs, phi_min, phi_max, theta_min, theta_max, bins)
	extent = [phi_min, phi_max, theta_min, theta_max]
	plt.xticks(np.arange(-180.0, 181.0, 90.0))
	im = plt.imshow(heatmap, extent=extent, cmap='viridis', vmin=ecc_min, vmax=ecc_max, aspect=2)
        plt.xlabel('$\phi$ [deg]',fontsize=22)
        plt.subplots_adjust(top=0.9, bottom=0.1, left=0.08, right=0.92, hspace=0.25, wspace=0.01)
        plt.subplots_adjust(wspace = 0.001) # shorten the width between left and right since there aren't tick marks
        plt.xlim(phi_min, phi_max)
        plt.ylim(theta_min, theta_max)

	plt.subplots_adjust(right=0.85)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.025, 0.7])
	cbar_ax.tick_params(labelsize=20)
	cbar = plt.colorbar(im, cax=cbar_ax)
	cbar.set_label('Eccentricity', fontsize=22)
	#plt.colorbar()
        #leg = plt.colorbar(loc='upper right', bbox_to_anchor=(0.9, 2))#loc = (1.5, 1))#, loc = 'upper right')
        #leg.get_frame().set_edgecolor('white')

        # make some plot_radius function to determine value for s for stars and planets...
        # stars and planets should look noticeably different in size, but stars/stars and planets/planets should be subtle
        plt.suptitle("Eccentricity versus Orientation", fontsize=35)
        print 'saving///' # just for the user's sake..
	plt.savefig('angle_maps/' + str(name) + '.png', dpi = 1000)
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
	# Parse each data file
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
