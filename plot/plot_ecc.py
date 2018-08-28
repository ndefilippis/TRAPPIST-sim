'''
Simple script to generate a line plot of eccentricity time series
'''
import matplotlib; matplotlib.use('agg')
import numpy as np, matplotlib.pyplot as plt
from amuse.lab import *
import matplotlib.patches as mpatches

filename1 = "../../distance000.txt"
filename2 = "../../distance001.txt"

def gen_time(n): # Assume simulation for 10,000 years
    dn = (10000.0) / n
    return np.arange(0, 10000.0, dn)

def parse_e(filename): # Parse eccentricities from file
    f = open(filename, "r")
    for _ in range(3): f.readline()

    eccs = f.readline()[2:-3]
    eccs = eccs.split("], [")
    return np.array(list(map(float, eccs[-1].split(", "))))

e1 = parse_e(filename1) 
e2 = parse_e(filename2)
t = gen_time(len(e1))

print len(e1), len(e2), len(t)

fig = plt.figure(figsize=(10, 6))

plt.plot(t, e1, "#000000", linewidth=6) # make an outline on the lighter color
plt.plot(t, e1, "#ffc600", linewidth=4)
plt.plot(t, e2, "#07294d", linewidth=4)

ax = plt.gca()
ax.tick_params(labelsize=20)
ax.set_facecolor('#d0d3d4')
plt.xlim(0, 10000.0)
plt.ylim(0, 0.040)

plt.suptitle("Eccentricity of TRAPPIST-1", fontsize=22)
plt.xlabel("Time [yr]", fontsize=22)
plt.ylabel("Eccentricity", fontsize=22)

perturbed_color = mpatches.Patch(color='#07294d', label='Perturbed System')
unperturbed_color = mpatches.Patch(facecolor='#ffc600', edgecolor='black', label='Unperturbed System')
leg = plt.legend(handles=[perturbed_color, unperturbed_color], loc='upper right', bbox_to_anchor=(1, 1), fontsize=20)#loc = (1.5, 1))#, loc = 'upper right')
leg.get_frame().set_edgecolor('white')

fig.subplots_adjust(right=0.95, bottom=0.15, left=0.15)
plt.savefig("ecc_map.png", dpi=1000)
