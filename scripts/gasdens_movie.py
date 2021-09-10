# plot multi-frame of desired gasdens*.dat output file from fargo as a
# FuncAnimation movie.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import os
import argparse

# number of files to include in the movie. Assumes indices start at 0.
parser = argparse.ArgumentParser()
parser.add_argument('-nfiles')
args = parser.parse_args()
nfiles = int(args.nfiles)

ifiles = np.arange(nfiles+1)

#initialise fig
fig = plt.figure()
rho = np.fromfile('gasdens0.dat').reshape(128,384)
im = plt.imshow(np.log10(rho),origin='lower',cmap=cm.Oranges_r,aspect='auto',
                animated=True)
#function for funcanimate animation
def animate(i):
    rho = np.fromfile('gasdens%i.dat' %(i*10)).reshape(128,384)
    print('Plotting dump %i...' %(i*10))
    im = plt.imshow(np.log10(rho),origin='lower',cmap=cm.Oranges_r,aspect='auto',
                    animated=True)
    return im

# call the animator
anim = animation.FuncAnimation(fig, animate, frames=nfiles, interval=200)
plt.show()
