# plot single frame of desired gasdens*.dat output file from fargo

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse

# name of the gasdens file to be plotted.
parser = argparse.ArgumentParser()
parser.add_argument('-file')
args = parser.parse_args()
file = args.file

rho = np.fromfile(file).reshape(128,384)

#plot figure
fig = plt.figure()
im = plt.imshow(np.log10(rho),origin='lower',cmap=cm.Oranges_r,aspect='auto')
plt.show()
