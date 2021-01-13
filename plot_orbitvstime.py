# plot orbital semi-major axis vs time from fargo orbit*.dat file output.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-file')
args = parser.parse_args()
file = args.file

orbit = np.genfromtxt(file)

plt.plot(orbit[:,0],orbit[:,2])
plt.ylabel('Semi-major axis, AU')
plt.xlabel('Time')
plt.show()
