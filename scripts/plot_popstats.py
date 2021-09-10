# plot population statistics of output files from fargo3d.

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nfol')
args = parser.parse_args()
nfol = int(args.nfol)

#first read data from binary runs
i=0
times_binary = []
radii_binary = []
times_single = []
radii_single = []
while i<nfol:
    try:
        bindat = np.genfromtxt('binarygpu%i/orbit0.dat' %i)
        times_binary.append([bindat[0,0],bindat[-1,0]])
        radii_binary.append([bindat[0,2],bindat[-1,2]])
        sindat = np.genfromtxt('singlegpu%i/orbit0.dat' %i)
        times_single.append([sindat[0,0],sindat[-1,0]])
        radii_single.append([sindat[0,2],sindat[-1,2]])
    except:
        # if the folder doesn't exist
        pass
    i+=1

times_binary = np.asarray(times_binary)
times_single = np.asarray(times_single)
radii_binary = np.asarray(radii_binary)
radii_single = np.asarray(radii_single)

dadt_single = (radii_single[:,0]-radii_single[:,1])/(times_single[:,1])
dadt_binary = (radii_binary[:,0]-radii_binary[:,1])/(times_binary[:,1])

print('Mean binary finish time: %.3f' %np.mean(times_binary[:,1]))
print('Mean solo finish time: %.3f' %np.mean(times_single[:,1]))
print('')
print('Mean da/dt single star system: %.3f' %np.mean(dadt_single))
print('Mean da/dt binary star system: %.3f' %np.mean(dadt_binary))

def bin_semimajor_axes(values):
    # manually bin planet semi-major axes between 0 and 100au
    radii = np.arange(1,100)
    binned_values = [sum((values < r) & (values >= r-1)) for r in radii]
    return radii, binned_values

plt.figure('rplot_binary')
radii, binned_a = bin_semimajor_axes(radii_binary[:,0])
plt.bar(radii,binned_a,label='initial',alpha=0.5,rwidth=0.9)
radii, binned_a = bin_semimajor_axes(radii_binary[:,1])
plt.bar(radii,binned_a,label='initial',alpha=0.5,rwidth=0.9)

#plt.hist(radii_binary[:,0],label='initial',alpha=0.5,rwidth=0.9)
#plt.hist(radii_binary[:,1],label='final',alpha=0.5,rwidth=0.9)
plt.legend()
plt.grid()
plt.xlabel('R (AU)')

plt.figure('rplot_single')
radii, binned_a = bin_semimajor_axes(radii_single[:,0])
plt.bar(radii,binned_a,label='initial',alpha=0.5,rwidth=0.9)
radii, binned_a = bin_semimajor_axes(radii_single[:,1])
plt.bar(radii,binned_a,label='initial',alpha=0.5,rwidth=0.9)

#plt.hist(radii_single[:,0],label='initial',alpha=0.5,rwidth=0.9)
#plt.hist(radii_single[:,1],label='final',alpha=0.5,rwidth=0.9)
plt.legend()
plt.grid()
plt.xlabel('R (AU)')

plt.figure('dadt')
plt.hist(dadt_single,label='single system',alpha=0.5,rwidth=0.9)
plt.hist(dadt_binary,label='binary system',alpha=0.5,rwidth=0.9)
plt.legend()
plt.grid()
plt.show()
