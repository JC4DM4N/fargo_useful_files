# perform latin hypercube sampling across the desired parameter space of disc properties.
# Specify desired params and their corresponding upper and lower bounds.
# Resultant file setups will be written into a 'fresh_batch/' subdirectory, along with a
# corresponding csv containing the setups params.

import numpy as np
import matplotlib.pyplot as plt
import lhsmdu
import os
import shutil

#set random seed
lhsmdu.setRandomSeed(111)

#param limits
Rp_lims = np.asarray([20.0,35.0])                   # initial planet semi-major axis, AU
iRp = 0
Mpl_lims = np.asarray([1.0,15.0])*0.0009543         # inital planet mass, Solar masses
iMpl = 1
a_lims = np.asarray([100.0,250.0])                  # initial binary separation, AU
ia = 2
Ms1_lims = np.asarray([0.85,1.2])                   # Primary star mass, Solar masses
iMs1 = 3
Ms2_lims = np.asarray([0.08,0.4])                   # Secondary star mass, Solar masses
iMs2 = 4
alpha_lims = np.asarray([1e-4,1e-2])                # alpha viscosity limits
ialpha = 5
pindex_lims = np.asarray([0.0,1.5])                 # Surface density profile exponent
ipindex = 6
H_lims = np.asarray([0.025,0.1])                    # Disc aspect ratio
iH = 7
Md_lims = np.asarray([0.01,0.1])                    # Disc mass
iMd = 8

headers = 'Rp,Mpl,a,Ms1,Ms2,alpha,pindex,H,Mdisc' #column headers for csv later

lims = [Rp_lims,Mpl_lims,a_lims,Ms1_lims,Ms2_lims,alpha_lims,pindex_lims,H_lims,Md_lims]
nparams = len(lims)

nsamp = 10

lhs = np.asarray(lhsmdu.sample(nparams,nsamp))   # perform latin hypercube sampling

# scale random variables to be within desired limits
for i,lim in enumerate(lims):
    lhs[i,:] = lhs[i,:]*(lim[1]-lim[0]) + lim[0]
    lhs[i,:] = lhs[i,:]

# calculate required grid size and timestep vals
grid_rad = lhs[iRp,:] + 10    # 10 au greater than the planet location
torb = 2*np.pi*grid_rad**(3./2.) # orbital period at grid outer radius
dt = torb/100.
tvisc = 4./9.*(grid_rad**2/lhs[ialpha,:]**2/lhs[iH,:]**2)*np.sqrt(grid_rad**3/lhs[iMs1,:])

# calc sigma0 from disc mass and pindex
# NOTE: this assumes this disc mass extent is 100au NOT just the size of the grid
sig0 = lhs[iMd,:]/np.pi/(100**(3-lhs[ipindex,:])/(3-lhs[ipindex,:]))

"""
Now write these params to the setup file fargo.par and planet file
"""
# first check the fresh_batch directory and subdirectories exist...
if os.path.exists('fresh_batch'):
    shutil.rmtree('fresh_batch')
# then setup desired directories
os.mkdir('fresh_batch')
os.mkdir('fresh_batch/setups')
os.mkdir('fresh_batch/planets')
for i in range(nsamp):
    os.mkdir('fresh_batch/planets/binary%i' %i)
    os.mkdir('fresh_batch/planets/single%i' %i)
    os.mkdir('fresh_batch/setups/binary%i' %i)
    os.mkdir('fresh_batch/setups/single%i' %i)

#now write setup files
for i in range(nsamp):
    # Start with planet.cfg file
    temp = open('setuptemplates/template.cfg','r')
    binfile = open('fresh_batch/planets/binary%i/binary%i.cfg' %(i,i), 'w')
    solofile = open('fresh_batch/planets/single%i/single%i.cfg' %(i,i), 'w')
    for j,line in enumerate(temp.readlines()):
        if j<10:
            # first 10 lines are header
            binfile.write(line)
            solofile.write(line)
        elif j==10:
            # first write line contating planet properties
            binfile.write('Planet   %.2f    %.3f    0.0     YES     YES \n' %(lhs[iRp,i],lhs[iMpl,i]))
            solofile.write('Planet   %.2f    %.3f    0.0     YES     YES \n' %(lhs[iRp,i],lhs[iMpl,i]))
        elif j==11:
            # now write line contating secondary star properties
            binfile.write('Companion   %.2f    %.3f    0.0     NO     NO \n' %(lhs[ia,i],lhs[iMs2,i]))
    temp.close()
    binfile.close()
    solofile.close()

    #now do .par files in setups directory
    binfile = open('fresh_batch/setups/binary%i/binary%i.par' %(i,i), 'w')
    solofile = open('fresh_batch/setups/single%i/single%i.par' %(i,i), 'w')
    for j,file in enumerate([binfile, solofile]):
        file.write('Setup			fargo \n')
        file.write('\n')
        file.write('### Disk parameters \n')
        file.write('\n')
        file.write('AspectRatio         %.3f            Thickness over Radius in the disc \n' %lhs[iH,i])
        file.write('Sigma0		        %.10f       	Surface Density at r=1 \n' %sig0[i])
        file.write('Alpha			    %.4f \n' %lhs[ialpha,i])
        file.write('SigmaSlope		    %.2f		    Slope for the surface density \n' %lhs[ipindex,i])
        file.write('FlaringIndex		0.0		        Slope for the aspect-ratio \n')
        file.write('\n')
        file.write('# Radial range for damping (in period-ratios). Values smaller than one \n')
        file.write('# prevent damping. \n')
        file.write('\n')
        file.write('DampingZone 1.15 \n')
        file.write('# Characteristic time for damping, in units of the inverse local \n')
        file.write('# orbital frequency. Higher values means lower damping \n')
        file.write('\n')
        file.write('TauDamp 0.3 \n')
        file.write('\n')
        file.write('### Planet parameters \n')
        file.write('\n')
        if j==0:
            file.write('PlanetConfig		planets/binary%i/binary%i.cfg \n' %(i,i))
        elif j==1:
            file.write('PlanetConfig		planets/single%i/single%i.cfg \n' %(i,i))
        file.write('ThicknessSmoothing      0.6 \n')
        file.write('RocheSmoothing          0.0 \n')
        file.write('Eccentricity		0.0 \n')
        file.write('ExcludeHill		no \n')
        file.write('IndirectTerm		Yes \n')
        file.write('\n')
        file.write('### Mesh parameters \n')
        file.write('\n')
        file.write('Nx			384		Azimuthal number of zones \n')
        file.write('Ny                      128		Radial number of zones \n')
        file.write('Xmin			-3.14159265358979323844 \n')
        file.write('Xmax			3.14159265358979323844 \n')
        file.write('Ymin			1.0   Inner boundary radius \n')
        file.write('Ymax			%.1f  Outer boundary radius \n' %grid_rad[i])
        file.write('OmegaFrame              1.0005		Angular velocity for the frame of reference (If Frame is F). \n')
        file.write('Frame			G		Method for moving the frame of reference \n')
        file.write('\n')
        file.write('### Output control parameters \n')
        file.write('\n')
        file.write('DT			%.5f		Physical time between fine-grain outputs \n' %dt[i])
        file.write('Ninterm	    100		    Number of DTs between scalar fields outputs \n')
        file.write('Ntot		%i		Total number of DTs \n' %int(tvisc[i]/dt[i]))
        file.write('\n')
        if j==0:
            file.write('OutputDir		@outputs/binarygpu%i \n' %i)
        elif j==1:
            file.write('OutputDir		@outputs/singlegpu%i \n' %i)
        file.write('\n')
        file.write('### Plotting parameters \n')
        file.write('\n')
        file.write('PlotLog			yes \n')
    binfile.close()
    solofile.close()

    # copy other required files from the setuptemplates/ into setups/ directory
    src_files = os.listdir('setuptemplates/templatedir/')
    for file_name in src_files:
        full_file_name = os.path.join('setuptemplates/templatedir/', file_name)
        if os.path.isfile(full_file_name):
            shutil.copy(full_file_name, 'fresh_batch/setups/single%i/' %i)
            shutil.copy(full_file_name, 'fresh_batch/setups/binary%i/' %i)

# also need to write slurm files
os.mkdir('fresh_batch/sbatchfiles/')
for i in range(nsamp):
    for str in ['single','binary']:
        file = open('fresh_batch/sbatchfiles/%s%i.sbatch' %(str,i), 'w')
        file.write('#!/bin/bash \n')
        file.write('\n')
        file.write('#SBATCH --partition=GPU \n')
        file.write('#SBATCH --gres=gpu:2080Ti \n')
        file.write('#SBATCH --time=12:00:00 \n')
        file.write('#SBATCH --ntasks=1 \n')
        file.write('#SBATCH -o stdout \n')
        file.write('\n')
        file.write('#SBATCH --job-name=%s%i \n' %(str,i))
        file.write('\n')
        file.write('# Where the errors go \n')
        file.write('#SBATCH -e log.error \n')
        file.write('# Where the messages go \n')
        file.write('#SBATCH -o log \n')
        file.write('\n')
        file.write('./fargo3d -m -D 0 setups/%s%i/%s%i.par >& %s%igpu.log \n' %(str,i,str,i,str,i))
        file.close()

# finally, write LHS parameters to a csv for reference
with open('fresh_batch/setups.csv', 'w') as f:
    np.savetxt(f, lhs.T, delimiter=',', header=headers)
