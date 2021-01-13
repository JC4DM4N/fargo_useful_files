#!/bin/bash 

#SBATCH --partition=GPU 
#SBATCH --gres=gpu:2080Ti 
#SBATCH --time=12:00:00 
#SBATCH --ntasks=1 
#SBATCH -o stdout 

#SBATCH --job-name=binary7 

# Where the errors go 
#SBATCH -e log.error 
# Where the messages go 
#SBATCH -o log 

./fargo3d -m -D 0 setups/binary7/binary7.par >& binary7gpu.log 
