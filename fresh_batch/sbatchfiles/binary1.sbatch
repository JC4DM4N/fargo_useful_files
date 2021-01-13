#!/bin/bash 

#SBATCH --partition=GPU 
#SBATCH --gres=gpu:2080Ti 
#SBATCH --time=12:00:00 
#SBATCH --ntasks=1 
#SBATCH -o stdout 

#SBATCH --job-name=binary1 

# Where the errors go 
#SBATCH -e log.error 
# Where the messages go 
#SBATCH -o log 

./fargo3d -m -D 0 setups/binary1/binary1.par >& binary1gpu.log 
