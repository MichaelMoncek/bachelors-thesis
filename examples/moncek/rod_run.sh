#!/bin/bash

#SBATCH --job-name=rod
#SBATCH --output=rod.txt
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=10:00
#SBATCH -p express3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michael.moncek@seznam.cz

# module add Arch/linux-ubuntu18.04-ivybridge
module add julia
julia -t 2 rod_run.jl
