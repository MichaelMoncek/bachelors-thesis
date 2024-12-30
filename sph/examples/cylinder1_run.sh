#!/bin/bash

#SBATCH --job-name=cylinder1
#SBATCH --output=cylinder1.txt
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --time=1:00:00
#SBATCH -p express3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=michael.moncek@seznam.cz

# module add Arch/linux-ubuntu18.04-ivybridge
module add julia
julia -t 36 cylinder1_run.jl
