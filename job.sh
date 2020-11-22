#!/bin/bash -l
#SBATCH --exclusive
#SBATCH -p aolin.q

hostname
echo

# Configure environment
module add gcc/9.2.0

# Read parameters
filename=$1
epochs=$2
seed=$3

# Compile
gcc -Ofast TSP.c -o $filename.exe

perf stat -d $filename.exe $epochs $seed 2>&1
perf record -o $filename.data $filename.exe $epochs $seed 2>&1
