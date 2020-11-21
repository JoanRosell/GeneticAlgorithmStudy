#!/bin/bash -l
#SBATCH --exclusive
#SBATCH -p aolin.q
# Parameters:
#   $1 Epochs
#   $2 seed

hostname
echo

# Exports
module add gcc/9.2.0
gcc -Ofast TSP.c -o TSP

perf stat -d TSP $1 $2 2>&1
perf record TSP $1 $2 2>&1
