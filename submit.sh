#!/bin/bash -l
#SBATCH --exclusive
#SBATCH -p aolin.q

hostname
echo

# Exports
module add gcc/9.2.0
gcc -Ofast TSP.c -o TSP

perf stat -d TSP 2>&1
perf record TSP 2>&1
