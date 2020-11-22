#!/bin/bash -l
# submit.h
#
# Description:
#   This script submits a batch job to a SLURM cluster.
#   The job script is defined in job.sh
#
# Parameters:
#   $1: Output filename
#   $2: First job parameter
#   $3: Second job parameter

# Read parameters
filename=$1
p1=$2
p2=$3

# Submit the job script
sbatch -o $filename.out -e $filename.err job.sh $filename $p1 $p2

# Delete empty error files
find -type f -name '*.err' -empty -delete

# Move output files and error files to output directory
mv *.out out/
mv *.err out/
