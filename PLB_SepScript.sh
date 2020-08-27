#!/bin/bash

#SBATCH --exclude=compute-1-[1-33]

cwd=$(pwd)
for i in $(seq 0 $4)
do
    sbatch --export=ALL,workdir=$cwd,dir=$1,bf=$2,bg=$3,spatialrundir=SpatialOut,rundir=SepOut,runnum=$i --job-name=Sep_$i -o FIXEDSep_$i.out plb_sep.job
done





