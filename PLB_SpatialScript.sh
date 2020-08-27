#!/bin/bash

cwd=$(pwd)
for i in $(seq 0 $4)
do
    sbatch --export=ALL,workdir=$cwd,dir=$1,bf=$2,bg=$3,rundir=SpatialOut,runnum=$i --job-name=Spatial_$i -o Spatial_$i.out plb_spatial.job
done


