#!/bin/bash


cwd=$(pwd)
for i in $(seq 0 $4)
do
    sbatch --export=ALL,workdir=$cwd,dir=$1,bf=$2,bg=$3,spatialrundir=SpatialOut,rundir=SepOut,runnum=$i --job-name=MS_$i -o MS_$i.out plb_makeseps.job
done





