#!/bin/bash
##SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --ntasks=93
#SBATCH --job-name=alfExtractMain
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=16

# CL variables passed take form, e.g., model.index=1, reps=1:200, years=1900:2100;
# currently, years defaults to 2008:2100, this will change;
# locgroup="Name String" is optional;
# doFire=TRUE, doAgeVeg=TRUE, and Rmpi=TRUE are defaults;
# Rmpi=FALSE means use R's parallel package on a single node;
# mpiBy="rep" is default. Alternative is "year". Called 'mpiBy' but also applies when uses the parallel package on a single node;
# repSample, optional, e.g., repSample=30, means randomly sample 30 replicates from the full set
mpirun -np 1 Rscript /workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/B1RERUN_alfExtractMain.R $1 $2 $3 $4 $5 $6 $7 $8 $9
