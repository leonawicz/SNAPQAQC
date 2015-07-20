#!/bin/bash
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --ntasks=200
####SBATCH --nodelist=atlas03,atlas04,atlas05,atlas06,atlas11,atlas12,atlas13,atlas14,atlas15
####SBATCH --nodelist=atlas05
#SBATCH --job-name=RmpiAlfStats
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=29

# CL variables passed take form, e.g., model.index=1, reps=1:200, years=1900:2100;
# locgroup=1:5 is optional;
# doFire=TRUE, doAgeVeg=TRUE, and Rmpi=TRUE are defaults;
# Rmpi=FALSE means use R's parallel package on a single node;
# repSample, optional, e.g., repSample=30, means randomly sample 30 replicates from the full set
mpirun -np 1 Rscript /workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/code/alfStatsByRep_Rmpi.R $1 $2 $3 $4 $5 $6 $7 $8
