#!/bin/bash
##SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --ntasks=93
#SBATCH --job-name=alfExtractMain
#SBATCH --nodes=4
#SBATCH --nodelist=atlas01,atlas13,atlas14,atlas15
#SBATCH --ntasks-per-node=24

# CL variables passed take form, e.g., model.index=1, reps=1:200, years=1900:2100;
# currently, years defaults to 2008:2100, this will change;
# locgroup="Name String" is optional;
# doFire=TRUE, doAgeVeg=TRUE, Rmpi=TRUE and useCRU=FALSE are defaults;

# Rmpi=FALSE means use R's parallel package on a single node;
# useCRU=TRUE means outputs are parsed identically on any given Alfresco "model", but are relabeled properly as historical CRU-based outputs afterward.
# This is intended for running on known historical years only, e.g., 1900-2007 for IEM Alfresco outputs.
# For known GCM-based years, leave useCRU=FALSE to avoid improper relabeling.

# mpiBy="rep" is default. Alternative is "year". Called 'mpiBy' but also applies when uses the parallel package on a single node;
# repSample, optional, e.g., repSample=30, means randomly sample 30 replicates from the full set
mpirun -np 1 Rscript /workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/B1RERUN_alfExtractMain.R $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
