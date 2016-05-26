#!/bin/bash
##SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
####SBATCH --reservation=shiny_4
#SBATCH -p main
####SBATCH --ntasks=93
#SBATCH --ntasks=108
#SBATCH --job-name=alfExtractMain
#SBATCH --nodes=4
####SBATCH --nodelist=atlas08,atlas09,atlas10
#SBATCH --ntasks-per-node=27

# CL variables passed take form, e.g., modelIndex=1, reps=1:200, years=1900:2100;
# currently, years defaults to 2008:2100, this will change;
# locgroup="Name String" is optional;
# doFire=TRUE, doAgeVeg=TRUE, Rmpi=TRUE and useCRU=FALSE are defaults;

# Rmpi=FALSE means use R's parallel package on a single node;
# useCRU=TRUE means outputs are parsed identically on any given Alfresco "model", but are relabeled properly as historical CRU-based outputs afterward.
# This is intended for running on known historical years only, e.g., 1900-2007 for IEM Alfresco outputs.
# For known GCM-based years, leave useCRU=FALSE to avoid improper relabeling.

# mpiBy="year" is default. Alternative is "rep". Called 'mpiBy' but also applies when uses the parallel package on a single node;
# repSample, optional, e.g., repSample=30, means randomly sample 30 replicates from the full set
# projectName, optional (recommended), defaults to "Unnamed_Project_Run_Extractions"
# readMethod, optional, defaults to "loop" for serial raster loading. Alternative is "stack" for read all files pre-loop.
# You will run out of memory with "stack" if you have too many years and reps.
mpirun --oversubscribe -np 1 Rscript /atlas_scratch/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfExtractMain.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}
