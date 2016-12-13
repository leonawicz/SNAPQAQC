#!/bin/bash
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
####SBATCH --reservation=shiny_4
#SBATCH -p main
#SBATCH --ntasks=32
####SBATCH --nodelist=atlas03,atlas04,atlas05,atlas06,atlas11,atlas12,atlas13,atlas14,atlas15
####SBATCH --nodelist=atlas05
#SBATCH --job-name=AlfPrep
#SBATCH --nodes=1

# CL variables passed take form, e.g., mainDir=\'/atlas_scratch/mfleonawicz/alfresco\', variable=\'veg\';
# mainDir is optional, default shown above;
# period, required, options are historical or projected;
# variable, required, options are age, veg, and fsv;
# projectName, optional (recommended), defaults to "Unnamed_Project_Run_Extractions"

# reps, optional (recommended), e.g., reps=1:200. You should know the number of total reps in your data and use it explicitly.
# reps defaults to NULL, in which case there is a slight possibility not all reps will be used.
# This could affect the accuracy of burn area, fire frequency, and vegetated area. Not applicable to age distributions.
# It is unlikely to be a problem, but more likely with small spatial domains and/or short time periods.
# Just be safe and specify your reps!

Rscript /atlas_scratch/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfPrep.R $1 $2 $3 $4 $5
