#!/bin/bash
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
#SBATCH --reservation=shiny_4
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
Rscript /atlas_scratch/mfleonawicz/projects/SNAPQAQC/code/alfresco/alfPrep.R $1 $2 $3 $4
