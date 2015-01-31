#!/bin/sh
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --job-name=GCM_clim_extract
#SBATCH --ntasks-per-node=12

#### $1 is batch=1; $2 is a domain reference, either akcan2km or world10min; $3, optional, is years=1870:2099, or similar.
Rscript /workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/code/AR4_AR5_extract.R $1 $2 $3

