#!/bin/sh
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --ntasks=1
#SBATCH --nodes=1
####SBATCH --nodelist=atlas09
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --job-name=CRU_clim_extract

#### $1 is a domain reference, either akcan2km or world10min; $2, optional, is years=1870:2099, or similar.
Rscript /workspace/UA/mfleonawicz/Leonawicz/Projects/2014/AR4_AR5_comparisons/code/CRU_extract.R $1 $2
