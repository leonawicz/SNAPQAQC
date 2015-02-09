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

#### batch argument not required for CRU as it is for GCM processing
#### $1 required, is a domain reference, e.g., domain=\'akcan2km\' or domain=\'world10min\';
#### $2/$3 are regions=TRUE and/or cities=TRUE, defaults to FALSE, at least one must be specified TRUE;
#### $4, optional, is years=1870:2099, or similar;
#### $5 optional, not generally needed for CRU (see GCM slurm script), cities.batch is an integer 1 through 10;
Rscript /workspace/UA/mfleonawicz/leonawicz/projects/AR4_AR5_comparisons/code/CRU_extract.R $1 $2 $3 $4 $5
