#!/bin/sh
####SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --ntasks=12
#SBATCH --nodes=1
####SBATCH --nodelist=atlas05
####SBATCH --exclude=atlas01,atlas03,atlas04,atlas05,atlas06
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --job-name=GCM_clim_extract
#SBATCH --ntasks-per-node=12

#### $1 required, is e.g., batch=1;
#### $2 required, is a domain reference, e.g., domain=\'akcan2km\' or domain=\'world10min\';
#### $3/$4 are regions=TRUE and/or cities=TRUE, defaults to FALSE, at least one must be specified TRUE;
#### $5 optional, is years=1870:2099, or similar;
#### $6 optional, recommended for large number of cities, integer 1-10, e.g., cities.batch=1, run ten times with regions FALSE;
Rscript /workspace/UA/mfleonawicz/projects/SNAPQAQC/code/climate/AR4_AR5_extract.R $1 $2 $3 $4 $5 $6

