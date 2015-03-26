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
#### $2 optional, is a CRU version reference, e.g., cru=32, defaults to 32 if not specified;
#### $3/$4 are regions=TRUE and/or cities=TRUE, defaults to FALSE, at least one must be specified TRUE;
#### $5, optional, is years=1870:2099, or similar;
#### $6 optional, not generally needed for CRU (see GCM slurm script), cities.batch is an integer 1 through 10;
Rscript /workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/code/CRU_extract.R $1 $2 $3 $4 $5 $6
