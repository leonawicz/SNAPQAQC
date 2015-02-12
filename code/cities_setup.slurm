#!/bin/sh
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --job-name=CitSet
#SBATCH --ntasks-per-node=32

#### $1 required, is a domain reference, e.g., domain=\'akcan2km\' or domain=\'world10min\';
#### $2 optional, recommended, integer 1-10, e.g., cities.batch=1, run ten times, will attempt all 10 batches together otherwise;
Rscript /workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/code/cities_setup.R $1 $2

