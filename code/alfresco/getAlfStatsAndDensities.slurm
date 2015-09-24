#!/bin/bash
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
#SBATCH -p main
####SBATCH --ntasks=32
####SBATCH --nodelist=atlas03,atlas04,atlas05,atlas06,atlas11,atlas12,atlas13,atlas14,atlas15
####SBATCH --nodelist=atlas05
#SBATCH --job-name=AlfStatsDens
#SBATCH --nodes=1
####SBATCH --ntasks-per-node=32

# CL variables passed take form, e.g., mainDir=\'/big_scratch/mfleonawicz/Rmpi/outputs\', variable=\'veg\'
# mainDir is optional, default shown above;
# variable is required, options are age, veg, abfc, and fsv;
# multiple variables can be processed serially in one job if passed with proper escaping as, e.g., variable=c\(\'age\',\'veg\'\);
# age input sample size currently known and hardcoded in script as n.samples.in=100
# density output sample size currently hardcoded in script as n.samples=1000
mpirun -np 1 Rscript /workspace/UA/mfleonawicz/projects/SNAPQAQC/code/alfresco/getAlfStatsAndDensities.R $1 $2
