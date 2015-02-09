#!/bin/bash
#SBATCH --exclusive
#SBATCH --mail-type=END
#SBATCH --mail-user=mfleonawicz@alaska.edu
#SBATCH --account=snap
#SBATCH -p main
#SBATCH --ntasks=200
#SBATCH --nodelist=atlas01,atlas03,atlas04,atlas05,atlas06,atlas11,atlas12,atlas13,atlas15
#SBATCH --job-name=RmpiAlfStats
#SBATCH --nodes=9
#SBATCH --ntasks-per-node=23

# CL variables passed take form, e.g., model.index=1, reps=1:200, years=1900:2100
mpirun -np 1 Rscript /workspace/UA/mfleonawicz/leonawicz/projects/SNAPQAQC/alfStatsByRep_Rmpi.R $1 $2 $3

# If mpirun call above is commented out and the lines below are not, pass any dummy model.index CL arg
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=1 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=2 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=3 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=4 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=5 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=6 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=7 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=8 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=9 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=10 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=11 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=12 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=13 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=14 $2 $3
#mpirun -np 1 Rscript /big_scratch/mfleonawicz/Rmpi/alfStatsByRep_Rmpi.R model.index=15 $2 $3

