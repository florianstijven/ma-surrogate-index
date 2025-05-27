#!/bin/bash
#SBATCH --nodes=1
#SBATCH --output=par-%J.out
#SBATCH --cpus-per-task=36
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT

ml fhR/4.4.0-foss-2023b

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"

make results/raw-results/simulations/ma-sim-results-proof-of-concept-large.rds