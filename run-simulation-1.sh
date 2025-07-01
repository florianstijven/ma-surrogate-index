#!/bin/bash
#SBATCH --output=par-%J.out
#SBATCH --ntasks=1 --cpus-per-task=36 --nodes=1
#SBATCH --time=60:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH -A lp_doctoralresearch

ml fhR/4.4.0-foss-2023b

export OMP_NUM_THREADS=1

Rscript -e "renv::restore()"

make results/raw-results/simulations/ma-sim-results-proof-of-concept-small.rds
make results/raw-results/simulations/ma-sim-results-vaccine-small.rds
