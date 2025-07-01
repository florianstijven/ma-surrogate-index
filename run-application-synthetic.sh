#!/bin/bash
#SBATCH --output=par-%J.out
#SBATCH --ntasks=1 --cpus-per-task=36 --nodes=1
#SBATCH --time=3:00:00
#SBATCH --cluster=wice
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH -A lp_doctoralresearch

ml fhR/4.4.0-foss-2023b
export OMP_NUM_THREADS=1

# Rscript -e "install.packages(c('BBmisc', 'delayed', 'caret')"
# Rscript -e "install.packages('tlverse-sl3-v1.4.4-243-g0e8f236.tar.gz', type = 'source', repos = NULL)"
Rscript -e "renv::restore()"

make application-synthetic