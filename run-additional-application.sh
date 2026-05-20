#!/bin/bash
#SBATCH --nodes=1
#SBATCH --output=par-%J.out
#SBATCH --cpus-per-task=36
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT

module load cluster/wice/batch
module load R/4.4.2-gfbf-2024a
module load R-bundle-CRAN/2024.11-foss-2024a

export OMP_NUM_THREADS=1

# Rscript -e "install.packages(c('BBmisc', 'delayed', 'caret')"
# Rscript -e "install.packages('tlverse-sl3-v1.4.4-243-g0e8f236.tar.gz', type = 'source', repos = NULL)"
Rscript -e "renv::restore()"

make additional-application