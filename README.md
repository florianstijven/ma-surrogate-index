# Meta-Analytic Evaluation of Surrogacy Using (Estimated) Surrogate Indices

This directory contains the code used for the data application and simualtions in the manuscript "Evaluation of Surrogate Endpoints Based on Meta-Analysis with Surrogate Indices" by Stijven, F. and Gilbert, P. B. (2025). This repo is an R-project and contains all files (except data) needed to reproduce the results from the manuscript. Files that are produced by running code (i.e., tables and figures) are not included in this repo. 

## Project Structure

The project is organized into the following directories:
* R/: Contains the R scripts and functions used for the analyses.
* data/: Directory where the data should be stored. This directory also 
  contains Rscripts to preprocess the original data and to generate a synthetic data set from the original data.
  Note that the original data are not included in this repo because they cannot be shared.
  The synthetic data are shared included, however.
* results/: When running the code, the results (including tables and figures) will be saved here.

These directories contain their own README files. 
The Makefile allows one to run all Rscripts in the correct order by running `make all` in the console. 

## Reproducibility

TO MODIFY

The results with the synthetic data can be reproduced by executing the Makefile
(by running `make all` in the command line). This will run all statistical
analyses (which produces results saved to `results/raw-results/`) and will
produce all plots and tables presented in the paper (which are saved to
`results/figures/` and `results/tables/`). Note that this requires the presence
of `data/processed_data_synthetic.csv` (which is the synthetic data set that was 
simulated to resemble the original data). This file is not included in this repo 
because it is too large, but is shared as a supplementary file to the manuscript.

Reproducibility is also facilitated through the use of the `renv` R package (see also below). 

TO MODIFY

Since the analyses are computationally intensive, they were run on a
Slurm-managed High Performance Computing (HPC) cluster, for which we used 
the `runanalysis.sh` batch script. If
these analyses are run on a personal computer, there should be enough RAM
available because the sandwich estimator based on the `geex` R package is greedy
in its use of RAM. Running the Makefile is expected to take 10 hours on a recent
personal computer.

## renv

The project uses the `renv` R package for dependency management. This helps to
ensure that the code can be run across devices using similar environments (i.e.,
using the same versions of the required R packages). More 
information about the use of `renv` can be found [here](https://rstudio.github.io/renv/articles/renv.html).
