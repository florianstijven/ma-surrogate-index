# Meta-Analytic Evaluation of Surrogacy Using (Estimated) Surrogate Indices

This repo contains the code used for the data application and simulations in the
manuscript "Evaluation of Surrogate Endpoints Based on Meta-Analysis with
Surrogate Indices" by Stijven, F. and Gilbert, P. B. (2025). This repo is an
R project and contains all files (except data) needed to reproduce the results
from the manuscript. Files that are produced by running code (i.e., tables and
figures) are not included in this repo.

## Project Structure

The project is organized into the following directories:
* `R/`: Contains the R scripts used for the analyses.
* `data/`: Directory where the data should be stored. This directory also 
  contains R scripts to preprocess the original data and to generate a synthetic 
  data set from the original data. Note that the original data are not included
  in this repo because they cannot be shared. The synthetic data are included, 
  however.
* `results/`: When running the code, the results (including unprocessed results,
  tables and figures) will be saved here.

These directories contain their own README files which explain the role of each
R script in the directory. The R scripts themselves also contain many comments
explaining the code in more detail. The README files are useful, however, to
better understand the overall structure of the code and how the different R
scripts interact with each other.

## Reproducibility

To improve reproducibility of the results reported on in the manuscript, we used
a Makefile for coordinating the order in which scripts should be run and with
which arguments, and we used renv for R package management.

### Makefile

To enhance reproducibility, we used a Makefile. This file contains all
information that is needed to run the R scripts in the correct order and with
the correct command-line arguments to reproduce the results reported on in the
manuscript. Specifically, all simulation and data-application results will be
produced by running `make all` in the command line. Note that the data
application code can only be run if the original data are present in the correct
location (i.e., `data/processed_data.csv` should exist).

The Makefile contains three phony targets which facilitate running subsets of
all code used in the manuscript:

1. `make simulation` will only run the simulations.
2. `make application` will only run the data application that uses the original
data. This recipe will not run correctly unless the original data are present in
the correct location.
4. `make application-synthetic` will run the data application that uses the
synthetic data. Because the synthetic data are included in this repo, this
recipe should run correctly after cloning the repo.

### renv

The project uses the `renv` R package for dependency management. This helps to
ensure that the code can be run across devices using similar environments (i.e.,
using the same versions of the required R packages). More information about the
use of `renv` can be found
[here](https://rstudio.github.io/renv/articles/renv.html).
