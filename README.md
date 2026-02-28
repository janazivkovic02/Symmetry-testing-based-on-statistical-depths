# Symmetry Testing Based on Statistical Depths

This repository contains the implementation and simulation study accompanying the seminar project on testing central symmetry using depth-based rank statistics.

## Project Overview

The project investigates a rank-based test for central symmetry in multivariate distributions. The performance of the test is evaluated via Monte Carlo simulations under different distributional settings and skewness levels.

Five depth notions are considered:

- Simplicial depth
- L_2 depth
- Mahalanobis depth
- Projection depth
- Zonoid depth

Empirical rejection rates are computed under:

1. Azzalini-type skewing
2. Azzalini skewing with contamination

### Description of Files

- `main.R`  
  Main simulation script. Runs Monte Carlo experiments and produces figures.

- `R/data_generation.R`  
  Generates samples under various symmetric kernels with Azzalini skewing.

- `R/data_generation2.R`  
  Generates contaminated samples.

- `R/rn_test_statistic.R`  
  Implements depth-based rank statistics and associated test procedures.

## Theoretical Background and Reference

This project is based on the methodology developed in:

R. Dyckerhoff, C. Ley, and D. Paindaveine (2014).  
"Depth-based Runs Tests for Bivariate Central Symmetry".  

The definition of the depth-based runs statistic, its asymptotic null distribution 
(with variance 11/3), and the simulation design follow the above reference.

The implementation of the simulations and numerical experiments was developed 
independently for the purposes of this seminar project.
