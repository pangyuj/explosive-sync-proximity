# Proximity to Explosive Synchronization Determines Critical Transitions

Code for the analyses presented in:

**Proximity to Explosive Synchronization Determines Network Collapse and Recovery Trajectories in Neural and Economic Crises**

U. Lee, H. Kim, M. Kim, G. Oh, P. Joo, A. Park, D. Pal, I. Tracey, C. E. Warnaby, J. Sleigh, and G. A. Mashour

*Proceedings of the National Academy of Sciences (PNAS)* (2025)
https://doi.org/10.1073/pnas.2505434122

## Overview

This repository provides the MATLAB scripts to reproduce the main findings of our study. The code covers three main analyses as described in the paper:
-   **/computational_model**: Contains scripts for the Stuart-Landau model simulations.
-   **/finite_size_scaling_analysis**: Includes a function and an example script to calculate susceptibility from weighted network for finite-size scaling analysis.
-   **/stock_market_analysis**: Includes scripts for analyzing stock market data and generating figures.

## How to Run

1. **Computational Model Analysis:**
- Navigate to `computational_model/` and run `SL_FindingCriticalPoint.m` to identify critical coupling strengths
- Run `SL_StimulationResponses.m` to simulate network responses to perturbations
- The `gong78.mat` file is a synthetic dataset with statistical properties similar to the original 78-node human brain network, which cannot be redistributed.
  
2. **Finite-Size Scaling Analysis:**
- Navigate to `finite_size_scaling_analysis/` and run the scripts to compute cluster-based susceptibility via percolation on weighted networks.

3. **Stock Market Analysis:**
- Navigate to `stock_market_analysis/` and run `Main_Make_Result_StockMarket_2511.m` to reproduce the 2008 crisis analysis

## Contact

For questions about the code, please contact Pangyu Joo at pangyuj@med.umich.edu.
