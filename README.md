# Explosive Synchronization Proximity and Critical Transitions

Code repository for the computational modeling analyses presented in:

**Lee, U.C., Kim, H., Kim, M., Oh, G., Park, A., Joo, P., Pal, D., Tracey, I., Warnaby, C.E., Sleigh, J., & Mashour, G.A.** (2024). Proximity to Explosive Synchronization Determines Network Collapse and Recovery Trajectories in Neural and Economic Crises. *bioRxiv*. https://doi.org/10.1101/2024.11.28.625924

Preprint available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC11642761/

## Description

This repository contains MATLAB code for simulating Stuart-Landau oscillator networks to investigate how proximity to explosive synchronization (ES) affects critical state transitions during perturbations.

## Requirements

- MATLAB R2018b or later
- No additional toolboxes required

## Main Scripts

- `SL_FindingCriticalPoint.m` - Identifies critical points in the network using pair correlation function (PCF)
- `SL_StimulationResponses.m` - Simulates network responses to external perturbations
- `SL_StimulationResponses_RobustnessTest.m` - Tests robustness across different perturbation strengths
- `IE_stuartlandau_distdelay_stim_af.m` - Core Stuart-Landau model with adaptive feedback and distance-dependent delays
- `OrderParameter_Comp.m` - Computes instantaneous order parameters
- `Ort2KPCF.m` - Calculates pair correlation function (PCF)
- `Ort2KACF.m` - Calculates autocorrelation function (ACF)

## Usage

1. Run `SL_FindingCriticalPoint.m` to identify critical coupling strengths for different ES proximities
2. Run `SL_StimulationResponses.m` to simulate perturbation responses at critical points

## Data

Network connectivity matrices (MAT) and distance matrices (Dmat) are not included due to third-party data restrictions. Synthetic networks with similar statistical properties can be generated for testing purposes.

## Citation

If you use this code, please cite:

Lee, U.C., Kim, H., Kim, M., Oh, G., Park, A., Joo, P., Pal, D., Tracey, I., Warnaby, C.E., Sleigh, J., & Mashour, G.A. (2024). Proximity to Explosive Synchronization Determines Network Collapse and Recovery Trajectories in Neural and Economic Crises. *bioRxiv*. https://doi.org/10.1101/2024.11.28.625924

## Contact

For questions, please contact: pangyuj@med.umich.edu
