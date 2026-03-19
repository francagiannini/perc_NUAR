
# Bestilling følsomhed på jordtype i NUAR


# Sensitivity of clay content and soil type category in the NUAR calculator  


## Overview
This advisory note evaluates the influence of soil physical properties on nitrogen (N) leaching predictions within the NUAR calculator, which is built upon the empirical NLES5 model. The analysis specifically addresses the sensitivity of nitrogen leaching to changes in clay content (S term) and soil classification (P term) to inform policy decisions regarding farmer objection rights to soil mapping.



## NLES5 Model Equations
Reference leaching (L) is calculated as a function of temporal trends, nitrogen inputs, crop type, and hydrology:

$$L = \tau \,(Y - 1991) + \left(\mu + {\theta_i}N + C\right)^{k} \,(P \cdot S) \,\rho$$

* **Nitrogen (N_func)**: Calculates inputs from topsoil N, mineral N, and organic fertilization.
* **Crop (C_func)**: Models crop-specific leaching based on current and previous crop history.
* **Soil (S_func)**: Calculates the soil term as an exponential function of clay content (CU).
* **Percolation (P_func)**: An indicator function using different parameters for JB 1-3 JB 4-7 soils.

## Script Workflow

The analysis is executed through the following structured steps in the R script perc_nles5_R:

* 0. NLES5 Functions
  - L function NLES5
  - N component
  - C component
  - S component
  - P component

* 1. Clay and Soil type exercise

  - 1.2 Overall Marginal effects NLES5
    - 1.2.1. Load and Prepare Clean Base Data
    - 1.2.2. Helper functions & Fixed Values
    - 1.2.3. Component and Color Definitions
    - 1.2.4. Perform Marginal Effect Analysis
      - 1.2.5.1 Plotting: Faceted Marginal Effects
      - 1.2.5.2 Plotting: Separated Percolation Marginal Effects

  - 1.3 Exploratory percolation distribution

  - 1.4 Marginal effects for percolation variables only, split by soil type (jbnr)
      - 1.4.1. Filter for Percolation Continuous Variables
      - 1.4.2. Point wise check: Is the model actually sensitive to jbnr?

  - 1.5 Interaction P and and JB categories
    - 1.5.1 Gap Analysis
    - 1.5.2 Slopes Analysis - Compare the "jumps" (Sensitivities)
    - 1.5.3 Maximum differences

* 2. Others

## Requirements
This analysis was performed in **R** using the following packages:
- `tidyverse`
- `readxl`
- `ggrepel`
- `patchwork`

## Methodology
The analysis employs a **Marginal Effects** approach to isolate the impact of individual parameters on leaching:

1. **Baseline**: Established using the mean for continuous variables and the mode for categorical variables.

2. **Variation**: Each input variable is tested across its observed range using 100 evenly spaced points.

3. **"Crossed"" Simulation**: Continuous percolation variables (AAa, AAb, APb) are analyzed in combination with JB 1-3 and JB 4-7 (JB > 3) soil types.

4. **Metrics**: The "JB 1-3 Soil Penalty" (absolute difference) and "Marginal Sensitivity" (slope) are calculated to determine soil behavior.

## References
* **Børgesen et al. (2020/2022)**: Parameterization and methodology of the NLES5 model[cite: 4, 146].
* **Eriksen et al. (2024)**: Generalized NUAR equations and reference leaching[cite: 149].