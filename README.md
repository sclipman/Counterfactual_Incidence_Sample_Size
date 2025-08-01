# Sample Size Calculation for Single-Arm HIV Trials with Counterfactual Comparator

## Overview

This R script calculates required sample sizes, statistical power, and variance components for single-arm HIV prevention trials using a counterfactual comparator derived from HIV-recency assays. The statistical methods closely follow Gao et al. (2021).

## Usage

- Ensure R is installed along with required packages (`stats`).
- Set trial parameters directly within the script (see examples at the end).
- Run the script in R or RStudio.

```R
source("sample_size_counterfactual_HIV_trial.R")
```

## Input Parameters

| Parameter    | Meaning                                                    |
|--------------|------------------------------------------------------------|
| `p`          | HIV prevalence at baseline                                 |
| `OmegaT`     | Mean duration of recent infection (MDRI), in years         |
| `betaT`      | False recent rate (FRR)                                    |
| `sigmaOmegaT`| Standard error of MDRI                                     |
| `sigmabetaT` | Standard error of FRR                                      |
| `R1`         | Alternative hypothesis incidence rate ratio                |
| `R0`         | Null hypothesis incidence rate ratio (default = 1)         |
| `tau`        | Follow-up time, in years                                   |
| `Time`       | Recency cut-off, in years                                  |
| `r`          | Proportion of HIV-negative individuals enrolled            |
| `alpha`      | Significance level (two-sided: 0.05)                       |
| `power`      | Desired statistical power (default = 0.8 or 0.9)           |

