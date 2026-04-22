---
title: "Military Deployment Intensity and Public Debt: An ARMA-X Analysis for the United States"
output: github_document
---

**HEC Lausanne – Macroeconometrics (Jean-Paul Renne)**

## Research Question

Does the intensity of US military engagement — measured by the number of troops deployed abroad — have a significant and persistent effect on the federal debt-to-GDP ratio? We estimate an ARMA-X model on a single quarterly US time series (post-WWII), treating troop deployment as a continuous exogenous regressor, and derive the implied impulse response functions (IRFs).

## Motivation

Prior work on war and public finances often relied on binary war dummies in panel settings. This conflates very different levels of military engagement and introduces cross-country heterogeneity that is hard to control for. By focusing on a single country and using a continuous deployment measure, we stay within a genuine time series framework and exploit the natural variation in engagement intensity across US military conflicts since 1950 (Korea, Vietnam, Gulf War, Iraq, Afghanistan).

## Project Structure

```
.
├── data/
│   ├── raw/                     # Original, unmodified data files
│   │   ├── fred_debt.csv        # Federal debt-to-GDP from FRED
│   │   └── troops_deployed.csv  # Quarterly US troops deployed abroad (DoD / SIPRI)
│   └── processed/
│       └── us_ts_data.csv       # Merged, cleaned quarterly time series
├── R/
│   ├── 00_data_prep.R           # Data download, merging, stationarity treatment
│   ├── 01_descriptive.R         # Time series plots, ACF/PACF, unit root tests
│   ├── 02_armax.R               # ARMA-X estimation (MLE), model selection (AIC/BIC)
│   ├── 03_irf.R                 # Impulse response functions from ARMA-X coefficients
│   ├── 04_robustness.R          # Lag sensitivity, sub-period analysis, alt. X vars
│   └── utils.R                  # Shared helper functions
├── output/
│   ├── figures/                 # All plots (PDF/PNG)
│   └── tables/                  # LaTeX/CSV tables
├── report/                      # Final report (PDF + .Rmd or .tex source)
├── README.Rmd                   # Edit this file, then knit to regenerate README.md
├── .gitignore
└── README.md
```

## Methodology

- **Model**: ARMA(p,q)-X — univariate time series with exogenous regressor
- **Dependent variable (y_t)**: First difference of US federal debt-to-GDP ratio (quarterly)
- **Exogenous variable (x_t)**: Number of US troops deployed abroad (in thousands, quarterly)
- **Estimation**: Maximum Likelihood Estimation (MLE)
- **Model selection**: AIC / BIC over a grid of (p, q) orders
- **Stationarity**: ADF and KPSS tests on both series; differencing applied as needed
- **IRFs**: Derived analytically from ARMA-X estimated coefficients to trace the dynamic debt response to a deployment shock
- **Robustness**:
  - Lagged deployment (x_{t-1}, x_{t-2}) to capture delayed fiscal effects
  - Sub-period analysis (Cold War vs. post-9/11)
  - Alternative intensity measure: defense expenditure growth (FRED: `FDEFX`)

## Model Specification

$$\Delta \text{Debt}_t = \mu + \sum_{i=1}^{p} \phi_i \Delta \text{Debt}_{t-i} + \sum_{j=1}^{q} \theta_j \varepsilon_{t-j} + \beta \cdot \text{Troops}_t + \varepsilon_t, \quad \varepsilon_t \sim \mathcal{N}(0, \sigma^2)$$

$\beta$ captures the contemporaneous effect of deployment intensity on the change in the debt ratio. The ARMA component controls for the intrinsic persistence in debt dynamics.

## Data Sources

| Variable | Dataset | Source | Frequency |
|---|---|---|---|
| Federal debt-to-GDP | `GFDEGDQ188S` | FRED (Federal Reserve Bank of St. Louis) | Quarterly |
| Troops deployed abroad | DMDC / SIPRI historical records | DoD Defense Manpower Data Center | Quarterly (constructed) |
| Defense expenditure (robustness) | `FDEFX` | FRED (Federal Reserve Bank of St. Louis) | Quarterly |

> **Note on troop deployment data**: The quarterly series is manually constructed from DoD/DMDC historical figures, covering Korea (1950–53), Vietnam (1964–73), Gulf War (1990–91), and Iraq/Afghanistan (2001–2011). Deployment is set to a non-zero baseline during peacetime to reflect the permanent non-combat overseas US military presence.

## Requirements

- R ≥ 4.2
- Key packages: `forecast`, `tseries`, `ggplot2`, `dplyr`, `readxl`, `tidyr`, `patchwork`

```r
install.packages(c("forecast", "tseries", "ggplot2",
                   "dplyr", "readxl", "tidyr", "patchwork"))
```

## Authors

- Martin Sinishta
- Natalia Wyszatycki

HEC Lausanne, Spring 2026
