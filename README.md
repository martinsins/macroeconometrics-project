# War Episodes and Public Debt: An ARMA-X Analysis

**HEC Lausanne – Macroeconometrics (Jean-Paul Renne)**

## Research Question

How do war episodes dynamically affect the debt-to-GDP ratio across advanced economies? We use an ARMA-X framework with war dummies as exogenous regressors to estimate impulse response functions (IRFs). Switzerland serves as a natural control group in a robustness analysis.

## Project Structure

```
.
├── data/
│   ├── raw/              # Original, unmodified data files
│   └── processed/        # Cleaned and merged datasets
├── R/
│   ├── 00_data_prep.R    # Data download and cleaning
│   ├── 01_descriptive.R  # Descriptive statistics and plots
│   ├── 02_armax.R        # ARMA-X estimation and IRF computation
│   ├── 03_robustness.R   # Switzerland control group analysis
│   └── utils.R           # Shared helper functions
├── output/
│   ├── figures/          # All plots (PDF/PNG)
│   └── tables/           # LaTeX/CSV tables
├── report/               # Final report (PDF + .tex or .Rmd)
├── .gitignore
└── README.md
```

## Methodology

- **Model**: ARMA-X with war episode dummies as exogenous regressors
- **Outcome variable**: Debt-to-GDP ratio
- **Sample**: Advanced economies (TBD), post-WWII
- **IRFs**: Derived from the ARMA-X estimated coefficients
- **Robustness**: Switzerland as a natural control (non-belligerent country)
- **Seasonality**: Addressed if using quarterly/monthly data

## Data Sources

| Dataset | Source | Frequency |
|---|---|---|
| Debt-to-GDP | IMF Historical Public Debt / Jordà-Schularick-Taylor Macrohistory | Annual |
| War episodes | Correlates of War (COW) / UCDP | Annual |
| GDP / Macro controls | World Bank / OECD | Annual |

## Requirements

- R ≥ 4.2
- Key packages: `forecast`, `tseries`, `vars`, `ggplot2`, `dplyr`, `readxl`

Install all dependencies:
```r
install.packages(c("forecast", "tseries", "vars", "ggplot2",
                   "dplyr", "readxl", "tidyr", "patchwork"))
```

## Authors

- [Your name]
- [Colleague's name]

HEC Lausanne, Spring 2026