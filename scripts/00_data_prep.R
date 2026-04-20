# =============================================================================
# 00_data_prep.R
# War Episodes & Debt/GDP – Data download and cleaning
# HEC Lausanne – Macroeconometrics project
# =============================================================================

library(dplyr)
library(readxl)
library(tidyr)

# --- 1. Load Debt/GDP data ---------------------------------------------------
# Suggested source: IMF Historical Public Debt Database
# or Jordà-Schularick-Taylor Macrohistory Database (https://www.macrohistory.net/database/)

# debt_raw <- read_csv("data/raw/imf_debt.csv")

# --- 2. Load War episode data ------------------------------------------------
# Suggested source: Correlates of War (COW) - https://correlatesofwar.org/
# or UCDP - https://ucdp.uu.se/

# war_raw <- read_csv("data/raw/cow_war.csv")

# --- 3. Create war dummies ---------------------------------------------------
# For each country-year, create a binary dummy = 1 if country was at war
# and a "post-war" dummy for the N years following conflict end

# --- 4. Merge datasets -------------------------------------------------------
# panel <- debt_clean %>%
#   left_join(war_dummies, by = c("country", "year")) %>%
#   left_join(macro_controls, by = c("country", "year"))

# --- 5. Save processed data --------------------------------------------------
# saveRDS(panel, "data/processed/panel_clean.rds")

cat("Data preparation script – fill in data loading steps above.\n")
