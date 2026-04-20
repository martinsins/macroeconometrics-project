# ==============================================================================
# 00_data_prep.R
# War Episodes and Public Debt — Data Preparation
# HEC Lausanne, Macroeconometrics, Spring 2026
# ==============================================================================
# This script:
#   1. Loads the raw JST Macrohistory database
#   2. Filters to the 5 countries of interest, 1946–2015
#   3. Selects and renames relevant variables
#   4. Runs basic sanity checks
#   5. Saves a clean panel to data/processed/panel_clean.rds
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 0. Paths
# ------------------------------------------------------------------------------

path_raw   <- "data/raw/JSTdatasetR6.xlsx"  
path_clean <- "data/processed/panel_clean.rds"

# ------------------------------------------------------------------------------
# 1. Load raw JST data
# ------------------------------------------------------------------------------

jst_raw <- read_excel(path_raw, sheet = 1)

cat("Dimensions of raw JST data:", dim(jst_raw), "\n")
cat("Years available:", range(jst_raw$year, na.rm = TRUE), "\n")

# ------------------------------------------------------------------------------
# 2. Filter countries and sample period
# ------------------------------------------------------------------------------

COUNTRIES <- c("USA", "GBR", "FRA", "AUS", "CHE")
YEAR_MIN  <- 1950
YEAR_MAX  <- 2020

panel <- jst_raw %>%
  filter(iso %in% COUNTRIES, year >= YEAR_MIN, year <= YEAR_MAX)

cat("\nRows after filter:", nrow(panel),
    "— expected:", length(COUNTRIES) * (YEAR_MAX - YEAR_MIN + 1), "\n")

# ------------------------------------------------------------------------------
# 3. Select and rename variables
# ------------------------------------------------------------------------------

panel <- panel %>%
  select(
    country, iso, year,
    debtgdp, gdp, cpi, iy, ca, stir, ltrate
  ) %>%
  rename(
    debt_gdp = debtgdp,
    inv_rate = iy,
    curr_acc = ca,
    rate_st  = stir,
    rate_lt  = ltrate
  )

panel <- panel %>%
  mutate(debt_gdp = debt_gdp * 100)
# ------------------------------------------------------------------------------
# 4. Derive additional variables
# ------------------------------------------------------------------------------

panel <- panel %>%
  arrange(iso, year) %>%
  group_by(iso) %>%
  mutate(
    inflation  = (cpi / lag(cpi) - 1) * 100,
    gdp_growth = (gdp / lag(gdp) - 1) * 100,
    d_debt_gdp = debt_gdp - lag(debt_gdp)
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 5. War dummies
# ------------------------------------------------------------------------------
# Common wars (at least 3 of 4 belligerents):
#   Korean War      1950-1953  USA, GBR, FRA (minor), AUS
#   Gulf War        1990-1991  USA, GBR, FRA, AUS
#   Afghanistan     2001-2014  USA, GBR, FRA, AUS
#
# Country-specific:
#   Indochina       FRA 1950-1954
#   Vietnam         USA + AUS 1965-1973
#   Algeria         FRA 1954-1962
#   Falklands       GBR 1982
#   Iraq War        USA + GBR + AUS 2003-2011
#
# CHE is never a belligerent — all dummies remain 0 by construction.
# ------------------------------------------------------------------------------

panel <- panel %>%
  mutate(
    
    war_common = case_when(
      iso %in% c("USA","GBR","FRA","AUS") & year %in% 1950:1953 ~ 1,
      iso %in% c("USA","GBR","FRA","AUS") & year %in% 1990:1991 ~ 1,
      iso %in% c("USA","GBR","FRA","AUS") & year %in% 2001:2014 ~ 1,
      TRUE ~ 0
    ),
    
    war_indochina = if_else(iso == "FRA" & year %in% 1950:1954, 1, 0),
    war_vietnam   = if_else(iso %in% c("USA","AUS") & year %in% 1965:1973, 1, 0),
    war_algeria   = if_else(iso == "FRA" & year %in% 1954:1962, 1, 0),
    war_falklands = if_else(iso == "GBR" & year == 1982, 1, 0),
    war_iraq      = if_else(iso %in% c("USA","GBR","AUS") & year %in% 2003:2011, 1, 0),
    
    war_any = if_else(
      war_common == 1 | war_indochina == 1 | war_vietnam == 1 |
        war_algeria == 1 | war_falklands == 1 | war_iraq == 1,
      1, 0
    )
  )

cat("\n--- War involvement (years at war) by country ---\n")
panel %>%
  group_by(iso) %>%
  summarise(
    common    = sum(war_common),
    vietnam   = sum(war_vietnam),
    algeria   = sum(war_algeria),
    falklands = sum(war_falklands),
    iraq      = sum(war_iraq),
    any_war   = sum(war_any)
  ) %>%
  print()

# ------------------------------------------------------------------------------
# 6. Sanity checks
# ------------------------------------------------------------------------------

cat("\n--- Missing values by variable ---\n")
panel %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0) %>%
  print()

cat("\n--- Debt/GDP summary by country ---\n")
panel %>%
  group_by(iso) %>%
  summarise(
    n      = n(),
    mean   = round(mean(debt_gdp, na.rm = TRUE), 1),
    sd     = round(sd(debt_gdp,   na.rm = TRUE), 1),
    min    = round(min(debt_gdp,  na.rm = TRUE), 1),
    max    = round(max(debt_gdp,  na.rm = TRUE), 1),
    n_miss = sum(is.na(debt_gdp))
  ) %>%
  print()

# ------------------------------------------------------------------------------
# 7. Visualisation — debt/GDP with war periods shaded
# ------------------------------------------------------------------------------

war_shade <- panel %>%
  filter(war_any == 1) %>%
  group_by(iso) %>%
  mutate(grp = cumsum(c(1, diff(year) > 1))) %>%
  group_by(iso, grp) %>%
  summarise(xmin = min(year) - 0.5, xmax = max(year) + 0.5, .groups = "drop")

p_debt <- ggplot(panel, aes(x = year, y = debt_gdp)) +
  geom_rect(
    data = war_shade,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  facet_wrap(~ iso, scales = "free_y") +
  labs(
    title    = "Government Debt / GDP (1950-2020)",
    subtitle = "Shaded areas = years at war (war_any = 1)",
    x        = NULL,
    y        = "Debt / GDP (%)",
    caption  = "Source: Jorda-Schularick-Taylor Macrohistory Database"
  ) +
  theme_minimal(base_size = 11)

ggsave("output/figures/debt_gdp_wars.png", p_debt,
       width = 10, height = 6, dpi = 150)
cat("\nPlot saved to output/figures/debt_gdp_wars.png\n")

# ------------------------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------------------------

saveRDS(panel, path_clean)
cat("\nClean panel saved to", path_clean, "\n")
cat("Final dimensions:", dim(panel), "\n")

panel <- readRDS("data/processed/panel_clean.rds")
glimpse(panel)

panel %>%
  group_by(iso) %>%
  summarise(across(c(debt_gdp, inflation, gdp_growth, d_debt_gdp),
                   ~ sum(is.na(.)), .names = "na_{.col}"))      

panel %>%
  filter(iso == "FRA", war_any == 1) %>%
  select(year, war_common, war_indochina, war_algeria) %>%
  print(n = 40)
