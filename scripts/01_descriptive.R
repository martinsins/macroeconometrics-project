# ==============================================================================
# 01_descriptive.R
# War Episodes and Public Debt — Descriptive Statistics & Pre-estimation Checks
# HEC Lausanne, Macroeconometrics, Spring 2026
# ==============================================================================
# This script:
#   1. Summary statistics table (for report Data section)
#   2. Publication-ready debt/GDP plot with war shading
#   3. First-difference of debt/GDP plot
#   4. ADF unit root tests on debt_gdp and d_debt_gdp
#   5. ACF/PACF plots of d_debt_gdp by country
#   6. Raw correlation between war_any and d_debt_gdp
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(tseries)   # adf.test
library(forecast)  # Acf, Pacf

# ------------------------------------------------------------------------------
# 0. Load cleaned panel
# ------------------------------------------------------------------------------

panel <- readRDS("data/processed/panel_clean.rds")

if (!dir.exists("output/figures")) dir.create("output/figures", recursive = TRUE)
if (!dir.exists("output/tables"))  dir.create("output/tables",  recursive = TRUE)

BELLIGERENTS <- c("USA", "GBR", "FRA", "AUS")

# ------------------------------------------------------------------------------
# 1. Summary statistics table
# ------------------------------------------------------------------------------

sumstats <- panel %>%
  filter(iso %in% BELLIGERENTS) %>%
  group_by(iso) %>%
  summarise(
    across(
      c(debt_gdp, d_debt_gdp, inflation, gdp_growth, rate_lt),
      list(
        mean = ~ round(mean(.x, na.rm = TRUE), 2),
        sd   = ~ round(sd(.x,   na.rm = TRUE), 2),
        min  = ~ round(min(.x,  na.rm = TRUE), 2),
        max  = ~ round(max(.x,  na.rm = TRUE), 2)
      ),
      .names = "{.col}__{.fn}"
    )
  ) %>%
  pivot_longer(-iso, names_to = c("variable", "stat"), names_sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value)

cat("\n--- Summary Statistics (belligerents only) ---\n")
print(sumstats, n = 50)

write.csv(sumstats, "output/tables/summary_stats.csv", row.names = FALSE)
cat("Summary stats saved to output/tables/summary_stats.csv\n")

# ------------------------------------------------------------------------------
# 2. Publication-ready debt/GDP plot with war shading
# ------------------------------------------------------------------------------

# Country labels for facets
country_labels <- c(
  AUS = "Australia", GBR = "United Kingdom",
  FRA = "France",    USA = "United States", CHE = "Switzerland"
)

war_shade <- panel %>%
  filter(war_any == 1) %>%
  group_by(iso) %>%
  mutate(grp = cumsum(c(1, diff(year) > 1))) %>%
  group_by(iso, grp) %>%
  summarise(xmin = min(year) - 0.5, xmax = max(year) + 0.5, .groups = "drop")

p1 <- ggplot(panel, aes(x = year, y = debt_gdp)) +
  geom_rect(
    data = war_shade,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  facet_wrap(~ iso, scales = "free_y",
             labeller = labeller(iso = country_labels)) +
  labs(
    title   = "Government Debt-to-GDP Ratio (1950-2020)",
    subtitle = "Shaded areas indicate years of war involvement",
    x       = NULL,
    y       = "Debt / GDP (%)",
    caption = "Source: Jorda-Schularick-Taylor Macrohistory Database"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave("output/figures/fig1_debt_gdp.png", p1,
       width = 10, height = 6, dpi = 150)
cat("Figure 1 saved to output/figures/fig1_debt_gdp.png\n")

# ------------------------------------------------------------------------------
# 3. First difference of debt/GDP
# ------------------------------------------------------------------------------

p2 <- ggplot(panel %>% filter(!is.na(d_debt_gdp)),
             aes(x = year, y = d_debt_gdp)) +
  geom_rect(
    data = war_shade,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "red", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.7, colour = "steelblue") +
  facet_wrap(~ iso, scales = "free_y",
             labeller = labeller(iso = country_labels)) +
  labs(
    title    = "First Difference of Debt-to-GDP Ratio (1950-2020)",
    subtitle = "Shaded areas indicate years of war involvement",
    x        = NULL,
    y        = "Change in Debt / GDP (pp)",
    caption  = "Source: Jorda-Schularick-Taylor Macrohistory Database"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave("output/figures/fig2_d_debt_gdp.png", p2,
       width = 10, height = 6, dpi = 150)
cat("Figure 2 saved to output/figures/fig2_d_debt_gdp.png\n")

# ------------------------------------------------------------------------------
# 4. ADF unit root tests
# ------------------------------------------------------------------------------
# Null hypothesis of ADF: series has a unit root (non-stationary)
# p < 0.05 -> reject unit root -> series is stationary
# We test both debt_gdp (levels) and d_debt_gdp (first difference)
# ------------------------------------------------------------------------------

cat("\n--- ADF Unit Root Tests ---\n")
cat(sprintf("%-6s %-15s %-10s %-10s %s\n",
            "ISO", "Variable", "ADF stat", "p-value", "Conclusion"))
cat(strrep("-", 60), "\n")

adf_results <- list()

for (ctry in BELLIGERENTS) {
  
  for (varname in c("debt_gdp", "d_debt_gdp")) {
    
    x <- panel %>%
      filter(iso == ctry, !is.na(.data[[varname]])) %>%
      pull(varname)
    
    test   <- adf.test(x)
    pval   <- round(test$p.value, 3)
    stat   <- round(test$statistic, 3)
    concl  <- if_else(pval < 0.05, "Stationary", "Non-stationary")
    
    cat(sprintf("%-6s %-15s %-10s %-10s %s\n",
                ctry, varname, stat, pval, concl))
    
    adf_results[[paste(ctry, varname, sep = "_")]] <- data.frame(
      iso = ctry, variable = varname,
      adf_stat = stat, p_value = pval, conclusion = concl
    )
  }
}

adf_table <- bind_rows(adf_results)
write.csv(adf_table, "output/tables/adf_tests.csv", row.names = FALSE)
cat("\nADF results saved to output/tables/adf_tests.csv\n")

# ------------------------------------------------------------------------------
# 5. ACF / PACF plots of d_debt_gdp by country
# ------------------------------------------------------------------------------
# These guide the choice of AR(p) and MA(q) orders in the ARMA-X.
# ACF  -> suggests MA order q (cuts off after lag q)
# PACF -> suggests AR order p (cuts off after lag p)
# ------------------------------------------------------------------------------

png("output/figures/fig3_acf_pacf.png",
    width = 1400, height = 900, res = 120)
par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))

for (ctry in BELLIGERENTS) {
  
  x <- panel %>%
    filter(iso == ctry, !is.na(d_debt_gdp)) %>%
    pull(d_debt_gdp)
  
  Acf(x,  main = paste(ctry, "— ACF of d_debt_gdp"),  lag.max = 15)
  Pacf(x, main = paste(ctry, "— PACF of d_debt_gdp"), lag.max = 15)
}

dev.off()
cat("Figure 3 (ACF/PACF) saved to output/figures/fig3_acf_pacf.png\n")

# ------------------------------------------------------------------------------
# 6. Raw correlation: war_any vs d_debt_gdp
# ------------------------------------------------------------------------------
# Simple preview of the relationship before modelling.
# Not causal — just checking the sign and magnitude are in the right ballpark.
# ------------------------------------------------------------------------------

cat("\n--- Raw correlation: war_any vs d_debt_gdp ---\n")
panel %>%
  filter(iso %in% BELLIGERENTS, !is.na(d_debt_gdp)) %>%
  group_by(iso) %>%
  summarise(
    mean_peace = round(mean(d_debt_gdp[war_any == 0], na.rm = TRUE), 2),
    mean_war   = round(mean(d_debt_gdp[war_any == 1], na.rm = TRUE), 2),
    diff       = round(mean(d_debt_gdp[war_any == 1], na.rm = TRUE) -
                         mean(d_debt_gdp[war_any == 0], na.rm = TRUE), 2),
    correlation = round(cor(war_any, d_debt_gdp, use = "complete.obs"), 3)
  ) %>%
  print()

cat("\n--- Descriptive script complete ---\n")
cat("Key outputs:\n")
cat("  output/tables/summary_stats.csv\n")
cat("  output/tables/adf_tests.csv\n")
cat("  output/figures/fig1_debt_gdp.png\n")
cat("  output/figures/fig2_d_debt_gdp.png\n")
cat("  output/figures/fig3_acf_pacf.png\n")
cat("\nNext step: inspect ADF results and ACF/PACF plots,\n")
cat("then use them to choose ARMA orders in 02_armax.R\n")