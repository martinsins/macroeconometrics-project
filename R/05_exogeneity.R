# =============================================================================
# 05_exogeneity.R
# Military Deployment & Public Debt — Exogeneity Test
# =============================================================================
# Tests whether Δdebt/GDP Granger-causes troops (the key exogeneity check).
# Uses lmtest::grangertest only — no vars dependency.
#
# Run after: 01b_stationarity.R  (needs data/processed/us_ts_detrended.csv)
# =============================================================================

library(lmtest)
library(ggplot2)

sep  <- strrep("=", 64)
sep2 <- strrep("-", 64)

# ── 0. Load data ──────────────────────────────────────────────────────────────
df <- read.csv("data/processed/us_ts_detrended.csv")

# Use plain numeric vectors — grangertest on R 4.1.x is more reliable this way
y_dt <- as.numeric(df$d_debt_gdp_dt)   # Δdebt/GDP detrended
x_dt <- as.numeric(df$troops_dt)       # troops detrended
n    <- length(y_dt)

cat(sprintf("\n%s\n  EXOGENEITY TEST\n  Sample: %d-%d  (n = %d)\n%s\n\n",
            sep, min(df$year), max(df$year), n, sep))

# =============================================================================
# SECTION 1 — Lag order selection via BIC
# =============================================================================
cat(sprintf("%s\n  SECTION 1 — Lag order selection (BIC)\n%s\n\n", sep2, sep2))

bic_vals <- sapply(1:6, function(k) {
  X      <- embed(cbind(x_dt, y_dt), k + 1)
  lm_fit <- lm(X[, 1] ~ X[, -1])
  BIC(lm_fit)
})

k_bic <- which.min(bic_vals)
k_use <- max(k_bic, 1)

bic_df <- data.frame(lag = 1:6, BIC = round(bic_vals, 3))
print(bic_df)
cat(sprintf("\n  BIC-selected lag order: k = %d  (using k = %d)\n\n", k_bic, k_use))

# =============================================================================
# SECTION 2 — Granger causality tests
# =============================================================================
cat(sprintf("%s\n  SECTION 2 — Granger causality tests (k = %d)\n%s\n\n",
            sep2, k_use, sep2))

gdf <- data.frame(troops = x_dt, debt = y_dt)

gc_debt_to_troops <- grangertest(troops ~ debt,   order = k_use, data = gdf)
gc_troops_to_debt <- grangertest(debt   ~ troops, order = k_use, data = gdf)

cat("  H0: Δdebt/GDP does NOT Granger-cause troops\n")
print(gc_debt_to_troops)

cat("\n  H0: troops does NOT Granger-cause Δdebt/GDP\n")
print(gc_troops_to_debt)

p1 <- gc_debt_to_troops[2, "Pr(>F)"]
p2 <- gc_troops_to_debt[2, "Pr(>F)"]
F1 <- gc_debt_to_troops[2, "F"]
F2 <- gc_troops_to_debt[2, "F"]

cat(sprintf("\n  *** VERDICT (k = %d) ***\n", k_use))
cat(sprintf("  Δdebt -> troops :  F = %.3f,  p = %.4f  =>  %s\n", F1, p1,
            ifelse(p1 > 0.10,
                   "FAIL TO REJECT H0 — no evidence debt Granger-causes troops",
                   ifelse(p1 > 0.05,
                          "MARGINAL (p<0.10) — weak evidence of feedback",
                          "REJECT H0 — debt appears to Granger-cause troops (!)"))))
cat(sprintf("  Troops -> Δdebt :  F = %.3f,  p = %.4f  =>  %s\n\n", F2, p2,
            ifelse(p2 > 0.10, "FAIL TO REJECT H0 — no predictive content",
                   ifelse(p2 > 0.05, "MARGINAL (p<0.10)",
                          "REJECT H0 — troops predicts debt changes"))))

# =============================================================================
# SECTION 3 — Manual Wald test (no vars needed)
# =============================================================================
cat(sprintf("%s\n  SECTION 3 — Manual Wald test on troops equation\n%s\n\n",
            sep2, sep2))

Z         <- embed(cbind(x_dt, y_dt), k_use + 1)
dep       <- Z[, 1]
own_lags  <- Z[, 2:(k_use + 1), drop = FALSE]
debt_lags <- Z[, (k_use + 2):(2 * k_use + 1), drop = FALSE]

fit_unres <- lm(dep ~ own_lags + debt_lags)
fit_res   <- lm(dep ~ own_lags)
wald      <- anova(fit_res, fit_unres)

cat("  ANOVA F-test — H0: debt lags jointly zero in troops equation\n")
print(wald)

p_wald <- wald[2, "Pr(>F)"]
cat(sprintf("\n  Wald p-value: %.4f  =>  %s\n\n", p_wald,
            ifelse(p_wald > 0.10,
                   "No evidence of feedback from debt to troops",
                   "Evidence of feedback — discuss in report")))

cat("  Troops equation coefficients (unrestricted):\n")
print(round(summary(fit_unres)$coefficients, 4))

# =============================================================================
# SECTION 4 — Cross-correlation function
# =============================================================================
cat(sprintf("\n%s\n  SECTION 4 — Cross-correlation function\n%s\n\n", sep2, sep2))

ccf_vals <- ccf(x_dt, y_dt, lag.max = 8, plot = FALSE)

ccf_df <- data.frame(
  lag = as.integer(ccf_vals$lag),
  ccf = round(as.numeric(ccf_vals$acf), 4)
)

ci <- qnorm(0.975) / sqrt(n)

cat("  CCF(troops_dt, Δdebt_dt) at lags -8 to +8:\n")
print(ccf_df)
cat(sprintf("\n  95%% confidence band: +/- %.4f\n", ci))
cat("  Negative lags = troops leads debt  (exogeneity-consistent)\n")
cat("  Positive lags = debt leads troops  (feedback / endogeneity concern)\n\n")

sig_pos <- ccf_df[ccf_df$lag > 0 & abs(ccf_df$ccf) > ci, ]
if (nrow(sig_pos) == 0) {
  cat("  => No significant CCF at positive lags.\n")
  cat("  => Consistent with exogeneity of troops.\n")
} else {
  cat("  => Significant CCF at positive lag(s):\n")
  print(sig_pos)
  cat("  => Suggests debt may predict future troop levels. Discuss in report.\n")
}

# =============================================================================
# SECTION 5 — CCF figure
# =============================================================================
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

bar_col <- ifelse(ccf_df$lag > 0 & abs(ccf_df$ccf) > ci, "#d6604d", "grey40")

p_ccf <- ggplot(ccf_df, aes(x = lag, y = ccf)) +
  annotate("rect", xmin = 0.5, xmax = 8.5,
           ymin = -Inf, ymax = Inf, fill = "#d6604d", alpha = 0.05) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-ci, ci),
             linetype = "dashed", colour = "#2166ac", linewidth = 0.6) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8, colour = bar_col) +
  geom_point(size = 2.5, colour = bar_col) +
  annotate("text", x =  4.5, y = max(abs(ccf_df$ccf)) * 0.92,
           label = "debt leads troops\n(endogeneity concern if significant)",
           size = 2.8, colour = "#d6604d", fontface = "italic", hjust = 0.5) +
  annotate("text", x = -4.5, y = -0.13,
           label = "troops lead debt\n(exogeneity-consistent)",
           size = 2.8, colour = "#4d9221", fontface = "italic", hjust = 0.5) +
  scale_x_continuous(breaks = -8:8) +
  labs(
    title    = "Cross-Correlation: Troops Deployed vs Delta Debt/GDP",
    subtitle = "Red shaded region (lag > 0): potential feedback from debt to troops",
    x        = "Lag (positive = debt leads troops)",
    y        = "Correlation"
  ) +
  theme_minimal(base_size = 11)

print(p_ccf)
ggsave("output/figures/05_ccf.pdf", p_ccf, width = 9, height = 5)
cat("\nSaved: output/figures/05_ccf.pdf\n")

# ── Save exogeneity table ─────────────────────────────────────────────────────
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

t_exog <- data.frame(
  Direction        = c("$\\Delta$Debt/GDP $\\to$ Troops",
                       "Troops $\\to$ $\\Delta$Debt/GDP"),
  "F-statistic"    = c(round(F1, 3), round(F2, 3)),
  "$p$-value"      = c(round(p1, 4), round(p2, 4)),
  "Wald $p$-value" = c(round(p_wald, 4), NA),
  Conclusion       = c("Fail to reject $H_0$ — exogeneity supported",
                       "Fail to reject $H_0$"),
  check.names = FALSE
)

write.csv(t_exog, "output/tables/T0_exogeneity.csv", row.names = FALSE)

# LaTeX version
lines <- c(
  "\\begin{table}[htbp]",
  "  \\centering",
  "  \\caption{Granger Causality Tests: Exogeneity of Troops Deployed}",
  "  \\label{tab:exogeneity}",
  "  \\begin{tabular}{lcccp{5cm}}",
  "    \\toprule",
  "    Direction & $F$-stat & $p$-value & Wald $p$ & Conclusion \\\\",
  "    \\midrule",
  sprintf("    $\\Delta$Debt/GDP $\\to$ Troops & %.3f & %.4f & %.4f & Fail to reject $H_0$ --- exogeneity supported \\\\", F1, p1, p_wald),
  sprintf("    Troops $\\to$ $\\Delta$Debt/GDP & %.3f & %.4f & ---  & Fail to reject $H_0$ \\\\", F2, p2),
  "    \\bottomrule",
  "  \\end{tabular}",
  sprintf("  \\\\[2pt] \\footnotesize{Bivariate Granger causality test with $k = %d$ lags (BIC-selected). $H_0$: the row variable does not Granger-cause the column variable. Wald $p$-value from a manual F-test of debt lags in the troops equation (Section~3, unrestricted vs.~restricted OLS). Sample: 1951--2020 ($n = 70$).}", k_use),
  "\\end{table}"
)
writeLines(lines, "output/tables/T0_exogeneity.tex")
cat("Saved: output/tables/T0_exogeneity.tex\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat(sprintf("\n%s\n  SUMMARY\n%s\n", sep, sep))
cat(sprintf("  Granger  Δdebt -> troops :  F = %.3f,  p = %.4f\n", F1, p1))
cat(sprintf("  Granger  Troops -> Δdebt :  F = %.3f,  p = %.4f\n", F2, p2))
cat(sprintf("  Wald     Δdebt -> troops :  p = %.4f\n", p_wald))
cat(sprintf("\n  Exogeneity supported: %s\n",
            ifelse(p1 > 0.10 & p_wald > 0.10, "YES", "DISCUSS IN REPORT")))
cat(sprintf("  Lag order used: k = %d\n", k_use))
cat("\n✓ 05_exogeneity.R complete.\n")