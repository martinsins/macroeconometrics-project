
# Plain-language summary
p1 <- gc_debt_to_troops[2, "Pr(>F)"]
p2 <- gc_troops_to_debt[2, "Pr(>F)"]

cat(sprintf("\n  *** VERDICT (k = %d) ***\n", k_use))
cat(sprintf("  Δdebt -> troops :  F = %.3f,  p = %.4f  =>  %s\n",
            gc_debt_to_troops[2, "F"],
            p1,
            ifelse(p1 > 0.10,
                   "FAIL TO REJECT H0 — no evidence debt Granger-causes troops",
                   ifelse(p1 > 0.05,
                          "MARGINAL (p<0.10) — weak evidence of feedback",
                          "REJECT H0 — debt appears to Granger-cause troops (!)"))))
cat(sprintf("  Troops -> Δdebt :  F = %.3f,  p = %.4f  =>  %s\n\n",
            gc_troops_to_debt[2, "F"],
            p2,
            ifelse(p2 > 0.10, "FAIL TO REJECT H0 — no predictive content",
                   ifelse(p2 > 0.05, "MARGINAL (p<0.10)",
                          "REJECT H0 — troops predicts debt changes"))))

# =============================================================================
# SECTION 3 — Robustness: Wald test on VAR troops equation
# =============================================================================
# We estimate the full bivariate VAR(k) and test whether lags of Δdebt are
# jointly zero in the troops equation.  This is equivalent to the Granger test
# but goes via the VAR object, which also lets us inspect all coefficients.
# =============================================================================

cat(sprintf("%s\n  SECTION 3 — VAR(%d) Wald test on troops equation\n%s\n\n",
            sep2, k_use, sep2))

var_fit <- VAR(biv_mat, p = k_use, type = "const")

# Causality test from vars package: tests Granger non-causality from
# "d_debt_dt" to "troops_dt" (H0: all lags of d_debt_dt = 0 in troops eq.)
wald_result <- causality(var_fit, cause = "d_debt_dt")

cat("  vars::causality — H0: Δdebt does NOT Granger-cause troops\n\n")
cat("  Granger test:\n")
print(wald_result$Granger)

cat("\n  Instantaneous causality test:\n")
print(wald_result$Instant)

# Coefficients of the troops equation
cat(sprintf("\n  Troops equation coefficients (VAR(%d)):\n", k_use))
print(round(coef(var_fit)$troops_dt, 4))

# =============================================================================
# SECTION 4 — Cross-correlation function (CCF)
# =============================================================================
# CCF(troops, Δdebt, lag h) = corr(troops_t, Δdebt_{t-h}).
# If troops_t correlates with future Δdebt (h > 0): troops -> debt (expected).
# If troops_t correlates with past  Δdebt (h < 0): debt -> troops (problem).
# =============================================================================

cat(sprintf("\n%s\n  SECTION 4 — Cross-correlation function\n%s\n\n", sep2, sep2))

ccf_vals <- ccf(as.numeric(x_dt), as.numeric(y_dt),
                lag.max = 8, plot = FALSE,
                ylab = "CCF(troops, Δdebt)")

cat("  CCF(troops_dt, Δdebt_dt) at lags -8 to +8:\n")
ccf_df <- data.frame(
  lag  = as.integer(ccf_vals$lag),
  ccf  = round(as.numeric(ccf_vals$acf), 4)
)
print(ccf_df)

ci <- qnorm(0.975) / sqrt(length(y_dt))
cat(sprintf("\n  95%% confidence band: ± %.4f\n", ci))
cat(sprintf("  Negative lags = troops leads debt (exogeneity-consistent)\n"))
cat(sprintf("  Positive lags = debt leads troops (feedback / endogeneity concern)\n\n"))

sig_pos <- ccf_df[ccf_df$lag > 0 & abs(ccf_df$ccf) > ci, ]
if (nrow(sig_pos) == 0) {
  cat("  => No significant CCF at positive lags (debt does not lead troops).\n")
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

p_ccf <- ggplot(ccf_df, aes(x = lag, y = ccf)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-ci, ci),
             linetype = "dashed", colour = "#2166ac", linewidth = 0.6) +
  geom_segment(aes(xend = lag, yend = 0), linewidth = 0.8,
               colour = ifelse(ccf_df$lag > 0 & abs(ccf_df$ccf) > ci,
                               "#d6604d", "grey40")) +
  geom_point(size = 2.5,
             colour = ifelse(ccf_df$lag > 0 & abs(ccf_df$ccf) > ci,
                             "#d6604d", "grey40")) +
  annotate("rect", xmin = 0.5, xmax = 8.5,
           ymin = -Inf, ymax = Inf,
           fill = "#d6604d", alpha = 0.04) +
  annotate("text", x = 4.5, y = max(abs(ccf_df$ccf)) * 0.95,
           label = "debt leads troops\n(endogeneity concern if significant)",
           size = 2.8, colour = "#d6604d", fontface = "italic",
           hjust = 0.5) +
  annotate("text", x = -4.5, y = max(abs(ccf_df$ccf)) * 0.95,
           label = "troops lead debt\n(exogeneity-consistent)",
           size = 2.8, colour = "#4d9221", fontface = "italic",
           hjust = 0.5) +
  scale_x_continuous(breaks = -8:8) +
  labs(
    title    = "Cross-Correlation: Troops Deployed vs Δ Debt/GDP",
    subtitle = "Bars in red shaded region (lag > 0): potential feedback from debt to troops",
    x        = "Lag (positive = debt leads troops)",
    y        = "Correlation"
  ) +
  theme_minimal(base_size = 11)

print(p_ccf)
ggsave("output/figures/05_ccf.pdf", p_ccf, width = 9, height = 5)
cat("Saved: output/figures/05_ccf.pdf\n")

cat(sprintf("\n%s\n  SUMMARY\n%s\n", sep, sep))
cat(sprintf("  Granger (debt -> troops):  p = %.4f\n", p1))
cat(sprintf("  Granger (troops -> debt):  p = %.4f\n", p2))
cat(sprintf("  Interpretation:\n"))
cat(sprintf("    If p(debt->troops) > 0.10: exogeneity assumption is supported.\n"))
cat(sprintf("    Report this alongside the main results as a specification check.\n"))

cat("\n✓ 05_exogeneity.R complete.\n")
