# =============================================================================
# 02_armax.R
# Military Deployment & Public Debt — ARMA-X Estimation
# =============================================================================
#
# Main spec: Delta_debt_dt_t = mu + AR(p) + MA(q) + beta*troops_dt_t
#                            + gamma1*d2009_t + gamma2*d2020_t + eps_t
#
# Uses detrended series from 01b_stationarity.R (us_ts_detrended.csv):
#   d_debt_gdp_dt : Δdebt/GDP, detrended if 01b found trend-stationarity
#   troops_dt     : troops deployed, detrended if 01b found trend-stationarity
#
# Crisis dummies absorb the GFC (2009) and COVID (2020) fiscal shocks.
# =============================================================================


# ── 0. Load data ──────────────────────────────────────────────────────────────
# us_ts_detrended.csv is produced by 01b_stationarity.R; already filtered
# (no NA rows) and contains _dt columns (detrended series where applicable).
df_clean <- read.csv("data/processed/us_ts_detrended.csv")

y <- ts(df_clean$d_debt_gdp_dt, start = 1951, frequency = 1)

# ── 1. Build xreg matrix ──────────────────────────────────────────────────────
xreg <- cbind(
  troops = df_clean$troops_dt,
  d2009  = as.integer(df_clean$year == 2009),
  d2020  = as.integer(df_clean$year == 2020)
)

cat("=== xreg matrix (first and last 5 rows) ===\n")
print(head(xreg, 5))
print(tail(xreg, 5))

# ── 2. Grid search over (p, q) — AIC / BIC ───────────────────────────────────
p_max <- 3
q_max <- 3

results <- expand.grid(p = 0:p_max, q = 0:q_max) |>
  filter(!(p == 0 & q == 0))

fit_armax <- function(p, q) {
  tryCatch(
    Arima(y, order = c(p, 0, q), xreg = xreg, include.mean = TRUE),
    error = function(e) NULL
  )
}

results$AIC <- mapply(function(p, q) {
  m <- fit_armax(p, q)
  if (is.null(m)) NA else AIC(m)
}, results$p, results$q)

results$BIC <- mapply(function(p, q) {
  m <- fit_armax(p, q)
  if (is.null(m)) NA else BIC(m)
}, results$p, results$q)

cat("\n=== Model Selection Grid (AIC / BIC) ===\n")
print(results |> arrange(AIC))

best_aic <- results[which.min(results$AIC), ]
best_bic <- results[which.min(results$BIC), ]

cat(sprintf("\nBest by AIC: ARMA(%d,%d)  AIC = %.2f\n", best_aic$p, best_aic$q, best_aic$AIC))
cat(sprintf("Best by BIC: ARMA(%d,%d)  BIC = %.2f\n", best_bic$p, best_bic$q, best_bic$BIC))

# If AIC and BIC select different orders, compare beta(troops) across both
if (!(best_aic$p == best_bic$p && best_aic$q == best_bic$q)) {
  fit_aic <- Arima(y, order = c(best_aic$p, 0, best_aic$q), xreg = xreg, include.mean = TRUE)
  fit_bic <- Arima(y, order = c(best_bic$p, 0, best_bic$q), xreg = xreg, include.mean = TRUE)
  b_aic   <- coef(fit_aic)["troops"]; se_aic <- sqrt(vcov(fit_aic)["troops", "troops"])
  b_bic   <- coef(fit_bic)["troops"]; se_bic <- sqrt(vcov(fit_bic)["troops", "troops"])
  cat(sprintf("\nbeta(troops): AIC ARMA(%d,%d) = %.5f (SE=%.5f) | BIC ARMA(%d,%d) = %.5f (SE=%.5f)\n",
              best_aic$p, best_aic$q, b_aic, se_aic,
              best_bic$p, best_bic$q, b_bic, se_bic))
}

# ── 3. Estimate selected model (BIC — penalises more on small samples) ────────
p_sel <- best_bic$p
q_sel <- best_bic$q

cat(sprintf("\n=== Estimating ARMA(%d,%d)-X ===\n", p_sel, q_sel))

fit <- Arima(y, order = c(p_sel, 0, q_sel), xreg = xreg, include.mean = TRUE)
print(summary(fit))

# ── 4. Coefficient table ──────────────────────────────────────────────────────
coef_est <- coef(fit)
coef_se  <- sqrt(diag(vcov(fit)))
coef_t   <- coef_est / coef_se
coef_p   <- 2 * pnorm(-abs(coef_t))

coef_table <- data.frame(
  estimate = round(coef_est, 5),
  std_err  = round(coef_se,  5),
  t_stat   = round(coef_t,   3),
  p_value  = round(coef_p,   4)
)

cat("\n=== Coefficient Table ===\n")
print(coef_table)

# ── 5. Residual diagnostics ───────────────────────────────────────────────────
res <- residuals(fit)

cat("\n=== Ljung-Box Test on Residuals ===\n")
lb <- Box.test(res, lag = 10, type = "Ljung-Box", fitdf = p_sel + q_sel)
cat(sprintf("  Ljung-Box (lag 10): stat = %.3f,  p = %.4f\n", lb$statistic, lb$p.value))
lb20 <- Box.test(res, lag = 20, type = "Ljung-Box", fitdf = p_sel + q_sel)
cat(sprintf("  Ljung-Box (lag 20): stat = %.3f,  p = %.4f\n", lb20$statistic, lb20$p.value))

res_df <- data.frame(year = df_clean$year, residual = as.numeric(res))

p_res <- ggplot(res_df, aes(x = year, y = residual)) +
  geom_line(colour = "#555555") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = c(-2, 2) * sd(res), linetype = "dotted", colour = "#d6604d") +
  labs(title = sprintf("Residuals — ARMA(%d,%d)-X", p_sel, q_sel),
       x = NULL, y = "Residual") +
  theme_minimal(base_size = 11)

p_acf_res <- forecast::ggAcf(res, lag.max = 20) +
  labs(title = "ACF of residuals") +
  theme_minimal(base_size = 11)

combo_diag <- p_res / p_acf_res
print(combo_diag)
ggsave("output/figures/02_residuals.pdf", combo_diag, width = 9, height = 7)

# ── 6. Fitted vs actual ───────────────────────────────────────────────────────
fit_df <- data.frame(
  year    = df_clean$year,
  actual  = as.numeric(y),
  fitted  = as.numeric(fitted(fit))
)

p_fit <- ggplot(fit_df, aes(x = year)) +
  geom_line(aes(y = actual), colour = "#d6604d", linewidth = 0.8) +
  geom_line(aes(y = fitted), colour = "#2166ac", linewidth = 0.7, linetype = "dashed") +
  labs(title = sprintf("Actual vs Fitted — ARMA(%d,%d)-X  (red=actual, blue=fitted)", p_sel, q_sel),
       x = NULL, y = "Delta Debt / GDP") +
  theme_minimal(base_size = 11)

print(p_fit)
ggsave("output/figures/02_fitted.pdf", p_fit, width = 9, height = 5)

# ── 7. Save model object ──────────────────────────────────────────────────────
dir.create("output", showWarnings = FALSE)
saveRDS(list(fit = fit, p = p_sel, q = q_sel, xreg = xreg, df = df_clean),
        "output/armax_fit.rds")

cat("\n✓ 02_armax.R complete. Model saved to output/armax_fit.rds\n")
