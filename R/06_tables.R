# =============================================================================
# 06_tables.R
# Military Deployment & Public Debt — Report Tables
# =============================================================================
#
# Generates all tables needed for the final report, in two formats:
#   - LaTeX (.tex)  : paste directly into the report
#   - CSV   (.csv)  : inspect / use in R Markdown via knitr::kable
#
# Tables produced:
#   T1 — Descriptive statistics
#   T2 — Stationarity test summary (ADF + KPSS)
#   T3 — Model selection grid (AIC / BIC over ARMA(p,q) orders)
#   T4 — Main ARMA(1,0)-X coefficient table
#   T5 — IRF point estimates + 95% CI (selected horizons)
#
# Run after: 01b_stationarity.R, 02_armax.R, 03_irf.R
# Requires: output/armax_fit.rds, data/processed/us_ts_detrended.csv
# =============================================================================

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

sep <- strrep("=", 64)

# ── Helper: write LaTeX table ─────────────────────────────────────────────────
write_latex <- function(df, file, caption, label,
                        align = NULL, hline_after = NULL, footnote = NULL) {
  
  nc   <- ncol(df)
  cols <- if (is.null(align)) paste0(c("l", rep("r", nc - 1)), collapse = "") else align
  
  lines <- c(
    "\\begin{table}[htbp]",
    "  \\centering",
    sprintf("  \\caption{%s}", caption),
    sprintf("  \\label{%s}", label),
    sprintf("  \\begin{tabular}{%s}", cols),
    "    \\toprule",
    paste0("    ", paste(names(df), collapse = " & "), " \\\\"),
    "    \\midrule"
  )
  
  for (i in seq_len(nrow(df))) {
    row_str <- paste(as.character(unlist(df[i, ])), collapse = " & ")
    lines   <- c(lines, paste0("    ", row_str, " \\\\"))
    if (!is.null(hline_after) && i %in% hline_after)
      lines <- c(lines, "    \\midrule")
  }
  
  lines <- c(lines, "    \\bottomrule")
  
  if (!is.null(footnote)) {
    lines <- c(lines,
               "  \\end{tabular}",
               sprintf("  \\\\[2pt] \\footnotesize{%s}", footnote),
               "\\end{table}")
  } else {
    lines <- c(lines, "  \\end{tabular}", "\\end{table}")
  }
  
  writeLines(lines, file)
  cat(sprintf("  Saved: %s\n", file))
}

# =============================================================================
# T1 — Descriptive statistics
# =============================================================================
cat(sprintf("\n%s\n  T1 — Descriptive statistics\n%s\n", sep, sep))

df <- read.csv("data/processed/us_ts_detrended.csv")

desc_stat <- function(x, name) {
  data.frame(
    Variable = name,
    N        = length(x),
    Mean     = round(mean(x),   4),
    SD       = round(sd(x),     4),
    Min      = round(min(x),    4),
    Median   = round(median(x), 4),
    Max      = round(max(x),    4)
  )
}

t1 <- rbind(
  desc_stat(df$d_debt_gdp,    "$\\Delta$ Debt/GDP (raw)"),
  desc_stat(df$d_debt_gdp_dt, "$\\Delta$ Debt/GDP (detrended)"),
  desc_stat(df$troops,        "Troops deployed (thousands, raw)"),
  desc_stat(df$troops_dt,     "Troops deployed (thousands, detrended)")
)

print(t1, row.names = FALSE)
write.csv(t1, "output/tables/T1_descriptive.csv", row.names = FALSE)
write_latex(
  t1,
  "output/tables/T1_descriptive.tex",
  caption   = "Descriptive Statistics",
  label     = "tab:descriptive",
  align     = "lrrrrrr",
  footnote  = paste0(
    "Annual data 1951--", max(df$year), ". ",
    "Detrended series are OLS residuals from regressing each variable on a ",
    "linear time trend (see Section~\\ref{sec:stationarity})."
  )
)

# =============================================================================
# T2 — Stationarity tests
# =============================================================================
cat(sprintf("\n%s\n  T2 — Stationarity test summary\n%s\n", sep, sep))

# Re-run tests to get exact numbers (mirrors 01b_stationarity.R logic)
run_tests_silent <- function(series) {
  adf_d  <- urca::ur.df(series, type = "drift", selectlags = "AIC")
  adf_t  <- urca::ur.df(series, type = "trend", selectlags = "AIC")
  kpss_l <- tseries::kpss.test(series, null = "Level")
  kpss_t <- tseries::kpss.test(series, null = "Trend")
  list(
    tau2    = adf_d@teststat[1, "tau2"],
    cv2_5   = adf_d@cval["tau2", "5pct"],
    kpss_l  = kpss_l$statistic,
    kpss_lp = kpss_l$p.value,
    tau3    = adf_t@teststat[1, "tau3"],
    cv3_5   = adf_t@cval["tau3", "5pct"],
    kpss_t  = kpss_t$statistic,
    kpss_tp = kpss_t$p.value
  )
}

ts_debt   <- ts(df$d_debt_gdp,    start = min(df$year), frequency = 1)
ts_ddbt   <- ts(df$d_debt_gdp_dt, start = min(df$year), frequency = 1)
ts_troops <- ts(df$troops,        start = min(df$year), frequency = 1)
ts_tdt    <- ts(df$troops_dt,     start = min(df$year), frequency = 1)

r1 <- run_tests_silent(ts_debt)
r2 <- run_tests_silent(ts_ddbt)
r3 <- run_tests_silent(ts_troops)
r4 <- run_tests_silent(ts_tdt)

stars <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                                                    ifelse(p < 0.10, "*", "")))

fmt_adf  <- function(tau, cv5) sprintf("%.3f [%.3f]", tau, cv5)
fmt_kpss <- function(stat, p)   sprintf("%.3f%s", stat, stars(p))

t2 <- data.frame(
  Variable = c("$\\Delta$ Debt/GDP (raw)",
               "$\\Delta$ Debt/GDP (detrended)",
               "Troops (raw)",
               "Troops (detrended)"),
  "ADF-drift $\\tau_2$ [CV5\\%]"  = c(fmt_adf(r1$tau2, r1$cv2_5),
                                      fmt_adf(r2$tau2, r2$cv2_5),
                                      fmt_adf(r3$tau2, r3$cv2_5),
                                      fmt_adf(r4$tau2, r4$cv2_5)),
  "KPSS Level"  = c(fmt_kpss(r1$kpss_l, r1$kpss_lp),
                    fmt_kpss(r2$kpss_l, r2$kpss_lp),
                    fmt_kpss(r3$kpss_l, r3$kpss_lp),
                    fmt_kpss(r4$kpss_l, r4$kpss_lp)),
  "ADF-trend $\\tau_3$ [CV5\\%]"  = c(fmt_adf(r1$tau3, r1$cv3_5),
                                      fmt_adf(r2$tau3, r2$cv3_5),
                                      fmt_adf(r3$tau3, r3$cv3_5),
                                      fmt_adf(r4$tau3, r4$cv3_5)),
  "KPSS Trend"  = c(fmt_kpss(r1$kpss_t, r1$kpss_tp),
                    fmt_kpss(r2$kpss_t, r2$kpss_tp),
                    fmt_kpss(r3$kpss_t, r3$kpss_tp),
                    fmt_kpss(r4$kpss_t, r4$kpss_tp)),
  Conclusion = c("Trend-stationary $\\Rightarrow$ detrend",
                 "Stationary around constant",
                 "Trend-stationary $\\Rightarrow$ detrend",
                 "Stationary around constant"),
  check.names = FALSE
)

print(t2, row.names = FALSE)
write.csv(t2, "output/tables/T2_stationarity.csv", row.names = FALSE)
write_latex(
  t2,
  "output/tables/T2_stationarity.tex",
  caption  = "Stationarity Tests: ADF and KPSS",
  label    = "tab:stationarity",
  align    = "lccccl",
  footnote = paste0(
    "ADF: Augmented Dickey--Fuller test. $H_0$: unit root. ",
    "Lag order selected by AIC. Brackets show 5\\% critical value; ",
    "reject $H_0$ if $\\tau < $ CV. ",
    "KPSS: Kwiatkowski--Phillips--Schmidt--Shin test. $H_0$: stationarity. ",
    "Asterisks denote rejection of $H_0$: ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$. ",
    "A series is classified as trend-stationary if ADF-trend rejects the unit root ",
    "and KPSS-trend does not reject stationarity."
  )
)

# =============================================================================
# T3 — Model selection grid (AIC / BIC)
# =============================================================================
cat(sprintf("\n%s\n  T3 — Model selection grid\n%s\n", sep, sep))

saved  <- readRDS("output/armax_fit.rds")
fit    <- saved$fit
p_sel  <- saved$p
q_sel  <- saved$q
df_est <- saved$df
xreg   <- saved$xreg

y_ts <- ts(df_est$d_debt_gdp_dt, start = min(df_est$year), frequency = 1)

grid <- expand.grid(p = 0:3, q = 0:3) |> dplyr::filter(!(p == 0 & q == 0))

aic_vals <- mapply(function(p, q) {
  m <- tryCatch(
    forecast::Arima(y_ts, order = c(p, 0, q), xreg = xreg, include.mean = TRUE),
    error = function(e) NULL)
  if (is.null(m)) NA else round(AIC(m), 2)
}, grid$p, grid$q)

bic_vals <- mapply(function(p, q) {
  m <- tryCatch(
    forecast::Arima(y_ts, order = c(p, 0, q), xreg = xreg, include.mean = TRUE),
    error = function(e) NULL)
  if (is.null(m)) NA else round(BIC(m), 2)
}, grid$p, grid$q)

grid$AIC <- aic_vals
grid$BIC <- bic_vals

# Format: bold the selected cell (LaTeX \textbf{})
t3 <- grid |>
  dplyr::arrange(p, q) |>
  dplyr::mutate(
    Spec = sprintf("ARMA(%d,%d)-X", p, q),
    AIC  = ifelse(p == which.min(aic_vals) %% 4 & q == ceiling(which.min(aic_vals) / 4),
                  sprintf("\\textbf{%.2f}", AIC), sprintf("%.2f", AIC)),
    BIC  = ifelse(BIC == min(BIC, na.rm = TRUE),
                  sprintf("\\textbf{%.2f}", BIC), sprintf("%.2f", BIC))
  ) |>
  dplyr::select(Spec, AIC, BIC)

# Mark selected models clearly
min_aic_spec <- sprintf("ARMA(%d,%d)-X", grid$p[which.min(aic_vals)], grid$q[which.min(aic_vals)])
min_bic_spec <- sprintf("ARMA(%d,%d)-X", grid$p[which.min(bic_vals)], grid$q[which.min(bic_vals)])

t3$AIC <- ifelse(t3$Spec == min_aic_spec,
                 paste0("\\textbf{", sub("\\\\textbf\\{|\\}", "", t3$AIC), "}"),
                 t3$AIC)
t3$BIC <- ifelse(t3$Spec == min_bic_spec,
                 paste0("\\textbf{", sub("\\\\textbf\\{|\\}", "", t3$BIC), "}"),
                 t3$BIC)

grid_print <- grid |>
  dplyr::arrange(p, q) |>
  dplyr::mutate(Spec = sprintf("ARMA(%d,%d)-X", p, q)) |>
  dplyr::select(Spec, AIC, BIC)
print(grid_print, row.names = FALSE)
cat(sprintf("\n  AIC selects: %s\n", min_aic_spec))
cat(sprintf("  BIC selects: %s\n", min_bic_spec))

write.csv(grid_print, "output/tables/T3_model_selection.csv", row.names = FALSE)
write_latex(
  t3,
  "output/tables/T3_model_selection.tex",
  caption  = "ARMA-X Model Selection: AIC and BIC",
  label    = "tab:model_selection",
  align    = "lrr",
  footnote = paste0(
    "AIC and BIC computed over a grid of ARMA$(p,q)$-X specifications ",
    "($p, q \\in \\{0,\\ldots,3\\}$, excluding the $(0,0)$ case). ",
    "The exogenous regressors (troops, crisis dummies) are held fixed across all models. ",
    "Bold entries indicate the model selected by each criterion. ",
    "BIC is used for the main specification."
  )
)

# =============================================================================
# T4 — Main coefficient table
# =============================================================================
cat(sprintf("\n%s\n  T4 — Main ARMA(%d,%d)-X coefficient table\n%s\n",
            sep, p_sel, q_sel, sep))

cf   <- coef(fit)
se   <- sqrt(diag(vcov(fit)))
tstat <- cf / se
pval <- 2 * pnorm(-abs(tstat))

stars_p <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                                                      ifelse(p < 0.10, "*", "")))

# Clean parameter names for the report
name_map <- c(
  ar1      = "$\\phi_1$ (AR 1)",
  ma1      = "$\\theta_1$ (MA 1)",
  intercept= "$\\mu$ (intercept)",
  troops   = "$\\beta$ (troops, detrended)",
  d2009    = "$\\gamma_{2009}$ (GFC dummy)",
  d2020    = "$\\gamma_{2020}$ (COVID dummy)"
)

param_names <- names(cf)
display_names <- ifelse(param_names %in% names(name_map),
                        name_map[param_names], param_names)

t4 <- data.frame(
  Parameter = display_names,
  Estimate  = sprintf("%.5f", cf),
  "Std. Error" = sprintf("%.5f", se),
  "$t$-stat"   = sprintf("%.3f",  tstat),
  "$p$-value"  = paste0(sprintf("%.4f", pval), stars_p(pval)),
  check.names  = FALSE
)

print(t4, row.names = FALSE)

cat(sprintf("\n  Log-likelihood: %.4f\n", logLik(fit)))
cat(sprintf("  AIC: %.4f\n", AIC(fit)))
cat(sprintf("  BIC: %.4f\n", BIC(fit)))
cat(sprintf("  sigma^2: %.6f\n", fit$sigma2))
cat(sprintf("  n = %d\n", length(residuals(fit))))

# Add model fit rows
t4_full <- rbind(
  t4,
  data.frame(Parameter    = "---",
             Estimate     = "---",
             "Std. Error" = "---",
             "$t$-stat"   = "---",
             "$p$-value"  = "---",
             check.names  = FALSE),
  data.frame(Parameter    = "$\\sigma^2$",
             Estimate     = sprintf("%.6f", fit$sigma2),
             "Std. Error" = "",
             "$t$-stat"   = "",
             "$p$-value"  = "",
             check.names  = FALSE),
  data.frame(Parameter    = "Log-likelihood",
             Estimate     = sprintf("%.3f", as.numeric(logLik(fit))),
             "Std. Error" = "",
             "$t$-stat"   = "",
             "$p$-value"  = "",
             check.names  = FALSE),
  data.frame(Parameter    = "AIC / BIC",
             Estimate     = sprintf("%.2f / %.2f", AIC(fit), BIC(fit)),
             "Std. Error" = "",
             "$t$-stat"   = "",
             "$p$-value"  = "",
             check.names  = FALSE),
  data.frame(Parameter    = "Observations",
             Estimate     = as.character(length(residuals(fit))),
             "Std. Error" = "",
             "$t$-stat"   = "",
             "$p$-value"  = "",
             check.names  = FALSE)
)

write.csv(t4, "output/tables/T4_coefficients.csv", row.names = FALSE)
write_latex(
  t4_full,
  "output/tables/T4_coefficients.tex",
  caption     = sprintf("ARMA(%d,%d)-X Estimation Results", p_sel, q_sel),
  label       = "tab:armax_results",
  align       = "lrrrr",
  hline_after = nrow(t4),
  footnote    = paste0(
    "Dependent variable: $\\Delta$Debt/GDP (detrended). ",
    "Exogenous regressor: troops deployed abroad (thousands, detrended). ",
    "Crisis dummies absorb the 2009 GFC and 2020 COVID fiscal shocks. ",
    "Estimated by maximum likelihood. Standard errors from the ",
    "observed information matrix. ",
    "* $p<0.10$, ** $p<0.05$, *** $p<0.01$."
  )
)

# =============================================================================
# T5 — IRF table (selected horizons)
# =============================================================================
cat(sprintf("\n%s\n  T5 — IRF table\n%s\n", sep, sep))

SHOCK_SIZE <- 100

cf_irf  <- coef(fit)
V       <- vcov(fit)
p_i     <- p_sel
q_i     <- q_sel

ar_names <- if (p_i > 0) paste0("ar", seq_len(p_i)) else character(0)
ma_names <- if (q_i > 0) paste0("ma", seq_len(q_i)) else character(0)
phi      <- if (p_i > 0) cf_irf[ar_names] else numeric(0)
theta    <- if (q_i > 0) cf_irf[ma_names] else numeric(0)
beta_irf <- cf_irf["troops"]

H    <- 20
psi  <- c(1, if (p_i > 0 || q_i > 0) ARMAtoMA(ar = phi, ma = theta, lag.max = H) else rep(0, H))
irf  <- psi * beta_irf * SHOCK_SIZE
cirf <- cumsum(irf)

# Delta method SE (reuse logic from 03_irf.R)
irf_params <- c(ar_names, ma_names, "troops")
irf_params <- irf_params[irf_params %in% colnames(V)]
params0    <- cf_irf[irf_params]
V_sub      <- V[irf_params, irf_params]
eps        <- 1e-6
n_par      <- length(params0)
jacobian   <- matrix(0, nrow = H + 1, ncol = n_par)

get_psi <- function(phi, theta, H) {
  if (length(phi) == 0 && length(theta) == 0) return(c(1, rep(0, H)))
  c(1, ARMAtoMA(ar = phi, ma = theta, lag.max = H))
}

for (k in seq_len(n_par)) {
  pup <- params0; pup[k] <- params0[k] + eps
  pdn <- params0; pdn[k] <- params0[k] - eps
  phi_up <- if (length(ar_names) > 0) pup[ar_names] else numeric(0)
  phi_dn <- if (length(ar_names) > 0) pdn[ar_names] else numeric(0)
  th_up  <- if (length(ma_names) > 0) pup[ma_names] else numeric(0)
  th_dn  <- if (length(ma_names) > 0) pdn[ma_names] else numeric(0)
  irf_up <- get_psi(phi_up, th_up, H) * pup["troops"] * SHOCK_SIZE
  irf_dn <- get_psi(phi_dn, th_dn, H) * pdn["troops"] * SHOCK_SIZE
  jacobian[, k] <- (irf_up - irf_dn) / (2 * eps)
}

irf_se      <- sqrt(pmax(diag(jacobian %*% V_sub %*% t(jacobian)), 0))
jac_cum     <- apply(jacobian, 2, cumsum)
cirf_se     <- sqrt(pmax(diag(jac_cum %*% V_sub %*% t(jac_cum)), 0))

# Select horizons 0-5, 10, 15, 20 for the table
horizons <- c(0:5, 10, 15, 20)
t5 <- data.frame(
  "$h$ (years)"       = horizons,
  IRF                  = sprintf("%.5f", irf[horizons + 1]),
  "SE"                 = sprintf("%.5f", irf_se[horizons + 1]),
  "95\\% CI lower"     = sprintf("%.5f", irf[horizons + 1] - 1.96 * irf_se[horizons + 1]),
  "95\\% CI upper"     = sprintf("%.5f", irf[horizons + 1] + 1.96 * irf_se[horizons + 1]),
  "Cum. IRF"           = sprintf("%.5f", cirf[horizons + 1]),
  "Cum. SE"            = sprintf("%.5f", cirf_se[horizons + 1]),
  check.names = FALSE
)

print(t5, row.names = FALSE)
cat(sprintf("\n  Long-run multiplier (sum over %d years): %.5f\n", H, sum(irf)))

write.csv(t5, "output/tables/T5_irf.csv", row.names = FALSE)
write_latex(
  t5,
  "output/tables/T5_irf.tex",
  caption  = "Impulse Response Function: Effect of $+100$k Troops on $\\Delta$ Debt/GDP",
  label    = "tab:irf",
  align    = "crrrrrr",
  footnote = paste0(
    "IRF computed from the MA($\\infty$) representation of the ARMA(",
    p_sel, ",", q_sel, ") part, scaled by $\\hat{\\beta}$ and the shock size. ",
    "Shock: a one-time increase of 100,000 troops. ",
    "Standard errors and confidence intervals via the delta method ",
    "(numerical Jacobian of the IRF with respect to all estimated parameters). ",
    sprintf("Long-run multiplier (sum over %d years): %.5f.", H, sum(irf))
  )
)

# =============================================================================
# Summary
# =============================================================================
cat(sprintf("\n%s\n  ALL TABLES SAVED\n%s\n", sep, sep))
cat("  output/tables/T1_descriptive.csv/.tex\n")
cat("  output/tables/T2_stationarity.csv/.tex\n")
cat("  output/tables/T3_model_selection.csv/.tex\n")
cat("  output/tables/T4_coefficients.csv/.tex\n")
cat("  output/tables/T5_irf.csv/.tex\n")
cat("\n✓ 06_tables.R complete.\n")