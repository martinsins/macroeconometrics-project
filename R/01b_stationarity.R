# =============================================================================
# 01b_stationarity.R
# Military Deployment & Public Debt — Trend Stationarity & Detrending
# =============================================================================
#
# For each key series (debt_gdp, d_debt_gdp, troops) we:
#   (1) Plot time series twice:
#         Panel A — series + horizontal mean line  (visual: constant stationarity?)
#         Panel B — series + OLS linear trend line (visual: trend stationarity?)
#   (2) Test stationarity AROUND A CONSTANT
#         ADF with drift  (H0: unit root; test eq. includes intercept only)
#         KPSS Level      (H0: stationary around constant)
#   (3) Test stationarity AROUND A LINEAR TREND
#         ADF with trend  (H0: unit root; test eq. includes intercept + trend)
#         KPSS Trend      (H0: stationary around linear trend)
#
# Decision rule (prof guidance):
#   Stationary around constant      => use as-is in ARMAX
#   Trend-stationary (not constant) => DETREND: OLS on time, keep residuals
#   Neither                         => ambiguous, flag for manual review
#
# Outputs:
#   output/figures/01b_all_trend.pdf      — combined 3×2 visual panels
#   output/figures/01b_detrended.pdf      — original vs detrended (if detrended)
#   data/processed/us_ts_detrended.csv   — data with *_dt columns added
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

# ── 0. Load data ──────────────────────────────────────────────────────────────
df       <- read.csv("data/processed/us_ts_data.csv")
df_clean <- df |> filter(!is.na(d_debt_gdp))

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

cat("Sample: debt_gdp =", nrow(df), "obs (",
    min(df$year), "-", max(df$year), ")\n")
cat("Sample: d_debt_gdp / troops =", nrow(df_clean), "obs (",
    min(df_clean$year), "-", max(df_clean$year), ")\n\n")

# =============================================================================
# SECTION 1 — Two-panel visualisation
# =============================================================================
# Panel A: series + mean dashed line   → "is the series flat around its mean?"
# Panel B: series + OLS trend line     → "does a linear trend fit the series?"
# =============================================================================

make_trend_panels <- function(df_in, varname, label, colour) {

  t_num  <- seq_len(nrow(df_in))
  y_vals <- df_in[[varname]]
  lm_fit <- lm(y_vals ~ t_num)

  df_plot <- data.frame(
    year  = df_in$year,
    y     = y_vals,
    trend = fitted(lm_fit),
    ymean = mean(y_vals)
  )

  yrange <- diff(range(y_vals))
  b <- coef(lm_fit)[2]
  a <- coef(lm_fit)[1]
  eq_lab <- sprintf("y = %.3f + (%.4f)×t", a, b)   # ×

  # Panel A — stationarity around constant
  pA <- ggplot(df_plot, aes(x = year, y = y)) +
    geom_line(colour = colour, linewidth = 0.8) +
    geom_hline(aes(yintercept = ymean),
               colour = "grey30", linetype = "dashed", linewidth = 0.7) +
    annotate("text",
             x = min(df_plot$year) + 1,
             y = df_plot$ymean[1] + 0.07 * yrange,
             label = sprintf("mean = %.3f", mean(y_vals)),
             hjust = 0, size = 3.2, colour = "grey30") +
    labs(subtitle = "(A) Around constant (mean)",
         x = NULL, y = label) +
    theme_minimal(base_size = 11)

  # Panel B — stationarity around linear trend
  pB <- ggplot(df_plot, aes(x = year, y = y)) +
    geom_line(colour = colour, linewidth = 0.8) +
    geom_line(aes(y = trend),
              colour = "grey30", linetype = "dashed", linewidth = 0.7) +
    annotate("text",
             x = min(df_plot$year) + 1,
             y = df_plot$trend[1] + 0.08 * yrange,
             label = eq_lab,
             hjust = 0, size = 3.2, colour = "grey30") +
    labs(subtitle = "(B) Around linear trend",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11)

  (pA | pB) +
    plot_annotation(
      title = label,
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )
}

p_debt   <- make_trend_panels(df,       "debt_gdp",   "Debt/GDP (level)",              "#2166ac")
p_ddbt   <- make_trend_panels(df_clean, "d_debt_gdp", "Delta Debt/GDP (first diff.)",  "#d6604d")
p_troops <- make_trend_panels(df_clean, "troops",     "Troops Deployed (thousands)",   "#4d9221")

combo_vis <- p_debt / p_ddbt / p_troops
print(combo_vis)
ggsave("output/figures/01b_all_trend.pdf", combo_vis, width = 11, height = 13)
cat("Saved: output/figures/01b_all_trend.pdf\n\n")

# =============================================================================
# SECTION 2 — Formal stationarity tests
# =============================================================================
# ADF via urca::ur.df
#   type = "drift" : H0 = unit root, test eq. = {intercept, Δlags, lagged level}
#   type = "trend" : H0 = unit root, test eq. = {intercept, trend, Δlags, lagged level}
#   lags selected by AIC
#   Reject H0 if tau < critical value (left-tail test, larger negative = more rejection)
#
# KPSS via tseries::kpss.test
#   null = "Level" : H0 = stationary around constant
#   null = "Trend" : H0 = stationary around linear trend
#   Reject H0 if stat > critical value (right-tail test, small p = reject)
# =============================================================================

run_tests_full <- function(series, label) {

  sep <- strrep("-", 64)
  cat(sprintf("\n%s\n  Variable: %s\n%s\n", sep, label, sep))

  # ─── (A) Around a constant ─────────────────────────────────────────────────
  cat("\n  [A] STATIONARITY AROUND A CONSTANT\n\n")

  adf_d        <- ur.df(series, type = "drift", selectlags = "AIC")
  tau2         <- adf_d@teststat[1, "tau2"]
  cv2          <- adf_d@cval["tau2", ]
  adf_d_reject <- tau2 < cv2["5pct"]   # reject unit root = evidence of stationarity

  cat(sprintf("  ADF (drift):   tau2 = %7.4f  |  CV: 1%%=%.3f  5%%=%.3f  10%%=%.3f\n",
              tau2, cv2["1pct"], cv2["5pct"], cv2["10pct"]))
  cat(sprintf("  => %s\n",
              ifelse(adf_d_reject,
                     "REJECT H0 (unit root) at 5% -- stationary",
                     "FAIL TO REJECT H0 (unit root) at 5%")))

  kpss_l        <- kpss.test(series, null = "Level")
  kpss_l_reject <- kpss_l$p.value < 0.05   # reject stationarity

  cat(sprintf("\n  KPSS (Level):  stat  = %7.4f  |  p-value = %.4f\n",
              kpss_l$statistic, kpss_l$p.value))
  cat(sprintf("  => %s\n",
              ifelse(kpss_l_reject,
                     "REJECT H0 (stationarity) at 5% -- NOT stationary around constant",
                     "FAIL TO REJECT H0 at 5% -- stationary around constant")))

  # ─── (B) Around a linear trend ─────────────────────────────────────────────
  cat("\n  [B] STATIONARITY AROUND A LINEAR TREND\n\n")

  adf_t        <- ur.df(series, type = "trend", selectlags = "AIC")
  tau3         <- adf_t@teststat[1, "tau3"]
  cv3          <- adf_t@cval["tau3", ]
  adf_t_reject <- tau3 < cv3["5pct"]

  cat(sprintf("  ADF (trend):   tau3 = %7.4f  |  CV: 1%%=%.3f  5%%=%.3f  10%%=%.3f\n",
              tau3, cv3["1pct"], cv3["5pct"], cv3["10pct"]))
  cat(sprintf("  => %s\n",
              ifelse(adf_t_reject,
                     "REJECT H0 (unit root) at 5% -- trend-stationary",
                     "FAIL TO REJECT H0 (unit root) at 5%")))

  kpss_t        <- kpss.test(series, null = "Trend")
  kpss_t_reject <- kpss_t$p.value < 0.05

  cat(sprintf("\n  KPSS (Trend):  stat  = %7.4f  |  p-value = %.4f\n",
              kpss_t$statistic, kpss_t$p.value))
  cat(sprintf("  => %s\n",
              ifelse(kpss_t_reject,
                     "REJECT H0 (stationarity around trend) at 5%",
                     "FAIL TO REJECT H0 at 5% -- stationary around linear trend")))

  # ─── Verdict ───────────────────────────────────────────────────────────────
  stat_const <- adf_d_reject & !kpss_l_reject
  stat_trend <- adf_t_reject & !kpss_t_reject

  cat("\n  *** VERDICT ***\n")
  if (stat_const) {
    cat("  => STATIONARY AROUND CONSTANT. Use as-is in ARMAX.\n")
  } else if (stat_trend) {
    cat("  => Trend-stationary (NOT stationary around constant).\n")
    cat("  => ACTION: DETREND -- regress on time, keep OLS residuals.\n")
  } else if (!adf_d_reject && !kpss_l_reject) {
    cat("  => Tests conflict (unit root in ADF, stationary in KPSS).\n")
    cat("  => Likely near-unit-root or structural break. Review visually.\n")
  } else {
    cat("  => NOT stationary around constant or trend. Consider differencing.\n")
    cat("  => If ADF (trend) is marginal, detrending may still be advisable.\n")
  }

  invisible(list(
    label        = label,
    stat_const   = stat_const,
    stat_trend   = stat_trend,
    adf_d_tau    = tau2,  adf_d_cv5 = cv2["5pct"],  adf_d_reject = adf_d_reject,
    adf_t_tau    = tau3,  adf_t_cv5 = cv3["5pct"],  adf_t_reject = adf_t_reject,
    kpss_l_p     = kpss_l$p.value,  kpss_l_reject = kpss_l_reject,
    kpss_t_p     = kpss_t$p.value,  kpss_t_reject = kpss_t_reject
  ))
}

# Build ts objects
ts_debt   <- ts(df$debt_gdp,         start = 1950, frequency = 1)
ts_ddbt   <- ts(df_clean$d_debt_gdp, start = 1951, frequency = 1)
ts_troops <- ts(df_clean$troops,     start = 1951, frequency = 1)

cat(strrep("=", 64), "\n")
cat("  STATIONARITY TESTS\n")
cat(strrep("=", 64), "\n")

res_debt   <- run_tests_full(ts_debt,   "Debt/GDP (level)")
res_ddbt   <- run_tests_full(ts_ddbt,   "Delta Debt/GDP (first diff.)")
res_troops <- run_tests_full(ts_troops, "Troops Deployed (thousands)")

# =============================================================================
# SECTION 3 — Summary table
# =============================================================================

cat(sprintf("\n%s\n  SUMMARY TABLE\n%s\n\n", strrep("=", 64), strrep("=", 64)))

summary_tbl <- data.frame(
  Variable      = c("Debt/GDP (level)", "Delta Debt/GDP", "Troops"),
  ADF_drift     = sprintf("%.3f [cv5=%.3f]",
                          c(res_debt$adf_d_tau, res_ddbt$adf_d_tau, res_troops$adf_d_tau),
                          c(res_debt$adf_d_cv5, res_ddbt$adf_d_cv5, res_troops$adf_d_cv5)),
  KPSS_Level_p  = round(c(res_debt$kpss_l_p, res_ddbt$kpss_l_p, res_troops$kpss_l_p), 4),
  ADF_trend     = sprintf("%.3f [cv5=%.3f]",
                          c(res_debt$adf_t_tau, res_ddbt$adf_t_tau, res_troops$adf_t_tau),
                          c(res_debt$adf_t_cv5, res_ddbt$adf_t_cv5, res_troops$adf_t_cv5)),
  KPSS_Trend_p  = round(c(res_debt$kpss_t_p, res_ddbt$kpss_t_p, res_troops$kpss_t_p), 4),
  Stat_Const    = c(res_debt$stat_const, res_ddbt$stat_const, res_troops$stat_const),
  Stat_Trend    = c(res_debt$stat_trend, res_ddbt$stat_trend, res_troops$stat_trend),
  Action        = c(
    ifelse(res_debt$stat_const,   "Use as-is", ifelse(res_debt$stat_trend,   "Detrend", "Review")),
    ifelse(res_ddbt$stat_const,   "Use as-is", ifelse(res_ddbt$stat_trend,   "Detrend", "Review")),
    ifelse(res_troops$stat_const, "Use as-is", ifelse(res_troops$stat_trend, "Detrend", "Review"))
  ),
  stringsAsFactors = FALSE
)

print(summary_tbl, row.names = FALSE)

# =============================================================================
# SECTION 4 — Detrend where indicated and save
# =============================================================================

cat(sprintf("\n%s\n  DETRENDING\n%s\n", strrep("=", 64), strrep("=", 64)))

detrend_ols <- function(y_vals, varname) {
  t_num  <- seq_len(length(y_vals))
  lm_obj <- lm(y_vals ~ t_num)
  s      <- summary(lm_obj)
  cat(sprintf("\n  Detrending %-15s  R2=%.4f  slope=%.5f  p(slope)=%.4f\n",
              varname,
              s$r.squared,
              coef(lm_obj)[2],
              s$coefficients[2, 4]))
  residuals(lm_obj)
}

df_out <- df_clean

# d_debt_gdp
if (!res_ddbt$stat_const && res_ddbt$stat_trend) {
  cat("\n[d_debt_gdp] Trend-stationary -> detrending\n")
  df_out$d_debt_gdp_dt <- detrend_ols(df_clean$d_debt_gdp, "d_debt_gdp")
} else {
  cat("\n[d_debt_gdp] Stationary around constant -> no detrending needed\n")
  df_out$d_debt_gdp_dt <- df_clean$d_debt_gdp
}

# troops
if (!res_troops$stat_const && res_troops$stat_trend) {
  cat("\n[troops] Trend-stationary -> detrending\n")
  df_out$troops_dt <- detrend_ols(df_clean$troops, "troops")
} else {
  cat("\n[troops] Stationary around constant -> no detrending needed\n")
  df_out$troops_dt <- df_clean$troops
}

write.csv(df_out, "data/processed/us_ts_detrended.csv", row.names = FALSE)
cat(sprintf("\nSaved: data/processed/us_ts_detrended.csv\n"))

cat("\nColumns in output:\n")
cat(sprintf("  d_debt_gdp_dt : %s\n",
            ifelse(isTRUE(all.equal(df_out$d_debt_gdp_dt, df_clean$d_debt_gdp)),
                   "identical to d_debt_gdp (no detrend)", "OLS residuals (detrended)")))
cat(sprintf("  troops_dt     : %s\n",
            ifelse(isTRUE(all.equal(df_out$troops_dt, df_clean$troops)),
                   "identical to troops (no detrend)", "OLS residuals (detrended)")))

# =============================================================================
# SECTION 5 — Plot original vs detrended (only if detrending was applied)
# =============================================================================

detrend_y <- !isTRUE(all.equal(df_out$d_debt_gdp_dt, df_clean$d_debt_gdp))
detrend_x <- !isTRUE(all.equal(df_out$troops_dt,     df_clean$troops))

if (detrend_y || detrend_x) {

  plots_dt <- list()

  if (detrend_y) {
    df_dt_y <- data.frame(
      year      = df_out$year,
      Original  = df_out$d_debt_gdp,
      Detrended = df_out$d_debt_gdp_dt
    ) |> pivot_longer(-year, names_to = "Series", values_to = "value")

    plots_dt[["y"]] <- ggplot(df_dt_y, aes(x = year, y = value, colour = Series)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      scale_colour_manual(values = c(Original = "#d6604d", Detrended = "#2166ac")) +
      labs(title = "Delta Debt/GDP: original vs detrended",
           x = NULL, y = "Delta Debt/GDP", colour = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
  }

  if (detrend_x) {
    df_dt_x <- data.frame(
      year      = df_out$year,
      Original  = df_out$troops,
      Detrended = df_out$troops_dt
    ) |> pivot_longer(-year, names_to = "Series", values_to = "value")

    plots_dt[["x"]] <- ggplot(df_dt_x, aes(x = year, y = value, colour = Series)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      scale_colour_manual(values = c(Original = "#4d9221", Detrended = "#762a83")) +
      labs(title = "Troops Deployed: original vs detrended",
           x = NULL, y = "Troops (thousands)", colour = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
  }

  combo_dt <- Reduce(`/`, plots_dt)
  print(combo_dt)
  ggsave("output/figures/01b_detrended.pdf", combo_dt,
         width = 9, height = 5 * length(plots_dt))
  cat("Saved: output/figures/01b_detrended.pdf\n")

} else {
  cat("\nNo detrending applied to any series -- 01b_detrended.pdf not produced.\n")
}

cat("\n✓ 01b_stationarity.R complete.\n")
