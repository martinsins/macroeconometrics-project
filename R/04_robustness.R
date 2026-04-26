# =============================================================================
# 04_robustness.R
# Military Deployment & Public Debt — Robustness Checks
# =============================================================================
#
# Loads the main ARMA(p,q)-X model estimated in 02_armax.R and re-estimates
# it under alternative specifications to assess coefficient stability.
# =============================================================================

library(forecast)
library(dplyr)

# ── 0. Load main model ────────────────────────────────────────────────────────
saved   <- readRDS("output/armax_fit.rds")
fit     <- saved$fit
p_sel   <- saved$p
q_sel   <- saved$q
df_clean <- saved$df

y <- ts(df_clean$d_debt_gdp, start = 1951, frequency = 1)

cat(sprintf("Loaded ARMA(%d,%d)-X from output/armax_fit.rds\n", p_sel, q_sel))

# ── 1. Detrended troops ───────────────────────────────────────────────────────
rtrend    <- time(y) - mean(time(y))
troops_dt <- residuals(lm(df_clean$troops ~ as.numeric(rtrend)))

cat(sprintf("\n[Check 1] Troops detrending: R2 = %.3f\n",
            summary(lm(df_clean$troops ~ as.numeric(rtrend)))$r.squared))

xreg_rob <- cbind(
  troops_dt = troops_dt,
  d2009     = as.integer(df_clean$year == 2009),
  d2020     = as.integer(df_clean$year == 2020)
)

fit_rob <- Arima(y, order = c(p_sel, 0, q_sel), xreg = xreg_rob, include.mean = TRUE)

coef_rob <- coef(fit_rob)
se_rob   <- sqrt(diag(vcov(fit_rob)))
t_rob    <- coef_rob / se_rob
p_rob    <- 2 * pnorm(-abs(t_rob))

rob_table <- data.frame(
  estimate = round(coef_rob, 5),
  std_err  = round(se_rob,   5),
  t_stat   = round(t_rob,    3),
  p_value  = round(p_rob,    4)
)

cat("\n=== [Check 1] Coefficient Table — Detrended Troops ===\n")
print(rob_table)

cat(sprintf("\n  Main spec   — beta(troops)    = %.5f  (SE = %.5f)\n",
            coef(fit)["troops"],        sqrt(vcov(fit)["troops",       "troops"])))
cat(sprintf("  Robustness  — beta(troops_dt) = %.5f  (SE = %.5f)\n",
            coef(fit_rob)["troops_dt"], sqrt(vcov(fit_rob)["troops_dt", "troops_dt"])))

# ── 2. Lagged troops (x_{t-1}, x_{t-2}) ─────────────────────────────────────
fit_lag <- function(k) {
  df_k      <- df_clean[(k + 1):nrow(df_clean), ]
  troops_k  <- df_clean$troops[1:(nrow(df_clean) - k)]
  y_k       <- ts(df_k$d_debt_gdp, start = min(df_k$year), frequency = 1)
  xreg_k    <- cbind(
    troops  = troops_k,
    d2009   = as.integer(df_k$year == 2009),
    d2020   = as.integer(df_k$year == 2020)
  )
  Arima(y_k, order = c(p_sel, 0, q_sel), xreg = xreg_k, include.mean = TRUE)
}

fit_l1 <- fit_lag(1)
fit_l2 <- fit_lag(2)

coef_summary <- function(fit, var = "troops") {
  b       <- coef(fit)[var]
  v       <- vcov(fit)[var, var]
  se      <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)  # non-PD Hessian → NA
  t       <- if (is.na(se)) NA_real_ else b / se
  p       <- if (is.na(t))  NA_real_ else 2 * pnorm(-abs(t))
  data.frame(estimate = round(b, 5), std_err = round(se, 5),
             t_stat = round(t, 3), p_value = round(p, 4))
}

cat("\n=== [Check 2] Lag Sensitivity — troops_{t-k} ===\n")
cat("\n  Lag 0 (main spec):\n");  print(coef_summary(fit,    "troops"))
cat("\n  Lag 1:\n");              print(coef_summary(fit_l1, "troops"))
cat("  (NA SE = non-PD Hessian at MLE; parameter near-unidentified)\n")
cat("\n  Lag 2:\n");              print(coef_summary(fit_l2, "troops"))

# ── 3. Defense expenditure as alternative X (FRED: FDEFX) ────────────────────
cat("\n[Check 3] Downloading FDEFX from FRED...\n")
tmp <- tempfile(fileext = ".csv")
download.file("https://fred.stlouisfed.org/graph/fredgraph.csv?id=FDEFX",
              tmp, quiet = TRUE)

fdefx_raw        <- read.csv(tmp)
names(fdefx_raw) <- c("date", "fdefx")
fdefx_raw$date   <- as.Date(fdefx_raw$date)
fdefx_raw$year   <- as.integer(format(fdefx_raw$date, "%Y"))
fdefx_raw$fdefx  <- suppressWarnings(as.numeric(fdefx_raw$fdefx))

fdefx_annual <- fdefx_raw |>
  filter(!is.na(fdefx)) |>
  group_by(year) |>
  summarise(fdefx = mean(fdefx), .groups = "drop")

df_def <- inner_join(df_clean, fdefx_annual, by = "year")
y_def  <- ts(df_def$d_debt_gdp, start = min(df_def$year), frequency = 1)

xreg_def <- cbind(
  fdefx = df_def$fdefx,
  d2009 = as.integer(df_def$year == 2009),
  d2020 = as.integer(df_def$year == 2020)
)

fit_def  <- Arima(y_def, order = c(p_sel, 0, q_sel), xreg = xreg_def, include.mean = TRUE)
cf_def   <- coef(fit_def)
se_def   <- sqrt(diag(vcov(fit_def)))

cat(sprintf("  Sample: %d-%d  (n=%d)\n", min(df_def$year), max(df_def$year), nrow(df_def)))
cat("\n=== [Check 3] Coefficient Table — Defense Expenditure (bn USD SAAR) ===\n")
print(data.frame(
  estimate = round(cf_def,                    6),
  std_err  = round(se_def,                    6),
  t_stat   = round(cf_def / se_def,           3),
  p_value  = round(2 * pnorm(-abs(cf_def / se_def)), 4)
))

# ── 4. Sub-period analysis: Cold War vs. Post-9/11 ───────────────────────────
run_subperiod <- function(label, year_min, year_max) {
  df_sub <- df_clean |> filter(year >= year_min, year <= year_max)
  cat(sprintf("\n  [%s]  n = %d (%d-%d)\n", label, nrow(df_sub), year_min, year_max))

  y_sub <- ts(df_sub$d_debt_gdp, start = year_min, frequency = 1)

  crisis_cols <- lapply(c(2009, 2020), function(cy) {
    if (cy >= year_min && cy <= year_max) as.integer(df_sub$year == cy) else NULL
  })
  crisis_mat <- do.call(cbind, Filter(Negate(is.null), crisis_cols))

  xreg_sub <- if (is.null(crisis_mat)) {
    cbind(troops = df_sub$troops)
  } else {
    cbind(troops = df_sub$troops, crisis_mat)
  }

  fit_sub <- tryCatch(
    Arima(y_sub, order = c(p_sel, 0, q_sel), xreg = xreg_sub, include.mean = TRUE),
    error = function(e) { cat("    Model failed to converge\n"); NULL }
  )
  if (!is.null(fit_sub)) print(coef_summary(fit_sub, "troops"))
  invisible(fit_sub)
}

cat("\n=== [Check 4] Sub-period Analysis — beta(troops) ===\n")
cat("\n  Full sample:"); print(coef_summary(fit, "troops"))
run_subperiod("Cold War   1951-1990", 1951, 1990)
run_subperiod("Post-9/11  2001-end",  2001, max(df_clean$year))

cat("\n✓ 04_robustness.R complete.\n")
