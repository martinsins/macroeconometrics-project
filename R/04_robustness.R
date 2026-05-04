# =============================================================================
# 04_robustness.R
# Military Deployment & Public Debt — Robustness Checks
# =============================================================================
#
# Loads the main ARMA(p,q)-X model estimated in 02_armax.R and re-estimates
# it under alternative specifications to assess coefficient stability.
# =============================================================================

sep  <- strrep("=", 64)
sep2 <- strrep("-", 64)

# ── 0. Load main model ────────────────────────────────────────────────────────
saved    <- readRDS("output/armax_fit.rds")
fit      <- saved$fit
p_sel    <- saved$p
q_sel    <- saved$q
df_clean <- saved$df

y <- ts(df_clean$d_debt_gdp_dt, start = 1951, frequency = 1)

cat(sprintf("\n%s\n  Loaded ARMA(%d,%d)-X  |  n = %d  (%d-%d)\n%s\n",
            sep, p_sel, q_sel, nrow(df_clean),
            min(df_clean$year), max(df_clean$year), sep))

# ── Helper: extract one-row summary for a coefficient ────────────────────────
coef_row <- function(fit, var = "troops") {
  b  <- coef(fit)[var]
  v  <- vcov(fit)[var, var]
  se <- if (is.na(v) || v < 0) NA_real_ else sqrt(v)
  t  <- if (is.na(se)) NA_real_ else b / se
  p  <- if (is.na(t))  NA_real_ else 2 * pnorm(-abs(t))
  data.frame(
    estimate = round(b,  5),
    std_err  = round(se, 5),
    t_stat   = round(t,  3),
    p_value  = round(p,  4),
    signif   = ifelse(is.na(p), "   ",
               ifelse(p < 0.01, "***",
               ifelse(p < 0.05, "** ",
               ifelse(p < 0.10, "*  ", "   "))))
  )
}

# ── Helper: full coefficient table ───────────────────────────────────────────
full_coef_table <- function(fit) {
  cf <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  t  <- cf / se
  p  <- 2 * pnorm(-abs(t))
  data.frame(
    estimate = round(cf, 5),
    std_err  = round(se, 5),
    t_stat   = round(t,  3),
    p_value  = round(p,  4),
    signif   = ifelse(p < 0.01, "***",
               ifelse(p < 0.05, "** ",
               ifelse(p < 0.10, "*  ", "   ")))
  )
}

# =============================================================================
# CHECK 1 — Lag sensitivity: troops_{t-k} for k = 0, 1, 2
# =============================================================================

fit_lag <- function(k) {
  df_k     <- df_clean[(k + 1):nrow(df_clean), ]
  troops_k <- df_clean$troops_dt[1:(nrow(df_clean) - k)]
  y_k      <- ts(df_k$d_debt_gdp_dt, start = min(df_k$year), frequency = 1)
  xreg_k   <- cbind(
    troops = troops_k,
    d2009  = as.integer(df_k$year == 2009),
    d2020  = as.integer(df_k$year == 2020)
  )
  Arima(y_k, order = c(p_sel, 0, q_sel), xreg = xreg_k, include.mean = TRUE)
}

fit_l1 <- fit_lag(1)
fit_l2 <- fit_lag(2)

lag_table <- rbind(
  cbind(spec = "Lag 0 — troops_t   (baseline)", coef_row(fit,    "troops")),
  cbind(spec = "Lag 1 — troops_{t-1}          ", coef_row(fit_l1, "troops")),
  cbind(spec = "Lag 2 — troops_{t-2}          ", coef_row(fit_l2, "troops"))
)

cat(sprintf("\n%s\n  CHECK 1 — Lag Sensitivity  |  beta(troops_{t-k})\n%s\n\n",
            sep, sep))
print(lag_table, row.names = FALSE)

# =============================================================================
# CHECK 2 — Alternative X: defense expenditure (FRED: FDEFX)
# =============================================================================

cat(sprintf("\n%s\n  CHECK 2 — Alternative X: Defense Expenditure (FDEFX)\n%s\n",
            sep, sep))
cat("  Downloading FDEFX from FRED...\n")

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

df_def   <- inner_join(df_clean, fdefx_annual, by = "year")
y_def    <- ts(df_def$d_debt_gdp_dt, start = min(df_def$year), frequency = 1)
xreg_def <- cbind(
  fdefx = df_def$fdefx,
  d2009 = as.integer(df_def$year == 2009),
  d2020 = as.integer(df_def$year == 2020)
)
fit_def <- Arima(y_def, order = c(p_sel, 0, q_sel),
                 xreg = xreg_def, include.mean = TRUE)

cat(sprintf("  Sample: %d-%d  (n=%d)\n\n", min(df_def$year), max(df_def$year), nrow(df_def)))
print(full_coef_table(fit_def))

# =============================================================================
# CHECK 3 — Sub-period analysis: Cold War vs. Post-9/11
# =============================================================================

run_subperiod <- function(label, year_min, year_max) {
  df_sub <- df_clean |> filter(year >= year_min, year <= year_max)

  y_sub <- ts(df_sub$d_debt_gdp_dt, start = year_min, frequency = 1)

  crisis_cols <- lapply(c(2009, 2020), function(cy) {
    if (cy >= year_min && cy <= year_max) as.integer(df_sub$year == cy) else NULL
  })
  crisis_mat <- do.call(cbind, Filter(Negate(is.null), crisis_cols))

  xreg_sub <- if (is.null(crisis_mat)) {
    cbind(troops = df_sub$troops_dt)
  } else {
    cbind(troops = df_sub$troops_dt, crisis_mat)
  }

  fit_sub <- tryCatch(
    Arima(y_sub, order = c(p_sel, 0, q_sel), xreg = xreg_sub, include.mean = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit_sub)) {
    return(cbind(
      spec = label,
      data.frame(n = nrow(df_sub), estimate = NA, std_err = NA,
                 t_stat = NA, p_value = NA, signif = "   ")
    ))
  }
  cbind(spec   = sprintf("%-30s (n=%d)", label, nrow(df_sub)),
        coef_row(fit_sub, "troops"))
}

sub_table <- rbind(
  cbind(spec = sprintf("%-30s (n=%d)", "Full sample 1951-end",
                       nrow(df_clean)), coef_row(fit, "troops")),
  run_subperiod("Cold War   1951-1990", 1951, 1990),
  run_subperiod("Post-9/11  2001-end",  2001, max(df_clean$year))
)

cat(sprintf("\n%s\n  CHECK 3 — Sub-period Analysis  |  beta(troops)\n%s\n\n",
            sep, sep))
print(sub_table, row.names = FALSE)

# =============================================================================
# FINAL SUMMARY — beta(X) across all specifications
# =============================================================================

summary_table <- rbind(
  cbind(spec = sprintf("Main spec   ARMA(%d,%d)-X  troops_dt lag0", p_sel, q_sel),
        coef_row(fit,    "troops")),
  cbind(spec = sprintf("Lag 1       ARMA(%d,%d)-X  troops_dt lag1", p_sel, q_sel),
        coef_row(fit_l1, "troops")),
  cbind(spec = sprintf("Lag 2       ARMA(%d,%d)-X  troops_dt lag2", p_sel, q_sel),
        coef_row(fit_l2, "troops")),
  cbind(spec = sprintf("Alt. X      ARMA(%d,%d)-X  FDEFX          ", p_sel, q_sel),
        coef_row(fit_def, "fdefx")),
  cbind(spec = "Sub-period  Cold War 1951-1990            ",
        run_subperiod("Cold War   1951-1990", 1951, 1990)[, -1]),
  cbind(spec = "Sub-period  Post-9/11 2001-end            ",
        run_subperiod("Post-9/11  2001-end", 2001, max(df_clean$year))[, -1])
)

cat(sprintf("\n%s\n  SUMMARY — beta(X) across all specifications\n%s\n\n", sep, sep))
print(summary_table, row.names = FALSE)
cat("\n  Signif. codes: *** p<0.01  ** p<0.05  * p<0.10\n")

# ── Save outputs ──────────────────────────────────────────────────────────────
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

# Table: CSV
write.csv(summary_table, "output/tables/T6_robustness.csv", row.names = FALSE)

# Coefficient plot across all specifications
rob_plot_df <- data.frame(
  spec     = summary_table$spec,
  estimate = as.numeric(summary_table$estimate),
  se       = as.numeric(summary_table$std_err)
) |>
  dplyr::mutate(
    lo95 = estimate - 1.96 * se,
    hi95 = estimate + 1.96 * se,
    spec = factor(spec, levels = rev(spec))
  )

p_rob <- ggplot(rob_plot_df, aes(x = estimate, y = spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0.2,
                 colour = "#2166ac", linewidth = 0.7) +
  geom_point(size = 3, colour = "#2166ac") +
  labs(
    title    = "Robustness: beta(X) across specifications",
    subtitle = "Point estimate + 95% CI  |  dependent variable: Delta Debt/GDP",
    x        = "Coefficient estimate",
    y        = NULL
  ) +
  theme_minimal(base_size = 11)

print(p_rob)
ggsave("output/figures/04_robustness.pdf", p_rob, width = 9, height = 5)
cat("Saved: output/tables/T6_robustness.csv\n")
cat("Saved: output/figures/04_robustness.pdf\n")

cat(sprintf("\n✓ 04_robustness.R complete.\n"))
