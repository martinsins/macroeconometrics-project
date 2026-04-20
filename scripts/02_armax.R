# =============================================================================
# 02_armax.R
# ARMA-X estimation and Impulse Response Functions
# HEC Lausanne – Macroeconometrics project
# =============================================================================

library(forecast)
library(ggplot2)

# --- Load processed data -----------------------------------------------------
# panel <- readRDS("data/processed/panel_clean.rds")

# --- Helper: Fit ARMA-X for one country --------------------------------------
fit_armax <- function(y, xreg, max_ar = 4, max_ma = 2) {
  # Auto-select ARMA orders by AIC
  best <- auto.arima(y, xreg = xreg, max.p = max_ar, max.q = max_ma,
                     seasonal = FALSE, stepwise = FALSE, ic = "aic")
  return(best)
}

# --- Helper: Compute IRF from ARMA-X coefficients ----------------------------
# An ARMA(p,q)-X model: y_t = mu + sum phi_i y_{t-i} + sum theta_j e_{t-j}
#                             + beta * war_t + e_t
# The IRF of a unit war shock traces the dynamic response of y over h horizons.

compute_irf <- function(model, shock = 1, horizon = 20) {
  phi   <- coef(model)[grepl("^ar", names(coef(model)))]
  theta <- coef(model)[grepl("^ma", names(coef(model)))]
  p <- length(phi)
  q <- length(theta)

  irf <- numeric(horizon)
  errors <- numeric(horizon)

  for (h in 1:horizon) {
    ar_part <- sum(phi[1:min(h-1, p)] * rev(irf[max(1,h-p):max(0,h-1)]))
    ma_part <- if (h <= q) theta[h] * shock else 0
    irf[h]  <- (if (h == 1) shock else 0) + ar_part + ma_part
  }
  return(irf)
}

# --- Plot IRF -----------------------------------------------------------------
plot_irf <- function(irf_vals, country = "Country", ci = NULL) {
  df <- data.frame(horizon = seq_along(irf_vals), irf = irf_vals)
  p <- ggplot(df, aes(x = horizon, y = irf)) +
    geom_line(linewidth = 1, color = "#2c7bb6") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = paste("IRF – War shock on Debt/GDP –", country),
         x = "Years after shock", y = "Response (pp of GDP)") +
    theme_minimal(base_size = 12)
  return(p)
}

# --- Main loop over countries ------------------------------------------------
# countries <- unique(panel$country)
# results <- list()
#
# for (ctry in countries) {
#   df_c  <- filter(panel, country == ctry)
#   y     <- ts(df_c$debt_gdp, start = min(df_c$year))
#   xreg  <- as.matrix(select(df_c, war_dummy))
#   model <- fit_armax(y, xreg)
#   irf   <- compute_irf(model, shock = 1, horizon = 20)
#   results[[ctry]] <- list(model = model, irf = irf)
# }

cat("ARMA-X estimation script ready – uncomment main loop when data is loaded.\n")
