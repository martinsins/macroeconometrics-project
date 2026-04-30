# =============================================================================
# 03_irf.R
# Military Deployment & Public Debt — Impulse Response Functions
# =============================================================================
#
# IRF interpretation: a one-time shock of +100k troops at time t increases
# Delta_debt_t by beta*100. Through AR dynamics this propagates:
#   IRF(h) = psi_h * beta * shock_size
# where psi_h are the MA(inf) coefficients of the ARMA(p,q) part.
#
# Confidence intervals via delta method (numerical Jacobian).
# =============================================================================


SHOCK_SIZE <- 100   # +100k troops

# ── 0. Load model ─────────────────────────────────────────────────────────────
saved  <- readRDS("output/armax_fit.rds")
fit    <- saved$fit
p_sel  <- saved$p
q_sel  <- saved$q

cf <- coef(fit)
V  <- vcov(fit)

# ── 1. Extract ARMA and troops parameters ────────────────────────────────────
ar_names  <- if (p_sel > 0) paste0("ar", seq_len(p_sel)) else character(0)
ma_names  <- if (q_sel > 0) paste0("ma", seq_len(q_sel)) else character(0)

phi   <- if (p_sel > 0) cf[ar_names] else numeric(0)
theta <- if (q_sel > 0) cf[ma_names] else numeric(0)
beta  <- cf["troops"]

# ── 2. MA(inf) coefficients ───────────────────────────────────────────────────
H <- 20

get_psi <- function(phi, theta, H) {
  if (length(phi) == 0 && length(theta) == 0) return(c(1, rep(0, H)))
  c(1, ARMAtoMA(ar = phi, ma = theta, lag.max = H))
}

psi <- get_psi(phi, theta, H)

# ── 3. IRF point estimates ────────────────────────────────────────────────────
irf_point  <- psi * beta * SHOCK_SIZE
cirf_point <- cumsum(irf_point)

# ── 4. Delta-method confidence intervals ─────────────────────────────────────
# Parameters that enter the IRF: AR, MA, and troops
# All must be present in vcov matrix
irf_params <- c(ar_names, ma_names, "troops")
irf_params <- irf_params[irf_params %in% colnames(V)]   # safety filter

params0 <- cf[irf_params]
V_sub   <- V[irf_params, irf_params]
n_par   <- length(params0)
eps     <- 1e-6

# Numerical Jacobian: rows = horizons 0..H, cols = parameters
jacobian <- matrix(0, nrow = H + 1, ncol = n_par)

for (k in seq_len(n_par)) {
  pup <- params0;  pup[k] <- params0[k] + eps
  pdn <- params0;  pdn[k] <- params0[k] - eps

  phi_up <- if (length(ar_names) > 0 && all(ar_names %in% names(pup))) pup[ar_names] else numeric(0)
  phi_dn <- if (length(ar_names) > 0 && all(ar_names %in% names(pdn))) pdn[ar_names] else numeric(0)
  th_up  <- if (length(ma_names) > 0 && all(ma_names %in% names(pup))) pup[ma_names] else numeric(0)
  th_dn  <- if (length(ma_names) > 0 && all(ma_names %in% names(pdn))) pdn[ma_names] else numeric(0)
  b_up   <- pup["troops"]
  b_dn   <- pdn["troops"]

  irf_up <- get_psi(phi_up, th_up, H) * b_up * SHOCK_SIZE
  irf_dn <- get_psi(phi_dn, th_dn, H) * b_dn * SHOCK_SIZE

  jacobian[, k] <- (irf_up - irf_dn) / (2 * eps)
}

irf_var  <- diag(jacobian %*% V_sub %*% t(jacobian))
irf_se   <- sqrt(pmax(irf_var, 0))

irf_lo95 <- irf_point - 1.96 * irf_se
irf_hi95 <- irf_point + 1.96 * irf_se
irf_lo68 <- irf_point - 1.00 * irf_se
irf_hi68 <- irf_point + 1.00 * irf_se

jac_cum   <- apply(jacobian, 2, cumsum)
cirf_var  <- diag(jac_cum %*% V_sub %*% t(jac_cum))
cirf_se   <- sqrt(pmax(cirf_var, 0))
cirf_lo95 <- cirf_point - 1.96 * cirf_se
cirf_hi95 <- cirf_point + 1.96 * cirf_se

# ── 5. Print table ────────────────────────────────────────────────────────────
irf_table <- data.frame(
  horizon = 0:H,
  irf     = round(irf_point, 5),
  se      = round(irf_se,    5),
  lo95    = round(irf_lo95,  5),
  hi95    = round(irf_hi95,  5)
)

cat("\n=== IRF table (shock = +100k troops) ===\n")
print(irf_table)
cat(sprintf("\nLong-run multiplier (sum IRF over %d years): %.5f\n", H, sum(irf_point)))

# ── 6. IRF plot ───────────────────────────────────────────────────────────────
irf_df <- data.frame(
  horizon = 0:H,
  irf     = irf_point,
  lo95    = irf_lo95, hi95 = irf_hi95,
  lo68    = irf_lo68, hi68 = irf_hi68
)

p_irf <- ggplot(irf_df, aes(x = horizon)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#2166ac", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo68, ymax = hi68), fill = "#2166ac", alpha = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = irf), colour = "#2166ac", linewidth = 1) +
  geom_point(aes(y = irf), colour = "#2166ac", size = 2) +
  scale_x_continuous(breaks = seq(0, H, 2)) +
  labs(
    title    = "IRF: Effect of +100k Troops on Annual Change in Debt/GDP",
    subtitle = sprintf("ARMA(%d,%d)-X  |  shaded: 68%% and 95%% confidence bands", p_sel, q_sel),
    x        = "Horizon (years)",
    y        = "Change in Delta Debt/GDP"
  ) +
  theme_minimal(base_size = 12)

# ── 7. Cumulative IRF plot ────────────────────────────────────────────────────
cirf_df <- data.frame(
  horizon = 0:H,
  cirf    = cirf_point,
  lo95    = cirf_lo95,
  hi95    = cirf_hi95
)

p_cirf <- ggplot(cirf_df, aes(x = horizon)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#d6604d", alpha = 0.20) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(y = cirf), colour = "#d6604d", linewidth = 1) +
  geom_point(aes(y = cirf), colour = "#d6604d", size = 2) +
  scale_x_continuous(breaks = seq(0, H, 2)) +
  labs(
    title    = "Cumulative IRF: Total Effect on Debt/GDP Level",
    subtitle = sprintf("ARMA(%d,%d)-X  |  95%% confidence band", p_sel, q_sel),
    x        = "Horizon (years)",
    y        = "Cumulative change in Debt/GDP"
  ) +
  theme_minimal(base_size = 12)

combo_irf <- p_irf / p_cirf
print(combo_irf)
ggsave("output/figures/03_irf.pdf", combo_irf, width = 9, height = 9)

cat("\n✓ 03_irf.R complete. Figures saved to output/figures/03_irf.pdf\n")
