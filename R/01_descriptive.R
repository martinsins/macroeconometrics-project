# ============================================================
# 01_descriptive.R
# Military Deployment & Public Debt — Descriptive Analysis
# ============================================================


# ── 0. Load data ────────────────────────────────────────────
df <- read.csv("data/processed/us_ts_data.csv")

# Drop first row (d_debt_gdp is NA for 1950)
df_clean <- df |> filter(!is.na(d_debt_gdp))

# Convert to ts objects (annual, start 1951)
y <- ts(df_clean$d_debt_gdp, start = 1951, frequency = 1)
x <- ts(df_clean$troops,     start = 1951, frequency = 1)

# ── 1. Time-series plots ─────────────────────────────────────

# Conflict periods (start, end, label)
conflicts <- data.frame(
  xmin  = c(1950, 1964, 1990, 2001),
  xmax  = c(1953, 1973, 1991, 2011),
  label = c("Korea", "Vietnam", "Gulf War", "Iraq/\nAfghanistan")
)

# Helper: add shaded conflict bands to any ggplot
add_conflicts <- function(p, ymin, ymax, label_y) {
  p +
    geom_rect(data = conflicts,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              fill = "grey70", alpha = 0.25) +
    geom_text(data = conflicts,
              aes(x = (xmin + xmax) / 2, y = label_y, label = label),
              inherit.aes = FALSE,
              size = 2.6, colour = "grey40", lineheight = 0.85,
              vjust = 1, fontface = "italic")
}

# 1a. Debt-to-GDP level
p1_base <- ggplot(df, aes(x = year, y = debt_gdp)) +
  geom_line(colour = "#2166ac", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "US Federal Debt-to-GDP (level)",
       x = NULL, y = "Debt / GDP") +
  theme_minimal(base_size = 11)
p1 <- add_conflicts(p1_base, ymin = -Inf, ymax = Inf, label_y = 1.27)

# 1b. First difference of Debt/GDP  (our dependent variable y_t)
p2_base <- ggplot(df_clean, aes(x = year, y = d_debt_gdp)) +
  geom_line(colour = "#d6604d", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "Delta Debt-to-GDP (first difference - dependent variable)",
       x = NULL, y = "Delta Debt / GDP") +
  theme_minimal(base_size = 11)
p2 <- add_conflicts(p2_base, ymin = -Inf, ymax = Inf, label_y = 0.20)

# 1c. Troops deployed abroad  (our exogenous variable x_t)
p3_base <- ggplot(df_clean, aes(x = year, y = troops)) +
  geom_area(fill = "#4d9221", alpha = 0.25) +
  geom_line(colour = "#4d9221", linewidth = 0.8) +
  labs(title = "US Troops Deployed Abroad (thousands - exogenous regressor)",
       x = NULL, y = "Troops (thousands)") +
  theme_minimal(base_size = 11)
p3 <- add_conflicts(p3_base, ymin = -Inf, ymax = Inf, label_y = 1060)

combo_ts <- p1 / p2 / p3
print(combo_ts)
ggsave("output/figures/01_time_series.pdf", combo_ts,
       width = 9, height = 10)


# ── 2. Summary statistics ────────────────────────────────────
cat("\n=== Summary Statistics ===\n")
print(summary(df_clean[, c("d_debt_gdp", "troops")]))

cat("\nCorrelation (Δdebt_gdp ~ troops):",
    round(cor(df_clean$d_debt_gdp, df_clean$troops), 3), "\n")


# ── 3. ACF / PACF ────────────────────────────────────────────

acf_pacf_plot <- function(series, title, lag_max = 20) {
  acf_vals  <- acf(series,  lag.max = lag_max, plot = FALSE)
  pacf_vals <- pacf(series, lag.max = lag_max, plot = FALSE)
  ci <- qnorm(0.975) / sqrt(length(series))
  
  acf_df  <- data.frame(lag = acf_vals$lag[-1],
                        value = acf_vals$acf[-1], type = "ACF")
  pacf_df <- data.frame(lag = pacf_vals$lag,
                        value = pacf_vals$acf,    type = "PACF")
  bind_rows(acf_df, pacf_df) |>
    ggplot(aes(x = lag, y = value)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-ci, ci),
               linetype = "dashed", colour = "steelblue") +
    geom_segment(aes(xend = lag, yend = 0), linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap(~type) +
    labs(title = title, x = "Lag", y = NULL) +
    theme_minimal(base_size = 11)
}

p_acf_y <- acf_pacf_plot(y, "ACF / PACF — Δ Debt-to-GDP")
p_acf_x <- acf_pacf_plot(x, "ACF / PACF — Troops Deployed")

combo_acf <- p_acf_y / p_acf_x
print(combo_acf)
ggsave("output/figures/01_acf_pacf.pdf", combo_acf,
       width = 9, height = 8)


# ── 4. Scatter: Δ debt vs troops ─────────────────────────────
p_scatter <- ggplot(df_clean, aes(x = troops, y = d_debt_gdp,
                                  label = year)) +
  geom_point(colour = "#762a83", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "#1b7837",
              linewidth = 0.8) +
  geom_text(check_overlap = TRUE, size = 2.5, vjust = -0.5,
            colour = "grey40") +
  labs(title = "Δ Debt-to-GDP vs Troops Deployed (scatter)",
       x = "Troops (thousands)", y = "Δ Debt / GDP") +
  theme_minimal(base_size = 11)

print(p_scatter)
ggsave("output/figures/01_scatter.pdf", p_scatter, width = 7, height = 5)

cat("\n✓ 01_descriptive.R complete. Figures saved to output/figures/\n")
cat("  -> Run 01b_stationarity.R for formal unit-root tests and detrending.\n")

