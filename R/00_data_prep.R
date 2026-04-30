# =============================================================================
# 00_data_prep.R
# Military Deployment & Public Debt — ARMA-X Project
# HEC Lausanne, Spring 2026
#
# Produces: data/processed/us_ts_data.csv
#   - year         : fiscal year (integer)
#   - debt_gdp     : Gross Federal Debt as % of GDP (JST database, "debtgdp")
#   - d_debt_gdp   : first difference of debt_gdp
#   - troops       : total US troops deployed overseas (thousands)
#
# Sources:
#   - Debt/GDP : Jordà-Schularick-Taylor Macrohistory Database (AEC package)
#   - Troops   : Google Sheets compilation of DMDC/Heritage Foundation data
#                (1950–2023, annual, by country)
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tseries)
library(urca)
library(forecast)

message("Loading JST data...")

# Download JST dataset if not already present
jst_path <- "data/raw/JSTdatasetR6.xlsx"
if (!file.exists(jst_path)) {
  dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
  download.file(
    url      = "https://www.macrohistory.net/app/download/9834519560/JSTdatasetR6.xlsx",
    destfile = jst_path,
    mode     = "wb"
  )
}

debt <- read_excel(jst_path) |>
  filter(iso == "USA") |>
  select(year, debt_gdp = debtgdp) |>
  filter(!is.na(debt_gdp)) |>
  arrange(year)

message(sprintf("  JST (USA): %d annual observations (%d - %d)",
                nrow(debt), min(debt$year), max(debt$year)))

# ── Paths ────────────────────────────────────────────────────────────────────
OUTPUT_PATH <- "data/processed/us_ts_data.csv"

TROOPS_URL  <- paste0(
  "https://docs.google.com/spreadsheets/d/",
  "1xCBWFXoNDnnnBMhXALmf-w_NoyPcIiqxUqQQBNceZY8",
  "/export?format=csv&gid=0"
)

# ── 2. Load troop deployment data ─────────────────────────────────────────────
# Google Sheets CSV export uses commas in large numbers ("1,082,777") —
# strip them before numeric conversion.

message("Downloading troop deployment data from Google Sheets...")

tmp <- tempfile(fileext = ".csv")
download.file(TROOPS_URL, tmp, quiet = TRUE)
raw_csv <- read.csv(tmp, header = FALSE, stringsAsFactors = FALSE,
                    fill = TRUE, check.names = FALSE)

# Locate the header row of the year matrix
# (col1 = empty, col2 = "1950-2023", col3 = "1950", ...)
header_row <- which(apply(raw_csv, 1, function(r) {
  trimws(r[2]) == "1950-2023" && trimws(r[3]) == "1950"
}))

if (length(header_row) == 0) stop("Could not locate the year-matrix header row.")
message(sprintf("  Header row found at row %d", header_row))

# Extract years (cols 3 onward; col1 = country name, col2 = cumulative total)
years_header <- suppressWarnings(
  as.integer(trimws(raw_csv[header_row, -(1:2)]))
)
message(sprintf("  Years: %d to %d (%d values, %d NAs)",
                min(years_header, na.rm = TRUE),
                max(years_header, na.rm = TRUE),
                length(years_header),
                sum(is.na(years_header))))

# Find and extract the "TOTAL (overseas)" row
data_rows <- raw_csv[(header_row + 1):nrow(raw_csv), ]
total_row  <- which(trimws(data_rows[, 1]) == "TOTAL (overseas)")

if (length(total_row) == 0) stop("Could not find 'TOTAL (overseas)' row.")
message(sprintf("  'TOTAL (overseas)' found at data row %d", total_row[1]))

raw_values    <- trimws(data_rows[total_row[1], -(1:2)])
troops_values <- suppressWarnings(as.numeric(gsub(",", "", raw_values)))

message(sprintf("  Values: %d total, %d NAs",
                length(troops_values), sum(is.na(troops_values))))

troops <- data.frame(
  year   = years_header,
  troops = troops_values / 1000    # convert to thousands
) |>
  filter(!is.na(year), !is.na(troops))

message(sprintf("  Troops: %d annual observations (%d - %d)",
                nrow(troops), min(troops$year), max(troops$year)))

# ── 3. Merge & compute first difference ───────────────────────────────────────
message("Merging series...")

# Common sample: 1950 to the last year available in both series
start_year <- 1950
end_year   <- min(max(debt$year), max(troops$year))
message(sprintf("  Common sample: %d - %d", start_year, end_year))

df <- inner_join(debt, troops, by = "year") |>
  arrange(year) |>
  filter(year >= start_year, year <= end_year) |>
  mutate(d_debt_gdp = debt_gdp - lag(debt_gdp)) |>
  select(year, debt_gdp, d_debt_gdp, troops)

message(sprintf("  Merged: %d observations (%d - %d)",
                nrow(df), min(df$year), max(df$year)))

# ── 4. Sanity checks ──────────────────────────────────────────────────────────
stopifnot(
  "NAs in debt_gdp"             = sum(is.na(df$debt_gdp))   == 0,
  "NAs in troops"               = sum(is.na(df$troops))      == 0,
  "d_debt_gdp NA only in row 1" = sum(is.na(df$d_debt_gdp)) == 1
)

# ── 5. Save ───────────────────────────────────────────────────────────────────
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
write.csv(df, OUTPUT_PATH, row.names = FALSE)

message(sprintf("\nSaved to %s", OUTPUT_PATH))

message("\nFirst 8 rows:")
print(head(df, 8))

message("\nLast 5 rows:")
print(tail(df, 5))

message("\nSummary:")
print(summary(df[, c("debt_gdp", "d_debt_gdp", "troops")]))