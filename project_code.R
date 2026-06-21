###########################################################################
# install packages
packages <- c("dplyr", "ggplot2", "lubridate", "readr", "tidyr", "purrr",
              "scales", "maps", "viridis")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}
#############################################################################
# read data
data <- readLines("SearchResults.txt", warn = FALSE)
data <- data[grepl("^\\d{4}/\\d{2}/\\d{2}", data)]

catalog <- read.table(
  text = data,
  header = FALSE,
  stringsAsFactors = FALSE
)

names(catalog) <- c(
  "date", "time", "event_type", "geo_type",
  "mag", "mag_type", "lat", "lon", "depth",
  "quality", "evid", "nph", "ngrm"
)

catalog <- catalog %>%
  mutate(
    datetime = ymd_hms(paste(date, time)),
    mag = as.numeric(mag),
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    depth = as.numeric(depth)
  ) %>%
  filter(!is.na(mag)) %>%
  arrange(datetime)
###############################################################################
# check data and some exploratory analysis
cat("Number of events:", nrow(catalog), "\n")
cat("Time range:", as.character(min(catalog$datetime)), "to",
    as.character(max(catalog$datetime)), "\n")
cat("Magnitude range:", min(catalog$mag), "to", max(catalog$mag), "\n")

head(catalog)
summary(catalog$mag)
###############################################################################
# ============================================================
# step 1: Estimate catalogue completeness magnitude Mc
# Methods:
# 1. FMD plot
# 2. MAXC
# 3. MBS-WW
# 4. GFT-95% / GFT-90%
# 5. Mc sensitivity thresholds
# ============================================================

# Frequency-Magnitude Distribution
bin_width <- 0.1

make_fmd <- function(mags, bin_width = 0.1) {
  mags <- mags[!is.na(mags)]
  
  min_m <- floor(min(mags) / bin_width) * bin_width
  max_m <- ceiling(max(mags) / bin_width) * bin_width
  
  bins <- seq(min_m, max_m, by = bin_width)
  
  fmd <- data.frame(mag_bin = bins) %>%
    mutate(
      count = sapply(mag_bin, function(m) {
        sum(mags >= m & mags < m + bin_width)
      }),
      cumulative_count = sapply(mag_bin, function(m) {
        sum(mags >= m)
      })
    ) %>%
    filter(cumulative_count > 0)
  
  return(fmd)
}

fmd <- make_fmd(catalog$mag, bin_width)

# Incremental FMD
p1 <- ggplot(fmd, aes(x = mag_bin, y = count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(
    title = "Incremental Frequency-Magnitude Distribution",
    x = "Magnitude",
    y = "Number of events"
  ) +
  theme_minimal()

# Cumulative FMD
p2 <- ggplot(fmd, aes(x = mag_bin, y = cumulative_count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(
    title = "Cumulative Frequency-Magnitude Distribution",
    x = "Magnitude",
    y = "Cumulative number of events"
  ) +
  theme_minimal()

print(p1)
print(p2)

# MAXC method

estimate_maxc <- function(mags, bin_width = 0.1) {
  fmd <- make_fmd(mags, bin_width)
  maxc <- fmd$mag_bin[which.max(fmd$count)]
  return(maxc)
}

mc_maxc <- estimate_maxc(catalog$mag, bin_width)

cat("MAXC Mc =", mc_maxc, "\n")

#  b-value MLE with magnitude bin correction （used for GFT and MBS-WW）

estimate_b_value <- function(mags, mc, bin_width = 0.1) {
  x <- mags[mags >= mc]
  n <- length(x)
  
  if (n < 30) {
    return(data.frame(
      mc = mc,
      n = n,
      b = NA,
      delta_b = NA
    ))
  }
  
  mean_mag <- mean(x)
  
  # Utsu / Aki MLE with bin correction
  b <- log10(exp(1)) / (mean_mag - (mc - bin_width / 2))
  
  # Shi & Bolt uncertainty
  delta_b <- 2.3 * b^2 * sqrt(sum((x - mean_mag)^2) / (n * (n - 1)))
  
  return(data.frame(
    mc = mc,
    n = n,
    b = b,
    delta_b = delta_b
  ))
}

# MBS-WW method

estimate_mbs_ww <- function(mags, bin_width = 0.1, delta_m = 0.5) {
  
  mags <- mags[!is.na(mags)]
  
  candidates <- seq(
    floor(min(mags) / bin_width) * bin_width,
    floor(max(mags) / bin_width) * bin_width - delta_m,
    by = bin_width
  )
  
  b_table <- purrr::map_dfr(
    candidates,
    ~ estimate_b_value(mags, .x, bin_width)
  )
  
  b_table <- b_table %>%
    filter(!is.na(b))
  
  mbs_table <- b_table %>%
    rowwise() %>%
    mutate(
      b_avg = mean(
        b_table$b[
          b_table$mc >= mc &
            b_table$mc <= mc + delta_m
        ],
        na.rm = TRUE
      ),
      abs_diff = abs(b - b_avg),
      stable = abs_diff <= delta_b
    ) %>%
    ungroup()
  
  mc_est <- mbs_table %>%
    filter(stable == TRUE) %>%
    slice(1) %>%
    pull(mc)
  
  if (length(mc_est) == 0) mc_est <- NA
  
  return(list(
    mc = mc_est,
    table = mbs_table
  ))
}

mbs <- estimate_mbs_ww(catalog$mag, bin_width, delta_m = 0.5)

mc_mbs <- mbs$mc

cat("MBS-WW Mc =", mc_mbs, "\n")

# Plot b-value stability
ggplot(mbs$table, aes(x = mc, y = b)) +
  geom_line() +
  geom_point() +
  geom_line(aes(y = b_avg), linetype = "dashed") +
  geom_errorbar(
    aes(ymin = b - delta_b, ymax = b + delta_b),
    width = 0.03,
    alpha = 0.4
  ) +
  geom_vline(xintercept = mc_mbs, linetype = "dashed") +
  labs(
    title = "MBS-WW b-value Stability Method",
    x = "Candidate Mc",
    y = "b-value"
  ) +
  theme_minimal()

# GFT check at MBS-WW Mc

gft_at_mc <- function(mags, mc, bin_width = 0.1) {
  
  mags <- mags[!is.na(mags)]
  
  # Keep events above candidate Mc
  x <- mags[mags >= mc]
  
  if (length(x) < 50) {
    warning("Too few events above this Mc.")
    return(data.frame(
      mc = mc,
      n = length(x),
      b = NA,
      gft = NA
    ))
  }
  
  # Estimate b-value at this Mc
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  
  if (is.na(b)) {
    warning("b-value could not be estimated.")
    return(data.frame(
      mc = mc,
      n = length(x),
      b = NA,
      gft = NA
    ))
  }
  
  # Magnitude bins above Mc
  mag_bins <- seq(mc, max(mags), by = bin_width)
  
  # Observed cumulative FMD
  obs <- sapply(mag_bins, function(m) {
    sum(mags >= m)
  })
  
  # Predicted cumulative FMD under Gutenberg-Richter law
  pred <- obs[1] * 10^(-b * (mag_bins - mc))
  
  # Goodness-of-fit value
  gft <- 100 - (sum(abs(obs - pred)) / sum(obs)) * 100
  
  # Return table for this Mc
  result <- data.frame(
    mc = mc,
    n = length(x),
    b = b,
    gft = gft
  )
  
  return(result)
}

# Use MBS-WW result
mc_mbs <- 1.7

gft_check_mbs <- gft_at_mc(
  mags = catalog$mag,
  mc = mc_mbs,
  bin_width = bin_width
)

print(gft_check_mbs)

cat("GFT value at MBS-WW Mc =", gft_check_mbs$gft, "\n")

if (gft_check_mbs$gft >= 95) {
  cat("Result: Mc = 1.7 passes the GFT-95% criterion.\n")
} else if (gft_check_mbs$gft >= 90) {
  cat("Result: Mc = 1.7 passes the GFT-90% criterion, but not GFT-95%.\n")
} else {
  cat("Result: Mc = 1.7 does not pass GFT-90%; consider using a higher Mc.\n")
}

# Plot observed vs predicted cumulative FMD at Mc = 1.7

plot_gft_at_mc <- function(mags, mc, bin_width = 0.1) {
  
  mags <- mags[!is.na(mags)]
  
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  
  mag_bins <- seq(mc, max(mags), by = bin_width)
  
  obs <- sapply(mag_bins, function(m) {
    sum(mags >= m)
  })
  
  pred <- obs[1] * 10^(-b * (mag_bins - mc))
  
  plot_data <- data.frame(
    mag_bin = mag_bins,
    observed = obs,
    predicted = pred
  )
  
  ggplot(plot_data, aes(x = mag_bin)) +
    geom_point(aes(y = observed), size = 2) +
    geom_line(aes(y = observed), linewidth = 0.8) +
    geom_line(aes(y = predicted), linetype = "dashed", linewidth = 0.9) +
    scale_y_log10() +
    geom_vline(xintercept = mc, linetype = "dotted") +
    labs(
      title = paste0("GFT Check at MBS-WW Mc = ", mc),
      subtitle = paste0(
        "Observed cumulative FMD vs Gutenberg-Richter prediction; b = ",
        round(b, 3)
      ),
      x = "Magnitude",
      y = "Cumulative number of events"
    ) +
    theme_minimal()
}

plot_gft_at_mc(
  mags = catalog$mag,
  mc = 1.7,
  bin_width = bin_width
)
##############################################################################
# ============================================================
# Monte Carlo envelope for GFT / cumulative FMD
# Test whether high-magnitude deviations are random tail effects
# ============================================================

gft_monte_carlo_envelope <- function(
    mags,
    mc = 1.7,
    bin_width = 0.1,
    n_sim = 1000,
    conf = 0.95,
    seed = 123
) {
  
  set.seed(seed)
  
  mags <- mags[!is.na(mags)]
  mags_mc <- mags[mags >= mc]
  
  n <- length(mags_mc)
  
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  
  mag_bins <- seq(mc, max(mags_mc), by = bin_width)
  
  observed <- sapply(mag_bins, function(m) {
    sum(mags_mc >= m)
  })
  
  # Theoretical expected cumulative number
  predicted <- n * 10^(-b * (mag_bins - mc))
  
  # Simulate magnitudes from GR distribution:
  # P(M >= m) = 10^(-b(m - Mc))
  # If U ~ Uniform(0,1), M = Mc - log10(U)/b
  sim_cum <- matrix(NA, nrow = length(mag_bins), ncol = n_sim)
  
  for (s in 1:n_sim) {
    u <- runif(n)
    sim_mags <- mc - log10(u) / b
    
    sim_cum[, s] <- sapply(mag_bins, function(m) {
      sum(sim_mags >= m)
    })
  }
  
  lower_prob <- (1 - conf) / 2
  upper_prob <- 1 - lower_prob
  
  envelope <- data.frame(
    mag_bin = mag_bins,
    observed = observed,
    predicted = predicted,
    lower = apply(sim_cum, 1, quantile, probs = lower_prob, na.rm = TRUE),
    upper = apply(sim_cum, 1, quantile, probs = upper_prob, na.rm = TRUE)
  )
  
  # Check whether observed values are inside the simulation envelope
  envelope <- envelope %>%
    mutate(
      inside_envelope = observed >= lower & observed <= upper,
      expected_count = predicted
    )
  
  summary <- envelope %>%
    summarise(
      mc = mc,
      n = n,
      b = b,
      n_bins = n(),
      bins_inside = sum(inside_envelope),
      proportion_inside = mean(inside_envelope),
      bins_outside = sum(!inside_envelope)
    )
  
  return(list(
    envelope = envelope,
    summary = summary,
    b_info = b_info
  ))
}


mc_env_1_7 <- gft_monte_carlo_envelope(
  mags = catalog$mag,
  mc = 1.7,
  bin_width = bin_width,
  n_sim = 1000,
  conf = 0.95,
  seed = 123
)

print(mc_env_1_7$summary)
##############################################################################
# ensure if the large number of earthquakes in 2019 will affect the choice of Mc
# ============================================================
# Split catalogue into 2019 and non-2019 groups
# Purpose: test whether 2019 earthquake sequence controls Mc/FMD
# ============================================================

catalog_fmd_grouped <- catalog %>%
  filter(!is.na(mag), !is.na(datetime)) %>%
  mutate(
    year = lubridate::year(datetime),
    period_group = ifelse(year == 2019, "2019", "non-2019")
  )

group_count_summary <- catalog_fmd_grouped %>%
  count(period_group) %>%
  mutate(percentage = 100 * n / sum(n))

print(group_count_summary)

write.csv(
  group_count_summary,
  "FMD_group_event_counts_2019_vs_non2019.csv",
  row.names = FALSE
)
# ============================================================
# Incremental FMD comparison: 2019 vs non-2019
# ============================================================

make_fmd_by_group <- function(df, bin_width = 0.1) {
  
  df <- df %>%
    filter(!is.na(mag), !is.na(period_group))
  
  mag_min <- floor(min(df$mag) / bin_width) * bin_width
  mag_max <- ceiling(max(df$mag) / bin_width) * bin_width
  
  mag_bins <- seq(mag_min, mag_max, by = bin_width)
  
  fmd_grouped <- df %>%
    mutate(
      mag_bin = floor(mag / bin_width) * bin_width
    ) %>%
    count(period_group, mag_bin, name = "incremental_count") %>%
    tidyr::complete(
      period_group,
      mag_bin = mag_bins,
      fill = list(incremental_count = 0)
    ) %>%
    group_by(period_group) %>%
    arrange(mag_bin, .by_group = TRUE) %>%
    mutate(
      cumulative_count = purrr::map_dbl(
        mag_bin,
        ~ sum(incremental_count[mag_bin >= .x])
      )
    ) %>%
    ungroup()
  
  return(fmd_grouped)
}

fmd_2019_compare <- make_fmd_by_group(
  df = catalog_fmd_grouped,
  bin_width = bin_width
)

p_incremental_2019 <- ggplot(
  fmd_2019_compare,
  aes(x = mag_bin, y = incremental_count, linetype = period_group)
) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(
    title = "Incremental FMD: 2019 vs non-2019",
    subtitle = "Dotted line indicates Mc = 1.7",
    x = "Magnitude",
    y = "Incremental number of events",
    linetype = "Period"
  ) +
  theme_minimal()

print(p_incremental_2019)

ggsave(
  "Incremental_FMD_2019_vs_non2019.png",
  p_incremental_2019,
  width = 8,
  height = 5,
  dpi = 300
)
# ============================================================
# Cumulative FMD comparison: 2019 vs non-2019
# ============================================================

p_cumulative_2019 <- ggplot(
  fmd_2019_compare,
  aes(x = mag_bin, y = cumulative_count, linetype = period_group)
) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(
    title = "Cumulative FMD: 2019 vs non-2019",
    subtitle = "Dotted line indicates Mc = 1.7",
    x = "Magnitude",
    y = "Cumulative number of events",
    linetype = "Period"
  ) +
  theme_minimal()

print(p_cumulative_2019)

ggsave(
  "Cumulative_FMD_2019_vs_non2019.png",
  p_cumulative_2019,
  width = 8,
  height = 5,
  dpi = 300
)
# ============================================================
# Normalized cumulative FMD comparison
# This compares shape rather than total number of events
# ============================================================

fmd_2019_compare_norm <- fmd_2019_compare %>%
  group_by(period_group) %>%
  mutate(
    cumulative_norm = cumulative_count / max(cumulative_count, na.rm = TRUE)
  ) %>%
  ungroup()

p_cumulative_norm_2019 <- ggplot(
  fmd_2019_compare_norm,
  aes(x = mag_bin, y = cumulative_norm, linetype = period_group)
) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(
    title = "Normalized Cumulative FMD: 2019 vs non-2019",
    subtitle = "Curves are normalized by each group's total event count",
    x = "Magnitude",
    y = "Normalized cumulative proportion",
    linetype = "Period"
  ) +
  theme_minimal()

print(p_cumulative_norm_2019)

ggsave(
  "Normalized_Cumulative_FMD_2019_vs_non2019.png",
  p_cumulative_norm_2019,
  width = 8,
  height = 5,
  dpi = 300
)
# ============================================================
# b-value and GFT comparison at Mc = 1.7
# ============================================================

mc_main <- 1.7

fmd_group_quality <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    
    mags_group <- .x$mag
    
    b_info <- estimate_b_value(
      mags = mags_group,
      mc = mc_main,
      bin_width = bin_width
    )
    
    gft_info <- gft_at_mc(
      mags = mags_group,
      mc = mc_main,
      bin_width = bin_width
    )
    
    data.frame(
      n_total = length(mags_group),
      n_above_mc = sum(mags_group >= mc_main, na.rm = TRUE),
      mc = mc_main,
      b = b_info$b,
      delta_b = b_info$delta_b,
      gft = gft_info$gft
    )
    
  }) %>%
  ungroup()

print(fmd_group_quality)

write.csv(
  fmd_group_quality,
  "FMD_quality_2019_vs_non2019_at_Mc_1_7.csv",
  row.names = FALSE
)
# ============================================================
# b-value and GFT comparison at Mc = 1.7
# ============================================================

mc_main <- 1.7

fmd_group_quality <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    
    mags_group <- .x$mag
    
    b_info <- estimate_b_value(
      mags = mags_group,
      mc = mc_main,
      bin_width = bin_width
    )
    
    gft_info <- gft_at_mc(
      mags = mags_group,
      mc = mc_main,
      bin_width = bin_width
    )
    
    data.frame(
      n_total = length(mags_group),
      n_above_mc = sum(mags_group >= mc_main, na.rm = TRUE),
      mc = mc_main,
      b = b_info$b,
      delta_b = b_info$delta_b,
      gft = gft_info$gft
    )
    
  }) %>%
  ungroup()

print(fmd_group_quality)

write.csv(
  fmd_group_quality,
  "FMD_quality_2019_vs_non2019_at_Mc_1_7.csv",
  row.names = FALSE
)

# ============================================================
# Estimate MBS-WW Mc separately for 2019 and non-2019
# ============================================================

mc_mbs_by_group <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    
    mags_group <- .x$mag
    
    mbs_group <- estimate_mbs_ww(
      mags = mags_group,
      bin_width = bin_width,
      delta_m = 0.5
    )
    
    data.frame(
      n_total = length(mags_group),
      mc_mbs = mbs_group$mc
    )
    
  }) %>%
  ungroup()

print(mc_mbs_by_group)

write.csv(
  mc_mbs_by_group,
  "MBS_WW_Mc_2019_vs_non2019.csv",
  row.names = FALSE
)
#################################################################
##############################################################################
# sensitivity analysis
mc_main <- 1.7

mc_sensitivity <- c(1.7, 2.0, 2.5)

for (mc in mc_sensitivity) {
  
  temp_cat <- catalog %>%
    filter(mag >= mc)
  
  out_name <- paste0(
    "catalog_for_etas_Mc_",
    gsub("\\.", "_", mc),
    ".csv"
  )
  
  write.csv(temp_cat, out_name, row.names = FALSE)
  
  cat(
    "Saved:", out_name,
    "with", nrow(temp_cat),
    "events; threshold: M >=", mc, "\n"
  )
}

# Summary table
mc_sensitivity_summary <- data.frame(
  mc_threshold = mc_sensitivity,
  n_events = sapply(mc_sensitivity, function(mc) {
    sum(catalog$mag >= mc)
  })
)

print(mc_sensitivity_summary)

write.csv(
  mc_sensitivity_summary,
  "mc_sensitivity_event_counts.csv",
  row.names = FALSE
)

# ============================================================
# STEP 2: Descriptive analysis for main catalogue
# Main catalogue: Mc = 1.7
# ============================================================

# Create main catalogue for descriptive analysis
mc_main <- 1.7

catalog_main <- catalog %>%
  filter(mag >= mc_main) %>%
  mutate(
    date_only = as.Date(datetime),
    year = year(datetime),
    month = floor_date(datetime, "month"),
    year_month = format(datetime, "%Y-%m")
  ) %>%
  arrange(datetime)

cat("Main catalogue Mc =", mc_main, "\n")
cat("Number of events in main catalogue:", nrow(catalog_main), "\n")
cat("Time range:", as.character(min(catalog_main$datetime)), "to",
    as.character(max(catalog_main$datetime)), "\n")
cat("Magnitude range:", min(catalog_main$mag), "to", max(catalog_main$mag), "\n")
cat("Depth range:", min(catalog_main$depth), "to", max(catalog_main$depth), "\n")

# Catalogue overview summary

catalog_summary <- data.frame(
  statistic = c(
    "Mc threshold",
    "Number of events",
    "Start time",
    "End time",
    "Minimum magnitude",
    "Maximum magnitude",
    "Mean magnitude",
    "Median magnitude",
    "Minimum depth",
    "Maximum depth",
    "Mean depth",
    "Median depth",
    "Minimum latitude",
    "Maximum latitude",
    "Minimum longitude",
    "Maximum longitude"
  ),
  value = c(
    mc_main,
    nrow(catalog_main),
    as.character(min(catalog_main$datetime, na.rm = TRUE)),
    as.character(max(catalog_main$datetime, na.rm = TRUE)),
    round(min(catalog_main$mag, na.rm = TRUE), 2),
    round(max(catalog_main$mag, na.rm = TRUE), 2),
    round(mean(catalog_main$mag, na.rm = TRUE), 2),
    round(median(catalog_main$mag, na.rm = TRUE), 2),
    round(min(catalog_main$depth, na.rm = TRUE), 2),
    round(max(catalog_main$depth, na.rm = TRUE), 2),
    round(mean(catalog_main$depth, na.rm = TRUE), 2),
    round(median(catalog_main$depth, na.rm = TRUE), 2),
    round(min(catalog_main$lat, na.rm = TRUE), 3),
    round(max(catalog_main$lat, na.rm = TRUE), 3),
    round(min(catalog_main$lon, na.rm = TRUE), 3),
    round(max(catalog_main$lon, na.rm = TRUE), 3)
  )
)

print(catalog_summary)

write.csv(
  catalog_summary,
  "descriptive_catalog_summary_Mc_1_7.csv",
  row.names = FALSE
)

# Annual event counts
annual_counts <- catalog_main %>%
  group_by(year) %>%
  summarise(
    n_events = n(),
    max_mag = max(mag, na.rm = TRUE),
    mean_mag = mean(mag, na.rm = TRUE),
    median_mag = median(mag, na.rm = TRUE),
    .groups = "drop"
  )

print(annual_counts)

write.csv(
  annual_counts,
  "annual_event_counts_Mc_1_7.csv",
  row.names = FALSE
)

p_annual <- ggplot(annual_counts, aes(x = year, y = n_events)) +
  geom_col() +
  labs(
    title = "Annual Number of Earthquakes",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Year",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_annual)

ggsave(
  "annual_event_counts_Mc_1_7.png",
  p_annual,
  width = 8,
  height = 5,
  dpi = 300
)

# Monthly event counts
monthly_counts <- catalog_main %>%
  group_by(month) %>%
  summarise(
    n_events = n(),
    max_mag = max(mag, na.rm = TRUE),
    mean_mag = mean(mag, na.rm = TRUE),
    .groups = "drop"
  )

print(head(monthly_counts))

write.csv(
  monthly_counts,
  "monthly_event_counts_Mc_1_7.csv",
  row.names = FALSE
)

p_monthly <- ggplot(monthly_counts, aes(x = month, y = n_events)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  labs(
    title = "Monthly Number of Earthquakes",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Time",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_monthly)

ggsave(
  "monthly_event_counts_Mc_1_7.png",
  p_monthly,
  width = 9,
  height = 5,
  dpi = 300
)

# Cumulative event count through time
catalog_cum <- catalog_main %>%
  arrange(datetime) %>%
  mutate(cumulative_events = row_number())

p_cumulative <- ggplot(catalog_cum, aes(x = datetime, y = cumulative_events)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Cumulative Number of Earthquakes over Time",
    subtitle = "Steeper segments indicate periods of increased seismic activity",
    x = "Time",
    y = "Cumulative number of events"
  ) +
  theme_minimal()

print(p_cumulative)

ggsave(
  "cumulative_event_count_Mc_1_7.png",
  p_cumulative,
  width = 9,
  height = 5,
  dpi = 300
)

#  Magnitude distribution
p_mag_hist <- ggplot(catalog_main, aes(x = mag)) +
  geom_histogram(binwidth = 0.1, boundary = 0) +
  labs(
    title = "Magnitude Distribution",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Magnitude",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_mag_hist)

ggsave(
  "magnitude_distribution_Mc_1_7.png",
  p_mag_hist,
  width = 8,
  height = 5,
  dpi = 300
)

# FMD and Gutenberg-Richter plot for main catalogue
fmd_main <- make_fmd(catalog_main$mag, bin_width)

# Incremental FMD
p_fmd_incremental <- ggplot(fmd_main, aes(x = mag_bin, y = count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  geom_vline(xintercept = mc_main, linetype = "dashed") +
  labs(
    title = "Incremental Frequency-Magnitude Distribution",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Magnitude",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_fmd_incremental)

ggsave(
  "incremental_FMD_Mc_1_7.png",
  p_fmd_incremental,
  width = 8,
  height = 5,
  dpi = 300
)

# Cumulative FMD
p_fmd_cumulative <- ggplot(fmd_main, aes(x = mag_bin, y = cumulative_count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  geom_vline(xintercept = mc_main, linetype = "dashed") +
  labs(
    title = "Cumulative Frequency-Magnitude Distribution",
    subtitle = "Gutenberg-Richter plot; main catalogue M >= 1.7",
    x = "Magnitude",
    y = "Cumulative number of events"
  ) +
  theme_minimal()

print(p_fmd_cumulative)

ggsave(
  "cumulative_FMD_Mc_1_7.png",
  p_fmd_cumulative,
  width = 8,
  height = 5,
  dpi = 300
)

# Gutenberg-Richter fit and b-value

b_info_main <- estimate_b_value(catalog_main$mag, mc_main, bin_width)

print(b_info_main)

write.csv(
  b_info_main,
  "b_value_Mc_1_7.csv",
  row.names = FALSE
)

b_main <- b_info_main$b

gr_data <- fmd_main %>%
  filter(mag_bin >= mc_main) %>%
  mutate(
    predicted_count = first(cumulative_count) *
      10^(-b_main * (mag_bin - mc_main))
  )

p_gr_fit <- ggplot(gr_data, aes(x = mag_bin)) +
  geom_point(aes(y = cumulative_count)) +
  geom_line(aes(y = cumulative_count), linewidth = 0.7) +
  geom_line(
    aes(y = predicted_count),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_y_log10() +
  geom_vline(xintercept = mc_main, linetype = "dotted") +
  labs(
    title = "Gutenberg-Richter Fit",
    subtitle = paste0(
      "Main catalogue: Mc = 1.7; b = ",
      round(b_main, 3)
    ),
    x = "Magnitude",
    y = "Cumulative number of events"
  ) +
  theme_minimal()

print(p_gr_fit)

ggsave(
  "GR_fit_Mc_1_7.png",
  p_gr_fit,
  width = 8,
  height = 5,
  dpi = 300
)

#  Depth distribution

depth_summary <- catalog_main %>%
  summarise(
    min_depth = min(depth, na.rm = TRUE),
    q25_depth = quantile(depth, 0.25, na.rm = TRUE),
    median_depth = median(depth, na.rm = TRUE),
    mean_depth = mean(depth, na.rm = TRUE),
    q75_depth = quantile(depth, 0.75, na.rm = TRUE),
    max_depth = max(depth, na.rm = TRUE)
  )

print(depth_summary)

write.csv(
  depth_summary,
  "depth_summary_Mc_1_7.csv",
  row.names = FALSE
)

p_depth <- ggplot(catalog_main, aes(x = depth)) +
  geom_histogram(binwidth = 1) +
  labs(
    title = "Depth Distribution",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Depth (km)",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_depth)

ggsave(
  "depth_distribution_Mc_1_7.png",
  p_depth,
  width = 8,
  height = 5,
  dpi = 300
)

# Magnitude over time

p_mag_time <- ggplot(catalog_main, aes(x = datetime, y = mag)) +
  geom_point(alpha = 0.4, size = 0.8) +
  labs(
    title = "Earthquake Magnitudes over Time",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Time",
    y = "Magnitude"
  ) +
  theme_minimal()

print(p_mag_time)

ggsave(
  filename = "magnitude_over_time_Mc_1_7.png",
  plot = p_mag_time,
  width = 9,
  height = 5,
  dpi = 300
)

#  Spatial distribution map

usa_map <- map_data("state")

california_map <- usa_map %>%
  filter(region == "california")

p_spatial_map <- ggplot() +
  geom_polygon(
    data = california_map,
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray60"
  ) +
  geom_point(
    data = catalog_main,
    aes(x = lon, y = lat, size = mag),
    alpha = 0.35
  ) +
  coord_fixed(
    xlim = c(min(catalog_main$lon, na.rm = TRUE) - 0.2,
             max(catalog_main$lon, na.rm = TRUE) + 0.2),
    ylim = c(min(catalog_main$lat, na.rm = TRUE) - 0.2,
             max(catalog_main$lat, na.rm = TRUE) + 0.2)
  ) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(
    title = "Spatial Distribution of Earthquakes",
    subtitle = "Point size represents magnitude; main catalogue M >= 1.7",
    x = "Longitude",
    y = "Latitude",
    size = "Magnitude"
  ) +
  theme_minimal()

print(p_spatial_map)

ggsave(
  "spatial_distribution_map_Mc_1_7.png",
  p_spatial_map,
  width = 7,
  height = 6,
  dpi = 300
)
#  Spatial distribution colored by magnitude

p_spatial_mag <- ggplot(catalog_main, aes(x = lon, y = lat)) +
  geom_point(aes(color = mag), alpha = 0.5, size = 1) +
  coord_fixed() +
  scale_color_viridis_c() +
  labs(
    title = "Spatial Distribution Colored by Magnitude",
    subtitle = "Main catalogue: M >= 1.7",
    x = "Longitude",
    y = "Latitude",
    color = "Magnitude"
  ) +
  theme_minimal()

print(p_spatial_mag)

ggsave(
  "spatial_distribution_colored_by_magnitude_Mc_1_7.png",
  p_spatial_mag,
  width = 7,
  height = 6,
  dpi = 300
)
#  Top 10 largest events

top_events <- catalog_main %>%
  arrange(desc(mag)) %>%
  select(datetime, mag, lat, lon, depth, quality, evid, nph, ngrm) %>%
  slice(1:10)

print(top_events)

write.csv(
  top_events,
  "top_10_largest_events_Mc_1_7.csv",
  row.names = FALSE
)

#  Event burst around the largest earthquake

largest_event <- catalog_main %>%
  arrange(desc(mag)) %>%
  slice(1)

print(largest_event)

mainshock_time <- largest_event$datetime[1]
window_days <- 30

catalog_around_mainshock <- catalog_main %>%
  filter(
    datetime >= mainshock_time - days(window_days),
    datetime <= mainshock_time + days(window_days)
  ) %>%
  mutate(
    days_from_mainshock = as.numeric(
      difftime(datetime, mainshock_time, units = "days")
    )
  )

p_around_mag <- ggplot(catalog_around_mainshock, aes(x = days_from_mainshock, y = mag)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Seismicity around the Largest Event",
    subtitle = paste0(
      "Largest event: M = ",
      largest_event$mag,
      "; window = ±",
      window_days,
      " days"
    ),
    x = "Days from largest event",
    y = "Magnitude"
  ) +
  theme_minimal()

print(p_around_mag)

ggsave(
  "seismicity_around_largest_event_Mc_1_7.png",
  p_around_mag,
  width = 8,
  height = 5,
  dpi = 300
)

daily_around_mainshock <- catalog_around_mainshock %>%
  mutate(relative_day = floor(days_from_mainshock)) %>%
  group_by(relative_day) %>%
  summarise(
    n_events = n(),
    max_mag = max(mag, na.rm = TRUE),
    .groups = "drop"
  )

p_around_daily <- ggplot(daily_around_mainshock, aes(x = relative_day, y = n_events)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Daily Event Counts around the Largest Event",
    subtitle = paste0(
      "Largest event: M = ",
      largest_event$mag
    ),
    x = "Days from largest event",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_around_daily)

ggsave(
  "daily_counts_around_largest_event_Mc_1_7.png",
  p_around_daily,
  width = 8,
  height = 5,
  dpi = 300
)

#  Inter-event time distribution （support for ETAS）

inter_event <- catalog_main %>%
  arrange(datetime) %>%
  mutate(
    inter_event_time_hours = as.numeric(
      difftime(datetime, lag(datetime), units = "hours")
    ),
    inter_event_time_days = inter_event_time_hours / 24
  ) %>%
  filter(!is.na(inter_event_time_hours))

inter_event_summary <- inter_event %>%
  summarise(
    min_hours = min(inter_event_time_hours, na.rm = TRUE),
    q25_hours = quantile(inter_event_time_hours, 0.25, na.rm = TRUE),
    median_hours = median(inter_event_time_hours, na.rm = TRUE),
    mean_hours = mean(inter_event_time_hours, na.rm = TRUE),
    q75_hours = quantile(inter_event_time_hours, 0.75, na.rm = TRUE),
    max_hours = max(inter_event_time_hours, na.rm = TRUE)
  )

print(inter_event_summary)

write.csv(
  inter_event_summary,
  "inter_event_time_summary_Mc_1_7.csv",
  row.names = FALSE
)

write.csv(
  inter_event,
  "inter_event_times_Mc_1_7.csv",
  row.names = FALSE
)

p_inter_event <- ggplot(inter_event, aes(x = inter_event_time_hours)) +
  geom_histogram(bins = 80) +
  scale_x_log10() +
  labs(
    title = "Inter-event Time Distribution",
    subtitle = "Short inter-event times indicate temporal clustering",
    x = "Inter-event time (hours, log scale)",
    y = "Number of event pairs"
  ) +
  theme_minimal()

print(p_inter_event)

ggsave(
  "inter_event_time_distribution_Mc_1_7.png",
  p_inter_event,
  width = 8,
  height = 5,
  dpi = 300
)

large_mag_threshold <- 5.0

large_events <- catalog_main %>%
  filter(mag >= large_mag_threshold)
#  Write short text summary

summary_text <- paste0(
  "Descriptive analysis summary\n",
  "============================\n",
  "Main catalogue threshold: M >= ", mc_main, "\n",
  "Number of events: ", nrow(catalog_main), "\n",
  "Time range: ", min(catalog_main$datetime), " to ", max(catalog_main$datetime), "\n",
  "Magnitude range: ", min(catalog_main$mag), " to ", max(catalog_main$mag), "\n",
  "Depth range: ", min(catalog_main$depth), " to ", max(catalog_main$depth), " km\n",
  "Estimated b-value at Mc = ", mc_main, ": ", round(b_main, 3), "\n",
  "Largest event magnitude: ", largest_event$mag, "\n",
  "Number of events with M >= ", large_mag_threshold, ": ", nrow(large_events), "\n"
)

cat(summary_text)

writeLines(
  summary_text,
  "descriptive_analysis_text_summary_Mc_1_7.txt"
)

# ============================================================
# STEP 3: Declustering and construct the initial background rate for ETAS
# Main catalogue: Mc = 1.7
# Methods:
# A. Window-based declustering
# B. Nearest-neighbour declustering
# C. Stochastic ETAS declustering
# ============================================================

# ============================================================
# A.Gardner-Knopoff type window-based declustering
# Reference logic: OpenQuake HMTK GardnerKnopoffType1
# Main catalogue: Mc = 1.7
# ============================================================

# Prepare main catalogue
mc_main <- 1.7

catalog_main <- catalog %>%
  filter(mag >= mc_main) %>%
  arrange(datetime) %>%
  mutate(
    original_order = row_number(),
    event_index = row_number()
  )

cat("Catalogue for Gardner-Knopoff declustering:", nrow(catalog_main), "events\n")


# Haversine distance function (calculating epicentral distance in km)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  
  earth_radius <- 6371.227
  
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  a <- sin(dlat / 2)^2 +
    cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  
  earth_radius * c
}


# Gardner-Knopoff time-distance windows
# sw_space = 10^(0.1238 * M + 0.983)
# sw_time  = 10^(0.032 * M + 2.7389) for M >= 6.5
# sw_time  = 10^(0.5409 * M - 0.547) for M < 6.5
# sw_time is in days and sw_space is in km.

gk_space_window_km <- function(mag) {
  10^(0.1238 * mag + 0.983)
}

gk_time_window_days <- function(mag) {
  ifelse(
    mag < 6.5,
    10^(0.5409 * mag - 0.547),
    10^(0.032 * mag + 2.7389)
  )
}


# Gardner-Knopoff declustering function

gardner_knopoff_declustering <- function(
    df,
    fs_time_prop = 0.5,
    time_cutoff_days = NULL
) {
  
  # df must contain datetime, mag, lat, lon, evid
  df_original <- df %>%
    arrange(datetime) %>%
    mutate(original_order = row_number())
  
  n <- nrow(df_original)
  
  # Calculate time and distance windows for each event
  df_original <- df_original %>%
    mutate(
      gk_space_window_km = gk_space_window_km(mag),
      gk_time_window_days = gk_time_window_days(mag)
    )
  
  # Optional cutoff: equivalent in spirit to OpenQuake time_cutoff
  if (!is.null(time_cutoff_days)) {
    df_original <- df_original %>%
      mutate(
        gk_time_window_days = pmin(gk_time_window_days, time_cutoff_days)
      )
  }
  
  # OpenQuake sorts events by decreasing magnitude before clustering
  df_work <- df_original %>%
    arrange(desc(mag)) %>%
    mutate(sorted_index = row_number())
  
  # Initialize cluster vectors
  df_work$gk_cluster_id <- 0L
  df_work$gk_flag <- 0L
  df_work$gk_parent_evid <- NA
  df_work$gk_parent_mag <- NA_real_
  
  cluster_id <- 0L
  
  for (i in seq_len(n - 1)) {
    
    if (i %% 1000 == 0) {
      cat("Processing event", i, "of", n, "\n")
    }
    
    # Only events not yet assigned to a cluster can start a new cluster
    if (df_work$gk_cluster_id[i] != 0L) next
    
    main_time <- df_work$datetime[i]
    main_lon <- df_work$lon[i]
    main_lat <- df_work$lat[i]
    main_evid <- df_work$evid[i]
    main_mag <- df_work$mag[i]
    
    sw_time <- df_work$gk_time_window_days[i]
    sw_space <- df_work$gk_space_window_km[i]
    
    # Time difference in days: event time - main event time
    dt_days <- as.numeric(
      difftime(df_work$datetime, main_time, units = "days")
    )
    
    # OpenQuake logic:
    # dt >= -sw_time * fs_time_prop and dt <= sw_time
    # plus only unclustered events
    inside_time <- df_work$gk_cluster_id == 0L &
      dt_days >= (-sw_time * fs_time_prop) &
      dt_days <= sw_time
    
    if (sum(inside_time) <= 1) next
    
    candidate_idx <- which(inside_time)
    
    distances <- haversine_km(
      lon1 = df_work$lon[candidate_idx],
      lat1 = df_work$lat[candidate_idx],
      lon2 = main_lon,
      lat2 = main_lat
    )
    
    inside_space_idx <- candidate_idx[distances <= sw_space]
    
    # Exclude the main event itself when checking if cluster exists
    cluster_members_without_main <- setdiff(inside_space_idx, i)
    
    if (length(cluster_members_without_main) > 0) {
      
      cluster_id <- cluster_id + 1L
      
      # Assign cluster id to all selected events, including main event
      df_work$gk_cluster_id[inside_space_idx] <- cluster_id
      
      # Default clustered members are aftershocks
      df_work$gk_flag[inside_space_idx] <- 1L
      
      # Events before the main event are foreshocks
      foreshock_idx <- inside_space_idx[dt_days[inside_space_idx] < 0]
      df_work$gk_flag[foreshock_idx] <- -1L
      
      # Main event itself is flag 0
      df_work$gk_flag[i] <- 0L
      
      # Parent information for foreshocks/aftershocks
      non_main_idx <- setdiff(inside_space_idx, i)
      df_work$gk_parent_evid[non_main_idx] <- main_evid
      df_work$gk_parent_mag[non_main_idx] <- main_mag
    }
  }
  
  # Re-sort to original chronological order
  result <- df_work %>%
    arrange(original_order) %>%
    mutate(
      gk_label = case_when(
        gk_cluster_id == 0 ~ "background",
        gk_cluster_id != 0 & gk_flag == 0 ~ "mainshock",
        gk_flag == 1 ~ "aftershock",
        gk_flag == -1 ~ "foreshock",
        TRUE ~ "unknown"
      ),
      # For ETAS initialization:
      # background/mainshock-like events get probability 1,
      # foreshock/aftershock clustered events get probability 0.
      gk_background_prob = ifelse(gk_flag == 0, 1, 0)
    )
  
  return(result)
}
# Run Gardner-Knopoff type declustering

declust_gk <- gardner_knopoff_declustering(
  df = catalog_main,
  fs_time_prop = 0.5,
  time_cutoff_days = NULL
)

# Summary
gk_summary <- declust_gk %>%
  count(gk_label) %>%
  mutate(percentage = 100 * n / sum(n))

print(gk_summary)

write.csv(
  declust_gk,
  "declustering_gardner_knopoff_Mc_1_7.csv",
  row.names = FALSE
)

write.csv(
  gk_summary,
  "declustering_gardner_knopoff_summary_Mc_1_7.csv",
  row.names = FALSE
)
# Plot Gardner-Knopoff declustering result

p_gk_time <- ggplot(declust_gk, aes(x = datetime, y = mag)) +
  geom_point(aes(color = gk_label), alpha = 0.5, size = 0.8) +
  labs(
    title = "Gardner-Knopoff Type Window-based Declustering",
    subtitle = "Labels: background, mainshock, foreshock, aftershock",
    x = "Time",
    y = "Magnitude",
    color = "Event label"
  ) +
  theme_minimal()

print(p_gk_time)

ggsave(
  "declustering_gardner_knopoff_time_Mc_1_7.png",
  p_gk_time,
  width = 9,
  height = 5,
  dpi = 300
)


p_gk_map <- ggplot(declust_gk, aes(x = lon, y = lat)) +
  geom_point(aes(color = gk_label), alpha = 0.5, size = 0.8) +
  coord_fixed() +
  labs(
    title = "Spatial Distribution after Gardner-Knopoff Declustering",
    subtitle = "Window-based declustering result",
    x = "Longitude",
    y = "Latitude",
    color = "Event label"
  ) +
  theme_minimal()

print(p_gk_map)

ggsave(
  "declustering_gardner_knopoff_map_Mc_1_7.png",
  p_gk_map,
  width = 7,
  height = 6,
  dpi = 300
)

# Initial background rate from Gardner-Knopoff declustering
# Monthly background count

bg_rate_gk_monthly <- declust_gk %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarise(
    total_events = n(),
    background_events = sum(gk_background_prob, na.rm = TRUE),
    background_rate = background_events,
    .groups = "drop"
  ) %>%
  mutate(method = "Gardner-Knopoff window-based")

print(head(bg_rate_gk_monthly))

write.csv(
  bg_rate_gk_monthly,
  "initial_background_rate_gardner_knopoff_monthly_Mc_1_7.csv",
  row.names = FALSE
)

p_bg_gk <- ggplot(bg_rate_gk_monthly, aes(x = month, y = background_rate)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Initial Background Rate from Gardner-Knopoff Declustering",
    subtitle = "Monthly background/mainshock-like event count",
    x = "Time",
    y = "Background events per month"
  ) +
  theme_minimal()

print(p_bg_gk)

ggsave(
  "initial_background_rate_gardner_knopoff_monthly_Mc_1_7.png",
  p_bg_gk,
  width = 9,
  height = 5,
  dpi = 300
)
# Export ETAS initialization file

etas_init_gk <- declust_gk %>%
  select(
    evid,
    datetime,
    lat,
    lon,
    depth,
    mag,
    gk_cluster_id,
    gk_flag,
    gk_label,
    gk_background_prob,
    gk_parent_evid,
    gk_parent_mag,
    gk_space_window_km,
    gk_time_window_days
  )

write.csv(
  etas_init_gk,
  "etas_initialization_gardner_knopoff_Mc_1_7.csv",
  row.names = FALSE
)

cat("Gardner-Knopoff declustering and ETAS initialization files saved.\n")

# ============================================================
# Nearest-neighbor declustering
# Based on Zaliapin & Ben-Zion nearest-neighbor proximity
# Main catalogue: Mc = 1.7
# ============================================================

mc_main <- 1.7

catalog_main <- catalog %>%
  filter(mag >= mc_main) %>%
  arrange(datetime) %>%
  mutate(
    event_index = row_number()
  )

cat("Catalogue for nearest-neighbor declustering:", nrow(catalog_main), "events\n")
cat("Magnitude threshold: M >=", mc_main, "\n")

# ============================================================
# Helper functions
# ============================================================

# Convert datetime to fractional years
decimal_year <- function(datetime) {
  y <- lubridate::year(datetime)
  start_year <- lubridate::ymd_hms(paste0(y, "-01-01 00:00:00"))
  start_next <- lubridate::ymd_hms(paste0(y + 1, "-01-01 00:00:00"))
  
  y + as.numeric(difftime(datetime, start_year, units = "secs")) /
    as.numeric(difftime(start_next, start_year, units = "secs"))
}

# Haversine distance in km
haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371.227
  
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  a <- sin(dlat / 2)^2 +
    cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  
  R * c
}

# Estimate b-value if not already available
estimate_b_value <- function(mags, mc, bin_width = 0.1) {
  x <- mags[mags >= mc]
  n <- length(x)
  
  if (n < 30) {
    return(data.frame(
      mc = mc,
      n = n,
      b = NA,
      delta_b = NA
    ))
  }
  
  mean_mag <- mean(x)
  
  b <- log10(exp(1)) / (mean_mag - (mc - bin_width / 2))
  
  delta_b <- 2.3 * b^2 *
    sqrt(sum((x - mean_mag)^2) / (n * (n - 1)))
  
  data.frame(
    mc = mc,
    n = n,
    b = b,
    delta_b = delta_b
  )
}
# ============================================================
# Compute nearest-neighbor parent and eta
# ============================================================

compute_nearest_neighbor_eta <- function(
    df,
    mc = 1.7,
    b_value = 0.833,
    d_f = 1.6,
    use_depth = FALSE,
    max_lookback_years = 10,
    max_previous_events = 5000
) {
  
  df <- df %>%
    arrange(datetime) %>%
    mutate(
      event_index = row_number(),
      decimal_year = decimal_year(datetime),
      nn_parent_index = NA_integer_,
      nn_parent_evid = NA,
      nn_time_years = NA_real_,
      nn_distance_km = NA_real_,
      nn_parent_mag = NA_real_,
      nn_log_eta = NA_real_
    )
  
  n <- nrow(df)
  
  for (j in 2:n) {
    
    if (j %% 1000 == 0) {
      cat("Processing event", j, "of", n, "\n")
    }
    
    # Previous events only
    previous <- 1:(j - 1)
    
    # Time difference in years
    dt_years_all <- df$decimal_year[j] - df$decimal_year[previous]
    
    # Limit lookback for computational feasibility
    previous <- previous[dt_years_all > 0 & dt_years_all <= max_lookback_years]
    
    if (length(previous) == 0) next
    
    # Limit number of previous events if needed
    if (length(previous) > max_previous_events) {
      previous <- tail(previous, max_previous_events)
    }
    
    dt_years <- df$decimal_year[j] - df$decimal_year[previous]
    dt_years[dt_years <= 0] <- 1 / (365.25 * 24 * 3600)
    
    epi_dist_km <- haversine_km(
      lon1 = df$lon[previous],
      lat1 = df$lat[previous],
      lon2 = df$lon[j],
      lat2 = df$lat[j]
    )
    
    epi_dist_km[epi_dist_km <= 0] <- 0.001
    
    if (use_depth) {
      dz <- df$depth[j] - df$depth[previous]
      r_km <- sqrt(epi_dist_km^2 + dz^2)
      r_km[r_km <= 0] <- 0.001
    } else {
      r_km <- epi_dist_km
    }
    
    parent_mag <- df$mag[previous]
    
    # Nearest-neighbor proximity:
    # eta_ij = T_ij * R_ij^df * 10^[-b(M_i - Mc)]
    # We use log10 eta for numerical stability.
    log_eta <- log10(dt_years) +
      d_f * log10(r_km) -
      b_value * (parent_mag - mc)
    
    min_id <- which.min(log_eta)
    parent_idx <- previous[min_id]
    
    df$nn_parent_index[j] <- parent_idx
    df$nn_parent_evid[j] <- df$evid[parent_idx]
    df$nn_time_years[j] <- dt_years[min_id]
    df$nn_distance_km[j] <- r_km[min_id]
    df$nn_parent_mag[j] <- parent_mag[min_id]
    df$nn_log_eta[j] <- log_eta[min_id]
  }
  
  return(df)
}
# ============================================================
# Classify events using log_eta distribution
# ============================================================

classify_nn_events <- function(
    df,
    stochastic = TRUE,
    alpha0 = 0,
    seed = 123
) {
  
  eta_values <- df$nn_log_eta[!is.na(df$nn_log_eta)]
  
  if (length(eta_values) < 50) {
    stop("Too few eta values for classification.")
  }
  
  set.seed(seed)
  
  km <- kmeans(eta_values, centers = 2, nstart = 100)
  
  centers <- as.numeric(km$centers)
  
  clustered_center <- min(centers)
  background_center <- max(centers)
  eta_threshold <- mean(centers)
  
  df$nn_eta_threshold <- eta_threshold
  
  # Deterministic classification
  df <- df %>%
    mutate(
      nn_is_background_det = ifelse(
        is.na(nn_log_eta),
        TRUE,
        nn_log_eta >= eta_threshold
      ),
      nn_is_clustered_det = !nn_is_background_det
    )
  
  # Stochastic thinning proxy:
  # Repository documentation says alternative realizations can be generated using
  # p = 10^(ad0 + alpha0), where p is used to retain background events.
  # Here alpha_proxy is centered log_eta relative to the threshold.
  df <- df %>%
    mutate(
      nn_alpha_proxy = nn_log_eta - eta_threshold,
      nn_background_prob = ifelse(
        is.na(nn_alpha_proxy),
        1,
        pmin(1, pmax(0, 10^(nn_alpha_proxy + alpha0)))
      )
    )
  
  if (stochastic) {
    set.seed(seed)
    df <- df %>%
      mutate(
        nn_is_background = runif(n()) < nn_background_prob,
        nn_is_clustered = !nn_is_background
      )
  } else {
    df <- df %>%
      mutate(
        nn_is_background = nn_is_background_det,
        nn_is_clustered = nn_is_clustered_det,
        nn_background_prob = ifelse(nn_is_background_det, 1, 0)
      )
  }
  
  attr(df, "eta_threshold") <- eta_threshold
  attr(df, "eta_centers") <- centers
  
  return(df)
}
# ============================================================
# Build clusters and event labels
# ============================================================

assign_nn_clusters <- function(df) {
  
  df <- df %>%
    arrange(datetime) %>%
    mutate(
      nn_cluster_id = 0L,
      nn_mainshock_evid = NA,
      nn_label = "background"
    )
  
  n <- nrow(df)
  
  # Build child list for clustered events
  children <- vector("list", n)
  
  for (i in seq_len(n)) {
    children[[i]] <- integer(0)
  }
  
  for (j in seq_len(n)) {
    p <- df$nn_parent_index[j]
    
    if (!is.na(p) && df$nn_is_clustered[j]) {
      children[[p]] <- c(children[[p]], j)
    }
  }
  
  visited <- rep(FALSE, n)
  cluster_id <- 0L
  
  for (i in seq_len(n)) {
    
    if (visited[i]) next
    
    # Find connected component by parent-child links
    stack <- i
    component <- integer(0)
    
    while (length(stack) > 0) {
      v <- stack[1]
      stack <- stack[-1]
      
      if (visited[v]) next
      
      visited[v] <- TRUE
      component <- c(component, v)
      
      # Add children
      stack <- c(stack, children[[v]])
      
      # Add parent if event is clustered
      p <- df$nn_parent_index[v]
      if (!is.na(p) && df$nn_is_clustered[v]) {
        stack <- c(stack, p)
      }
    }
    
    if (length(component) == 1) {
      idx <- component[1]
      df$nn_label[idx] <- ifelse(
        df$nn_is_background[idx],
        "background",
        "clustered"
      )
      next
    }
    
    cluster_id <- cluster_id + 1L
    
    # Mainshock = largest magnitude event in component
    main_idx <- component[which.max(df$mag[component])]
    main_time <- df$datetime[main_idx]
    main_evid <- df$evid[main_idx]
    
    df$nn_cluster_id[component] <- cluster_id
    df$nn_mainshock_evid[component] <- main_evid
    
    for (idx in component) {
      if (idx == main_idx) {
        df$nn_label[idx] <- "mainshock"
      } else if (df$datetime[idx] < main_time) {
        df$nn_label[idx] <- "foreshock"
      } else {
        df$nn_label[idx] <- "aftershock"
      }
    }
  }
  
  # For ETAS initialization:
  # background + mainshock are background-like;
  # foreshock + aftershock are clustered/triggered.
  df <- df %>%
    mutate(
      nn_background_like = nn_label %in% c("background", "mainshock"),
      nn_init_background_prob = ifelse(nn_background_like, 1, 0)
    )
  
  return(df)
}
# ============================================================
# Run nearest-neighbor declustering
# ============================================================

# Use b-value from your GFT check if available
if (exists("gft_check_mbs")) {
  b_for_nn <- gft_check_mbs$b
} else {
  b_for_nn <- estimate_b_value(catalog_main$mag, mc_main)$b
}

cat("b-value used for NN declustering:", b_for_nn, "\n")

# Parameter df:
# df = 1.6 is commonly used for epicentral distributions.
# If you want to follow the repository note for epicenters, keep use_depth = FALSE.
df_for_nn <- 1.6

nn_raw <- compute_nearest_neighbor_eta(
  df = catalog_main,
  mc = mc_main,
  b_value = b_for_nn,
  d_f = df_for_nn,
  use_depth = FALSE,
  max_lookback_years = 10,
  max_previous_events = 5000
)

nn_classified <- classify_nn_events(
  df = nn_raw,
  stochastic = FALSE,
  alpha0 = 0,
  seed = 123
)

declust_nn <- assign_nn_clusters(nn_classified)

nn_summary <- declust_nn %>%
  count(nn_label) %>%
  mutate(percentage = 100 * n / sum(n))

print(nn_summary)

cat("Nearest-neighbor eta threshold:",
    attr(nn_classified, "eta_threshold"), "\n")

write.csv(
  declust_nn,
  "declustering_nearest_neighbor_ZB_style_Mc_1_7.csv",
  row.names = FALSE
)

write.csv(
  nn_summary,
  "declustering_nearest_neighbor_ZB_style_summary_Mc_1_7.csv",
  row.names = FALSE
)
# ============================================================
# Diagnostic plot: log_eta distribution
# ============================================================

eta_threshold <- attr(nn_classified, "eta_threshold")

p_eta <- ggplot(declust_nn, aes(x = nn_log_eta)) +
  geom_histogram(bins = 80) +
  geom_vline(xintercept = eta_threshold, linetype = "dashed") +
  labs(
    title = "Nearest-neighbor Declustering Diagnostic",
    subtitle = "Distribution of log10(eta); dashed line is classification threshold",
    x = "log10(eta)",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_eta)

ggsave(
  "nearest_neighbor_log_eta_distribution_Mc_1_7.png",
  p_eta,
  width = 8,
  height = 5,
  dpi = 300
)
# ============================================================
# Time-magnitude plot
# ============================================================

p_nn_time <- ggplot(declust_nn, aes(x = datetime, y = mag)) +
  geom_point(aes(color = nn_label), alpha = 0.5, size = 0.8) +
  labs(
    title = "Nearest-neighbor Declustering",
    subtitle = "Zaliapin-Ben-Zion style nearest-neighbor proximity",
    x = "Time",
    y = "Magnitude",
    color = "Event label"
  ) +
  theme_minimal()

print(p_nn_time)

ggsave(
  "nearest_neighbor_declustering_time_Mc_1_7.png",
  p_nn_time,
  width = 9,
  height = 5,
  dpi = 300
)
# ============================================================
# Spatial plot
# ============================================================

p_nn_map <- ggplot(declust_nn, aes(x = lon, y = lat)) +
  geom_point(aes(color = nn_label), alpha = 0.5, size = 0.8) +
  coord_fixed() +
  labs(
    title = "Spatial Distribution after Nearest-neighbor Declustering",
    subtitle = "Zaliapin-Ben-Zion style nearest-neighbor proximity",
    x = "Longitude",
    y = "Latitude",
    color = "Event label"
  ) +
  theme_minimal()

print(p_nn_map)

ggsave(
  "nearest_neighbor_declustering_map_Mc_1_7.png",
  p_nn_map,
  width = 7,
  height = 6,
  dpi = 300
)

# ============================================================
# Export ETAS initialization file
# ============================================================

etas_init_nn <- declust_nn %>%
  select(
    evid,
    datetime,
    lat,
    lon,
    depth,
    mag,
    nn_parent_evid,
    nn_parent_index,
    nn_time_years,
    nn_distance_km,
    nn_parent_mag,
    nn_log_eta,
    nn_cluster_id,
    nn_mainshock_evid,
    nn_label,
    nn_background_prob,
    nn_init_background_prob
  )

write.csv(
  etas_init_nn,
  "etas_initialization_nearest_neighbor_Mc_1_7.csv",
  row.names = FALSE
)

cat("Nearest-neighbor declustering and ETAS initialization files saved.\n")

if (exists("declust_nn")) {
  
  bg_rate_nn_monthly <- declust_nn %>%
    mutate(month = floor_date(datetime, "month")) %>%
    group_by(month) %>%
    summarise(
      total_events = n(),
      background_events = sum(nn_init_background_prob, na.rm = TRUE),
      background_rate = background_events,
      .groups = "drop"
    ) %>%
    mutate(method = "Nearest-neighbour")
  
} else {
  
  message("declust_nn does not exist. Run nearest-neighbour declustering first.")
  
}

# ============================================================
# Reasenberg-style declustering
# R adaptation based on Reasenberg (1985) / USGS CLUSTER2000 logic
# Main catalogue: Mc = 1.7
# ============================================================

mc_main <- 1.7

catalog_main <- catalog %>%
  filter(!is.na(datetime), !is.na(mag), !is.na(lat), !is.na(lon)) %>%
  filter(mag >= mc_main) %>%
  arrange(datetime) %>%
  mutate(
    event_index = row_number(),
    original_order = row_number()
  )

cat("Catalogue for Reasenberg-style declustering:",
    nrow(catalog_main), "events\n")
# ============================================================
# Haversine distance in km
# ============================================================

haversine_km <- function(lon1, lat1, lon2, lat2) {
  
  R <- 6371.227
  
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  a <- sin(dlat / 2)^2 +
    cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  
  R * c
}
# ============================================================
# Reasenberg-style interaction radius
# Approximate source/interaction radius in km
# ============================================================

reasenberg_radius_km <- function(mag, rfact = 10) {
  
  # Approximate rupture/source radius in km.
  # The multiplier rfact expands this into an interaction zone.
  source_radius <- 10^(0.4 * mag - 1.2)
  
  rfact * source_radius
}
# ============================================================
# Reasenberg-style look-ahead time window
# Bounded between tau_min and tau_max
# ============================================================

reasenberg_tau_days <- function(
    mag,
    tau_min = 1,
    tau_max = 10,
    xmeff = 4.0,
    xk = 0.5
) {
  
  # Effective magnitude increases with event magnitude.
  # This mimics the idea that larger events can keep a cluster active longer.
  m_eff <- pmax(mag, xmeff - xk)
  
  # Smoothly scale tau between tau_min and tau_max.
  # Larger events receive longer interaction times.
  scale_val <- pmin(1, pmax(0, (m_eff - mc_main) / (xmeff - mc_main)))
  
  tau <- tau_min + scale_val * (tau_max - tau_min)
  
  pmin(tau_max, pmax(tau_min, tau))
}
# ============================================================
# Reasenberg-style declustering function
# ============================================================

reasenberg_declustering <- function(
    df,
    tau_min = 1,
    tau_max = 10,
    p1 = 0.95,
    xmeff = 4.0,
    xk = 0.5,
    rfact = 10
) {
  
  df <- df %>%
    arrange(datetime) %>%
    mutate(
      event_index = row_number(),
      reasenberg_cluster_id = 0L,
      reasenberg_label = "background",
      reasenberg_parent_index = NA_integer_,
      reasenberg_parent_evid = NA,
      reasenberg_distance_km = NA_real_,
      reasenberg_time_days = NA_real_
    )
  
  n <- nrow(df)
  
  # Precompute event-level interaction radius and look-ahead time
  df <- df %>%
    mutate(
      reasenberg_radius_km = reasenberg_radius_km(mag, rfact = rfact),
      reasenberg_tau_days = reasenberg_tau_days(
        mag = mag,
        tau_min = tau_min,
        tau_max = tau_max,
        xmeff = xmeff,
        xk = xk
      )
    )
  
  # Union-find structure for linking events into clusters
  parent <- seq_len(n)
  
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  
  union_events <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) {
      parent[rb] <<- ra
    }
  }
  
  # Main linking loop
  for (j in 2:n) {
    
    if (j %% 1000 == 0) {
      cat("Processing event", j, "of", n, "\n")
    }
    
    t_j <- df$datetime[j]
    
    # Only search previous events within tau_max days
    dt_all <- as.numeric(difftime(t_j, df$datetime[1:(j - 1)], units = "days"))
    
    candidate_idx <- which(
      dt_all > 0 &
        dt_all <= tau_max
    )
    
    if (length(candidate_idx) == 0) next
    
    # For each candidate previous event, check event-specific look-ahead time
    candidate_idx <- candidate_idx[
      dt_all[candidate_idx] <= df$reasenberg_tau_days[candidate_idx]
    ]
    
    if (length(candidate_idx) == 0) next
    
    # Compute distance to candidate previous events
    distances <- haversine_km(
      lon1 = df$lon[candidate_idx],
      lat1 = df$lat[candidate_idx],
      lon2 = df$lon[j],
      lat2 = df$lat[j]
    )
    
    # Spatial interaction condition:
    # distance must be less than the previous event's interaction radius
    inside_space <- distances <= df$reasenberg_radius_km[candidate_idx]
    
    if (!any(inside_space)) next
    
    linked_idx <- candidate_idx[inside_space]
    linked_dist <- distances[inside_space]
    
    # Link current event to all compatible previous events
    for (k in seq_along(linked_idx)) {
      union_events(linked_idx[k], j)
    }
    
    # Store nearest compatible previous event as parent for diagnostics
    nearest_id <- which.min(linked_dist)
    nearest_parent <- linked_idx[nearest_id]
    
    df$reasenberg_parent_index[j] <- nearest_parent
    df$reasenberg_parent_evid[j] <- df$evid[nearest_parent]
    df$reasenberg_distance_km[j] <- linked_dist[nearest_id]
    df$reasenberg_time_days[j] <- as.numeric(
      difftime(df$datetime[j], df$datetime[nearest_parent], units = "days")
    )
  }
  
  # Convert union-find roots to cluster ids
  roots <- sapply(seq_len(n), find_root)
  root_table <- table(roots)
  
  clustered_roots <- names(root_table[root_table > 1])
  
  cluster_map <- data.frame(
    root = as.integer(clustered_roots),
    reasenberg_cluster_id = seq_along(clustered_roots)
  )
  
  df$root <- roots
  
  df <- df %>%
    left_join(cluster_map, by = "root", suffix = c("", "_new")) %>%
    mutate(
      reasenberg_cluster_id = ifelse(
        is.na(reasenberg_cluster_id_new),
        0L,
        reasenberg_cluster_id_new
      )
    ) %>%
    select(-reasenberg_cluster_id_new)
  
  # Label background / mainshock / foreshock / aftershock
  df <- df %>%
    group_by(reasenberg_cluster_id) %>%
    group_modify(~ {
      
      group_df <- .x
      
      if (.y$reasenberg_cluster_id == 0) {
        group_df$reasenberg_label <- "background"
        group_df$reasenberg_mainshock_evid <- NA
        return(group_df)
      }
      
      main_idx_local <- which.max(group_df$mag)
      main_time <- group_df$datetime[main_idx_local]
      main_evid <- group_df$evid[main_idx_local]
      
      group_df <- group_df %>%
        mutate(
          reasenberg_mainshock_evid = main_evid,
          reasenberg_label = case_when(
            row_number() == main_idx_local ~ "mainshock",
            datetime < main_time ~ "foreshock",
            datetime > main_time ~ "aftershock",
            TRUE ~ "clustered"
          )
        )
      
      group_df
    }) %>%
    ungroup()
  
  # For ETAS initialization:
  # background + mainshock are background-like,
  # foreshock + aftershock are clustered/triggered.
  df <- df %>%
    mutate(
      reasenberg_background_prob = ifelse(
        reasenberg_label %in% c("background", "mainshock"),
        1,
        0
      )
    ) %>%
    arrange(original_order)
  
  return(df)
}
# ============================================================
# Run Reasenberg-style declustering
# Standard parameters based on commonly used Reasenberg settings
# ============================================================

declust_reasenberg <- reasenberg_declustering(
  df = catalog_main,
  tau_min = 1,
  tau_max = 10,
  p1 = 0.95,
  xmeff = 4.0,
  xk = 0.5,
  rfact = 10
)

reasenberg_summary <- declust_reasenberg %>%
  count(reasenberg_label) %>%
  mutate(
    percentage = 100 * n / sum(n)
  )

print(reasenberg_summary)

write.csv(
  declust_reasenberg,
  "declustering_reasenberg_style_Mc_1_7.csv",
  row.names = FALSE
)

write.csv(
  reasenberg_summary,
  "declustering_reasenberg_style_summary_Mc_1_7.csv",
  row.names = FALSE
)
# ============================================================
# Initial background rate from Reasenberg-style declustering
# ============================================================

bg_rate_reasenberg_monthly <- declust_reasenberg %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarise(
    total_events = n(),
    background_events = sum(reasenberg_background_prob, na.rm = TRUE),
    background_rate = background_events,
    .groups = "drop"
  ) %>%
  mutate(method = "Reasenberg-style")

print(head(bg_rate_reasenberg_monthly))

write.csv(
  bg_rate_reasenberg_monthly,
  "initial_background_rate_reasenberg_style_monthly_Mc_1_7.csv",
  row.names = FALSE
)

p_bg_reasenberg <- ggplot(
  bg_rate_reasenberg_monthly,
  aes(x = month, y = background_rate)
) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Initial Background Rate from Reasenberg-style Declustering",
    subtitle = "Monthly background/mainshock-like event count",
    x = "Time",
    y = "Background events per month"
  ) +
  theme_minimal()

print(p_bg_reasenberg)

ggsave(
  "initial_background_rate_reasenberg_style_monthly_Mc_1_7.png",
  p_bg_reasenberg,
  width = 9,
  height = 5,
  dpi = 300
)

# ============================================================
# Compare initial background rates from three declustering methods
# ============================================================

bg_rate_three_methods <- bind_rows(
  bg_rate_gk_monthly %>%
    select(month, total_events, background_events, background_rate, method),
  
  bg_rate_nn_monthly %>%
    select(month, total_events, background_events, background_rate, method),
  
  bg_rate_reasenberg_monthly %>%
    select(month, total_events, background_events, background_rate, method)
)

p_bg_three <- ggplot(
  bg_rate_three_methods,
  aes(x = month, y = background_rate, linetype = method)
) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Comparison of Initial Background Rates",
    subtitle = "Gardner-Knopoff vs nearest-neighbour vs Reasenberg-style declustering",
    x = "Time",
    y = "Background events per month",
    linetype = "Declustering method"
  ) +
  theme_minimal()

print(p_bg_three)

ggsave(
  "initial_background_rate_three_declustering_methods_Mc_1_7.png",
  p_bg_three,
  width = 9,
  height = 5,
  dpi = 300
)

write.csv(
  bg_rate_three_methods,
  "initial_background_rate_three_declustering_methods_Mc_1_7.csv",
  row.names = FALSE
)

# ============================================================
# STEP 4: Fit temporal non-stationary ETAS model
# Using declustering-based initial background rates as starting values
# ============================================================

library(dplyr)
library(lubridate)
library(splines)
library(ggplot2)
library(purrr)

# ------------------------------------------------------------
# 1. Prepare full catalogue for ETAS fitting
# IMPORTANT: "full catalogue" here means full Mc-complete catalogue,
# not declustered catalogue.
# ------------------------------------------------------------

mc_main <- 1.7

etas_cat <- catalog %>%
  filter(!is.na(datetime), !is.na(mag)) %>%
  filter(mag >= mc_main) %>%
  arrange(datetime) %>%
  mutate(
    event_id = row_number(),
    t_days = as.numeric(difftime(datetime, min(datetime), units = "days")),
    mag_excess = mag - mc_main
  )

t_start <- min(etas_cat$t_days)
t_end <- max(etas_cat$t_days)
T_days <- t_end - t_start

cat("ETAS fitting catalogue:", nrow(etas_cat), "events\n")
cat("Time span:", round(T_days, 2), "days\n")


# ------------------------------------------------------------
# 2. Spline basis for non-stationary background rate mu(t)
# log(mu(t)) = B(t) %*% beta
# ------------------------------------------------------------

mu_df <- 8   # number of spline basis functions; can test 6, 8, 10

B_event <- splines::ns(
  etas_cat$t_days,
  df = mu_df,
  intercept = TRUE,
  Boundary.knots = c(t_start, t_end)
)

# Daily grid for numerical integration of background rate
t_grid <- seq(t_start, t_end, by = 1)

B_grid <- predict(B_event, newx = t_grid)

grid_dt <- c(diff(t_grid), 1)


# ------------------------------------------------------------
# 3. Convert monthly declustering background counts into
# starting beta values for log(mu(t))
# ------------------------------------------------------------

make_initial_beta_from_bg_rate <- function(
    bg_rate_monthly,
    B_event,
    t_origin,
    t_end_datetime,
    ridge = 1e-4
) {
  
  bg <- bg_rate_monthly %>%
    mutate(
      month_start = as.POSIXct(month),
      month_end = month_start %m+% months(1),
      month_start_clip = pmax(month_start, t_origin),
      month_end_clip = pmin(month_end, t_end_datetime),
      exposure_days = as.numeric(
        difftime(month_end_clip, month_start_clip, units = "days")
      ),
      month_mid = month_start_clip + (month_end_clip - month_start_clip) / 2,
      t_mid_days = as.numeric(difftime(month_mid, t_origin, units = "days")),
      # small offset avoids log(0)
      initial_rate_per_day = (background_events + 0.1) / pmax(exposure_days, 1),
      log_initial_rate = log(initial_rate_per_day)
    ) %>%
    filter(
      is.finite(t_mid_days),
      t_mid_days >= 0,
      is.finite(log_initial_rate)
    )
  
  B_mid <- predict(B_event, newx = bg$t_mid_days)
  
  XtX <- t(B_mid) %*% B_mid
  Xty <- t(B_mid) %*% bg$log_initial_rate
  
  beta_init <- solve(
    XtX + ridge * diag(ncol(B_mid)),
    Xty
  )
  
  as.numeric(beta_init)
}


t_origin <- min(etas_cat$datetime)
t_end_datetime <- max(etas_cat$datetime)

beta_init_gk <- make_initial_beta_from_bg_rate(
  bg_rate_gk_monthly,
  B_event = B_event,
  t_origin = t_origin,
  t_end_datetime = t_end_datetime
)

beta_init_nn <- make_initial_beta_from_bg_rate(
  bg_rate_nn_monthly,
  B_event = B_event,
  t_origin = t_origin,
  t_end_datetime = t_end_datetime
)

beta_init_reasenberg <- make_initial_beta_from_bg_rate(
  bg_rate_reasenberg_monthly,
  B_event = B_event,
  t_origin = t_origin,
  t_end_datetime = t_end_datetime
)


# ------------------------------------------------------------
# 4. Temporal ETAS negative log-likelihood
# ------------------------------------------------------------
# lambda(t_i) = mu(t_i) + sum_j K exp(alpha(M_j-Mc)) (t_i-t_j+c)^(-p)
#
# Parameter vector:
# par = c(logK, log_alpha, log_c, log_p_minus_1, beta_1, ..., beta_df)
#
# K     = exp(logK)
# alpha = exp(log_alpha)
# c     = exp(log_c)
# p     = 1 + exp(log_p_minus_1)
#
# p > 1 is imposed for finite triggering integral.
# ------------------------------------------------------------

etas_negloglik <- function(
    par,
    etas_cat,
    B_event,
    B_grid,
    grid_dt,
    t_grid,
    max_trigger_days = 3650,
    penalty_lambda = 1
) {
  
  logK <- par[1]
  log_alpha <- par[2]
  log_c <- par[3]
  log_p_minus_1 <- par[4]
  beta <- par[-c(1:4)]
  
  K <- exp(logK)
  alpha <- exp(log_alpha)
  c_par <- exp(log_c)
  p_par <- 1 + exp(log_p_minus_1)
  
  t <- etas_cat$t_days
  m_excess <- etas_cat$mag_excess
  n <- length(t)
  
  # Background intensity at event times
  mu_event <- as.numeric(exp(B_event %*% beta))
  
  # Background integral
  mu_grid <- as.numeric(exp(B_grid %*% beta))
  bg_integral <- sum(mu_grid * grid_dt)
  
  # Triggered intensity at each event
  trig_event <- numeric(n)
  
  for (i in 2:n) {
    
    dt <- t[i] - t[1:(i - 1)]
    
    keep <- dt > 0 & dt <= max_trigger_days
    
    if (!any(keep)) {
      trig_event[i] <- 0
    } else {
      prev_idx <- which(keep)
      
      trig_event[i] <- sum(
        K *
          exp(alpha * m_excess[prev_idx]) *
          (dt[prev_idx] + c_par)^(-p_par)
      )
    }
  }
  
  lambda_event <- mu_event + trig_event
  
  # Avoid log(0)
  if (any(!is.finite(lambda_event)) || any(lambda_event <= 0)) {
    return(1e100)
  }
  
  loglik_events <- sum(log(lambda_event))
  
  # Triggering integral
  # Integral from each event time to min(T, t_i + max_trigger_days)
  T_end <- max(t)
  
  trigger_integral_each <- numeric(n)
  
  for (j in 1:n) {
    
    upper_time <- min(T_end, t[j] + max_trigger_days)
    upper_dt <- upper_time - t[j]
    
    if (upper_dt <= 0) {
      trigger_integral_each[j] <- 0
    } else {
      trigger_integral_each[j] <-
        K *
        exp(alpha * m_excess[j]) *
        (
          c_par^(1 - p_par) -
            (upper_dt + c_par)^(1 - p_par)
        ) / (p_par - 1)
    }
  }
  
  trigger_integral <- sum(trigger_integral_each)
  
  loglik <- loglik_events - bg_integral - trigger_integral
  
  # Smoothness penalty for non-stationary background
  # This prevents monthly/spline background from becoming too wiggly.
  penalty <- penalty_lambda * sum(diff(beta, differences = 2)^2)
  
  negloglik <- -loglik + penalty
  
  if (!is.finite(negloglik)) {
    return(1e100)
  }
  
  return(negloglik)
}


# ------------------------------------------------------------
# 5. Fit ETAS model from one declustering-based starting value
# ------------------------------------------------------------

fit_etas_from_start <- function(
    beta_init,
    method_name,
    etas_cat,
    B_event,
    B_grid,
    grid_dt,
    t_grid,
    max_trigger_days = 3650,
    penalty_lambda = 1,
    maxit = 300
) {
  
  # Starting values for triggering parameters
  # These are starting values only; likelihood optimization updates them.
  K0 <- 0.02
  alpha0 <- 1.0
  c0 <- 0.01
  p0 <- 1.1
  
  par0 <- c(
    log(K0),
    log(alpha0),
    log(c0),
    log(p0 - 1),
    beta_init
  )
  
  lower <- c(
    log(1e-6),     # K
    log(0.05),     # alpha
    log(1e-4),     # c
    log(0.001),    # p - 1
    rep(-20, length(beta_init))
  )
  
  upper <- c(
    log(10),       # K
    log(5),        # alpha
    log(10),       # c
    log(5),        # p - 1
    rep(5, length(beta_init))
  )
  
  cat("\nFitting ETAS model using starting background from:",
      method_name, "\n")
  
  fit <- optim(
    par = par0,
    fn = etas_negloglik,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      maxit = maxit,
      trace = 1,
      REPORT = 5
    ),
    etas_cat = etas_cat,
    B_event = B_event,
    B_grid = B_grid,
    grid_dt = grid_dt,
    t_grid = t_grid,
    max_trigger_days = max_trigger_days,
    penalty_lambda = penalty_lambda
  )
  
  par_hat <- fit$par
  
  result <- list(
    method = method_name,
    fit = fit,
    par_hat = par_hat,
    K = exp(par_hat[1]),
    alpha = exp(par_hat[2]),
    c = exp(par_hat[3]),
    p = 1 + exp(par_hat[4]),
    beta = par_hat[-c(1:4)],
    negloglik = fit$value,
    convergence = fit$convergence,
    message = fit$message
  )
  
  return(result)
}


# ------------------------------------------------------------
# 6. Fit three ETAS models using three different initial values
# ------------------------------------------------------------

fit_gk <- fit_etas_from_start(
  beta_init = beta_init_gk,
  method_name = "Gardner-Knopoff",
  etas_cat = etas_cat,
  B_event = B_event,
  B_grid = B_grid,
  grid_dt = grid_dt,
  t_grid = t_grid,
  max_trigger_days = 3650,
  penalty_lambda = 1,
  maxit = 300
)

fit_nn <- fit_etas_from_start(
  beta_init = beta_init_nn,
  method_name = "Nearest-neighbour",
  etas_cat = etas_cat,
  B_event = B_event,
  B_grid = B_grid,
  grid_dt = grid_dt,
  t_grid = t_grid,
  max_trigger_days = 3650,
  penalty_lambda = 1,
  maxit = 300
)

fit_reasenberg <- fit_etas_from_start(
  beta_init = beta_init_reasenberg,
  method_name = "Reasenberg-style",
  etas_cat = etas_cat,
  B_event = B_event,
  B_grid = B_grid,
  grid_dt = grid_dt,
  t_grid = t_grid,
  max_trigger_days = 3650,
  penalty_lambda = 1,
  maxit = 300
)


# ------------------------------------------------------------
# 7. Compare convergence and fitted ETAS parameters
# ------------------------------------------------------------

etas_fit_summary <- data.frame(
  method = c(
    fit_gk$method,
    fit_nn$method,
    fit_reasenberg$method
  ),
  convergence = c(
    fit_gk$convergence,
    fit_nn$convergence,
    fit_reasenberg$convergence
  ),
  negloglik = c(
    fit_gk$negloglik,
    fit_nn$negloglik,
    fit_reasenberg$negloglik
  ),
  K = c(fit_gk$K, fit_nn$K, fit_reasenberg$K),
  alpha = c(fit_gk$alpha, fit_nn$alpha, fit_reasenberg$alpha),
  c = c(fit_gk$c, fit_nn$c, fit_reasenberg$c),
  p = c(fit_gk$p, fit_nn$p, fit_reasenberg$p)
)

print(etas_fit_summary)

write.csv(
  etas_fit_summary,
  "nonstationary_ETAS_fit_summary_three_initializations_Mc_1_7.csv",
  row.names = FALSE
)
