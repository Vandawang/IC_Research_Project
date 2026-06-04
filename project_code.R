packages <- c("dplyr", "ggplot2", "lubridate", "readr", "tidyr", "purrr",
              "scales", "maps", "viridis")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

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

# check data and some exploratory analysis
cat("Number of events:", nrow(catalog), "\n")
cat("Time range:", as.character(min(catalog$datetime)), "to",
    as.character(max(catalog$datetime)), "\n")
cat("Magnitude range:", min(catalog$mag), "to", max(catalog$mag), "\n")

head(catalog)
summary(catalog$mag)

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
  "magnitude_over_time_Mc_1_7.png",
  p_mag_time,
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
