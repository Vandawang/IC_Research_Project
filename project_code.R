###########################################################################
# install packages
packages <- c("dplyr", "ggplot2", "lubridate", "readr", "tidyr", "purrr",
              "scales", "maps", "viridis","SAPP","mclust")
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)}}

#############################################################################
# read data
data <- readLines("SearchResults.txt", warn = FALSE)
data <- data[grepl("^\\d{4}/\\d{2}/\\d{2}", data)]
catalog <- read.table(text = data,header = FALSE,stringsAsFactors = FALSE)
names(catalog) <- c("date", "time", "event_type", "geo_type","mag", "mag_type", 
                    "lat", "lon", "depth","quality", "evid", "nph", "ngrm")
catalog <- catalog %>%
  mutate(datetime = ymd_hms(paste(date, time)),mag = as.numeric(mag),
         lat = as.numeric(lat),lon = as.numeric(lon),depth = as.numeric(depth)) %>%
  filter(!is.na(mag)) %>%
  arrange(datetime)
# check data and some exploratory analysis
cat("Number of events:", nrow(catalog), "\n")
cat("Time range:", as.character(min(catalog$datetime)), "to",
    as.character(max(catalog$datetime)), "\n")
cat("Magnitude range:", min(catalog$mag), "to", max(catalog$mag), "\n")

head(catalog)
summary(catalog$mag)
###############################################################################
# step 1: Estimate catalogue completeness magnitude Mc
# Methods:
# 1. FMD plot
# 2. MAXC
# 3. MBS-WW
# 4. GFT-95% / GFT-90%
# 5. Mc sensitivity thresholds

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
        sum(mags >= m & mags < m + bin_width)}),
      cumulative_count = sapply(mag_bin, function(m) {
        sum(mags >= m)})) %>%
    filter(cumulative_count > 0)
  return(fmd)}

fmd <- make_fmd(catalog$mag, bin_width)

# Incremental FMD
p1 <- ggplot(fmd, aes(x = mag_bin, y = count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(title = "Incremental Frequency-Magnitude Distribution",
       x = "Magnitude",y = "Number of events") +
  theme_minimal()
# Cumulative FMD
p2 <- ggplot(fmd, aes(x = mag_bin, y = cumulative_count)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  labs(title = "Cumulative Frequency-Magnitude Distribution",
       x = "Magnitude",y = "Cumulative number of events") +
  theme_minimal()

print(p1)
print(p2)

# MAXC method (reference)
estimate_maxc <- function(mags, bin_width = 0.1) {
  fmd <- make_fmd(mags, bin_width)
  maxc <- fmd$mag_bin[which.max(fmd$count)]
  return(maxc)}
mc_maxc <- estimate_maxc(catalog$mag, bin_width)
cat("MAXC Mc =", mc_maxc, "\n")

#  b-value MLE with magnitude bin correction （used for GFT and MBS-WW）
estimate_b_value <- function(mags, mc, bin_width = 0.1) {
  x <- mags[mags >= mc]
  n <- length(x)
  if (n < 30) {
    return(data.frame(mc = mc,n = n,b = NA,delta_b = NA))}
  mean_mag <- mean(x)
  # Utsu / Aki MLE with bin correction
  b <- log10(exp(1)) / (mean_mag - (mc - bin_width / 2))
  # Shi & Bolt uncertainty
  delta_b <- 2.3 * b^2 * sqrt(sum((x - mean_mag)^2) / (n * (n - 1)))
  return(data.frame(mc = mc,n = n,b = b,delta_b = delta_b))}

# MBS-WW method
estimate_mbs_ww <- function(mags, bin_width = 0.1, delta_m = 0.5) {
  mags <- mags[!is.na(mags)]
  candidates <- seq(floor(min(mags) / bin_width) * bin_width,
                    floor(max(mags) / bin_width) * bin_width - delta_m,
                    by = bin_width)
  b_table <- purrr::map_dfr(
    candidates,
    ~ estimate_b_value(mags, .x, bin_width))
  b_table <- b_table %>%
    filter(!is.na(b))
  mbs_table <- b_table %>%
    rowwise() %>%
    mutate(b_avg = mean(b_table$b[b_table$mc >= mc & b_table$mc <= mc + delta_m],
                        na.rm = TRUE),
           abs_diff = abs(b - b_avg),
           stable = abs_diff <= delta_b) %>%
    ungroup()

  mc_est <- mbs_table %>%
    filter(stable == TRUE) %>%
    slice(1) %>%
    pull(mc)
  
  if (length(mc_est) == 0) mc_est <- NA
  
  return(list(mc = mc_est,table = mbs_table))}

mbs <- estimate_mbs_ww(catalog$mag, bin_width, delta_m = 0.5)
mc_mbs <- mbs$mc
cat("MBS-WW Mc =", mc_mbs, "\n")

# Plot b-value stability
ggplot(mbs$table, aes(x = mc, y = b)) +
  geom_line() +
  geom_point() +
  geom_line(aes(y = b_avg), linetype = "dashed") +
  geom_errorbar(aes(ymin = b - delta_b, ymax = b + delta_b),width = 0.03,alpha = 0.4) +
  geom_vline(xintercept = mc_mbs, linetype = "dashed") +
  labs(title = "MBS-WW b-value Stability Method",x = "Candidate Mc",y = "b-value") +
  theme_minimal()

# GFT check at MBS-WW Mc
gft_at_mc <- function(mags, mc, bin_width = 0.1) {
  mags <- mags[!is.na(mags)]
  x <- mags[mags >= mc]
  if (length(x) < 50) {
    warning("Too few events above this Mc.")
    return(data.frame(mc = mc,n = length(x),b = NA,gft = NA))}
  
  # Estimate b-value at this Mc
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  
  if (is.na(b)) {
    warning("b-value could not be estimated.")
    return(data.frame(mc = mc,n = length(x),b = NA,gft = NA))}
  
  # Magnitude bins above Mc
  mag_bins <- seq(mc, max(mags), by = bin_width)
  # Observed cumulative FMD
  obs <- sapply(mag_bins, function(m) {sum(mags >= m)})
  # Predicted cumulative FMD under Gutenberg-Richter law
  pred <- obs[1] * 10^(-b * (mag_bins - mc))
  # Goodness-of-fit value
  gft <- 100 - (sum(abs(obs - pred)) / sum(obs)) * 100
  # Return table for this Mc
  result <- data.frame(mc = mc,n = length(x),b = b,gft = gft)
  
  return(result)
}

# Use MBS-WW result
mc_mbs <- 1.7
gft_check_mbs <- gft_at_mc(mags = catalog$mag,mc = mc_mbs,bin_width = bin_width)
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
  obs <- sapply(mag_bins, function(m) {sum(mags >= m)})
  pred <- obs[1] * 10^(-b * (mag_bins - mc))
  plot_data <- data.frame(mag_bin = mag_bins,observed = obs,predicted = pred)
  
  ggplot(plot_data, aes(x = mag_bin)) +
    geom_point(aes(y = observed), size = 2) +
    geom_line(aes(y = observed), linewidth = 0.8) +
    geom_line(aes(y = predicted), linetype = "dashed", linewidth = 0.9) +
    scale_y_log10() +
    geom_vline(xintercept = mc, linetype = "dotted") +
    labs(title = paste0("GFT Check at MBS-WW Mc = ", mc),
         subtitle = paste0("Observed cumulative FMD vs Gutenberg-Richter prediction; b = ",
        round(b, 3)),
        x = "Magnitude",y = "Cumulative number of events") +
    theme_minimal()
  }

plot_gft_at_mc(mags = catalog$mag,mc = 1.7,bin_width = bin_width)
##############################################################################
# Monte Carlo envelope for GFT / cumulative FMD
# Test whether high-magnitude deviations are random tail effects

gft_monte_carlo_envelope <- function(mags,mc = 1.7,bin_width = 0.1,n_sim = 1000,
                                     conf = 0.95,seed = 123) {
  set.seed(seed)
  mags <- mags[!is.na(mags)]
  mags_mc <- mags[mags >= mc]
  n <- length(mags_mc)
  b_info <- estimate_b_value(mags, mc, bin_width)
  b <- b_info$b
  mag_bins <- seq(mc, max(mags_mc), by = bin_width)
  observed <- sapply(mag_bins, function(m) {sum(mags_mc >= m)})
  
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
  
  envelope <- data.frame(mag_bin = mag_bins,observed = observed,
                         predicted = predicted,
                         lower = apply(sim_cum, 1, quantile, probs = lower_prob, na.rm = TRUE),
                         upper = apply(sim_cum, 1, quantile, probs = upper_prob, na.rm = TRUE))
  
  # Check whether observed values are inside the simulation envelope
  envelope <- envelope %>%
    mutate(inside_envelope = observed >= lower & observed <= upper,
           expected_count = predicted)
  
  summary <- envelope %>%
    summarise(mc = mc,n = n,b = b,n_bins = n(),bins_inside = sum(inside_envelope),
              proportion_inside = mean(inside_envelope),bins_outside = sum(!inside_envelope))
  
  return(list(envelope = envelope,summary = summary,b_info = b_info))
}


mc_env_1_7 <- gft_monte_carlo_envelope(mags = catalog$mag,mc = 1.7,bin_width = bin_width,
                                       n_sim = 1000,conf = 0.95,seed = 123)
print(mc_env_1_7$summary)
##############################################################################
# ensure if the large number of earthquakes in 2019 will affect the choice of Mc
# ============================================================
# Split catalogue into 2019 and non-2019 groups
# Purpose: test whether 2019 earthquake sequence controls Mc/FMD

catalog_fmd_grouped <- catalog %>%
  filter(!is.na(mag), !is.na(datetime)) %>%
  mutate(year = lubridate::year(datetime),
         period_group = ifelse(year == 2019, "2019", "non-2019"))
group_count_summary <- catalog_fmd_grouped %>%
  count(period_group) %>%
  mutate(percentage = 100 * n / sum(n))
print(group_count_summary)

write.csv(group_count_summary,"FMD_group_event_counts_2019_vs_non2019.csv",
          row.names = FALSE)
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
    mutate(mag_bin = floor(mag / bin_width) * bin_width) %>%
    count(period_group, mag_bin, name = "incremental_count") %>%
    tidyr::complete(period_group,mag_bin = mag_bins,fill = list(incremental_count = 0)) %>%
    group_by(period_group) %>%
    arrange(mag_bin, .by_group = TRUE) %>%
    mutate(cumulative_count = purrr::map_dbl( mag_bin,
                                              ~ sum(incremental_count[mag_bin >= .x]))) %>%
    ungroup()
  return(fmd_grouped)
}

fmd_2019_compare <- make_fmd_by_group(df = catalog_fmd_grouped,bin_width = bin_width)
p_incremental_2019 <- ggplot(fmd_2019_compare,
                             aes(x = mag_bin, y = incremental_count, linetype = period_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(title = "Incremental FMD: 2019 vs non-2019",
       subtitle = "Dotted line indicates Mc = 1.7",x = "Magnitude",
       y = "Incremental number of events",linetype = "Period") +
  theme_minimal()
print(p_incremental_2019)

ggsave("Incremental_FMD_2019_vs_non2019.png",p_incremental_2019,width = 8,height = 5,
       dpi = 300)
# ============================================================
# Cumulative FMD comparison: 2019 vs non-2019
# ============================================================
p_cumulative_2019 <- ggplot(fmd_2019_compare,
                            aes(x = mag_bin, y = cumulative_count, linetype = period_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(title = "Cumulative FMD: 2019 vs non-2019",
       subtitle = "Dotted line indicates Mc = 1.7",x = "Magnitude",
       y = "Cumulative number of events",linetype = "Period") +
  theme_minimal()
print(p_cumulative_2019)
ggsave("Cumulative_FMD_2019_vs_non2019.png",p_cumulative_2019,width = 8,height = 5,dpi = 300)
# ============================================================
# Normalized cumulative FMD comparison
# This compares shape rather than total number of events
# ============================================================
fmd_2019_compare_norm <- fmd_2019_compare %>%
  group_by(period_group) %>%
  mutate(cumulative_norm = cumulative_count / max(cumulative_count, na.rm = TRUE)) %>%
  ungroup()

p_cumulative_norm_2019 <- ggplot(fmd_2019_compare_norm,
                                 aes(x = mag_bin, y = cumulative_norm, linetype = period_group)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 1.7, linetype = "dotted") +
  scale_y_log10() +
  labs(title = "Normalized Cumulative FMD: 2019 vs non-2019",
       subtitle = "Curves are normalized by each group's total event count",
       x = "Magnitude",y = "Normalized cumulative proportion",linetype = "Period") +
  theme_minimal()
print(p_cumulative_norm_2019)
ggsave("Normalized_Cumulative_FMD_2019_vs_non2019.png",p_cumulative_norm_2019,width = 8,
       height = 5,dpi = 300)
# ============================================================
# b-value and GFT comparison at Mc = 1.7
# ============================================================
mc_main <- 1.7
fmd_group_quality <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    mags_group <- .x$mag
    b_info <- estimate_b_value(mags = mags_group, mc = mc_main,bin_width = bin_width)
    gft_info <- gft_at_mc(mags = mags_group, mc = mc_main,bin_width = bin_width)
    data.frame(n_total = length(mags_group),
               n_above_mc = sum(mags_group >= mc_main, na.rm = TRUE),
               mc = mc_main,b = b_info$b,delta_b = b_info$delta_b,gft = gft_info$gft)
  }) %>%
  ungroup()

print(fmd_group_quality)
write.csv(fmd_group_quality,"FMD_quality_2019_vs_non2019_at_Mc_1_7.csv",row.names = FALSE)

# ============================================================
# Estimate MBS-WW Mc separately for 2019 and non-2019
# ============================================================
mc_mbs_by_group <- catalog_fmd_grouped %>%
  group_by(period_group) %>%
  group_modify(~ {
    mags_group <- .x$mag
    mbs_group <- estimate_mbs_ww(mags = mags_group,bin_width = bin_width,delta_m = 0.5)
    data.frame(n_total = length(mags_group), mc_mbs = mbs_group$mc)
  }) %>%
  ungroup()
print(mc_mbs_by_group)
write.csv( mc_mbs_by_group,"MBS_WW_Mc_2019_vs_non2019.csv",row.names = FALSE)

# ============================================================
# Temporal Mc analysis using rolling event windows
# check if mc changed a lot during 2019
# ============================================================
window_size <- 1000     # number of earthquakes per window
step_size   <- 250      # move forward by 250 events each time
bin_width   <- 0.1

# Optional fixed Mc values for comparison
mc_check_values <- c(1.7, 2.0)

catalog_temporal <- catalog %>%
  filter(!is.na(datetime), !is.na(mag)) %>%
  arrange(datetime)

n_events <- nrow(catalog_temporal)
window_starts <- seq(from = 1,to = n_events - window_size + 1,by = step_size)
# Estimate Mc and b-value in each rolling window
temporal_mc_results <- purrr::map_dfr(window_starts,function(start_idx) {
    end_idx <- start_idx + window_size - 1
    temp <- catalog_temporal[start_idx:end_idx, ]
    mags <- temp$mag
    start_time  <- min(temp$datetime)
    end_time    <- max(temp$datetime)
    center_time <- start_time + (end_time - start_time) / 2
    mc_maxc_temp <- estimate_maxc(mags = mags,bin_width = bin_width)
    mbs_temp <- estimate_mbs_ww(mags = mags,bin_width = bin_width,delta_m = 0.5)
    mc_mbs_temp <- mbs_temp$mc
    
    b_17 <- estimate_b_value(mags = mags,mc = 1.7, bin_width = bin_width)

    b_20 <- estimate_b_value(mags = mags,mc = 2.0,bin_width = bin_width)
    
    data.frame(start_time = start_time,center_time = center_time,end_time = end_time,
               n_total = length(mags), mc_maxc = mc_maxc_temp, mc_mbs = mc_mbs_temp,
               n_above_1_7 = b_17$n, b_at_1_7 = b_17$b,delta_b_1_7 = b_17$delta_b,
               n_above_2_0 = b_20$n,b_at_2_0 = b_20$b, delta_b_2_0 = b_20$delta_b)
  }
)

print(head(temporal_mc_results))

# Plot temporal Mc
temporal_mc_long <- temporal_mc_results %>%
  select( center_time,mc_maxc, mc_mbs) %>%
  pivot_longer( cols = c(mc_maxc, mc_mbs),names_to = "method",values_to = "mc")

p_temporal_mc <- ggplot(temporal_mc_long,
                        aes(x = center_time, y = mc, linetype = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 1.7,linetype = "dotted" ) +
  geom_hline( yintercept = 2.0, linetype = "dashed") +
  labs(title = "Temporal Variation of Catalogue Completeness Magnitude",
       subtitle = paste0("Rolling window: ", window_size," events; step: ", step_size,
                         " events" ),
       x = "Time",
       y = "Estimated Mc",
       linetype = "Method") +
  theme_minimal()

print(p_temporal_mc)

p_temporal_mc_2019 <- ggplot(temporal_mc_long,
                             aes(x = center_time,y = mc,linetype = method )) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  annotate( "rect",xmin = as.POSIXct("2019-01-01", tz = "UTC"),
            xmax = as.POSIXct("2019-12-31 23:59:59", tz = "UTC"), ymin = -Inf, ymax = Inf,
            alpha = 0.08) +
  geom_hline(yintercept = 1.7,linetype = "dotted") +
  geom_hline( yintercept = 2.0, linetype = "dashed" ) +
  labs(title = "Temporal Mc with 2019 Highlighted",
       subtitle = "Shaded interval indicates 2019",
       x = "Time",y = "Estimated Mc",linetype = "Method" ) +
  theme_minimal()

print(p_temporal_mc_2019)

# Temporal b-value plot
p_b_temporal <- ggplot(temporal_mc_results,
                       aes( x = center_time,y = b_at_2_0)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = b_at_2_0 - delta_b_2_0, ymax = b_at_2_0 + delta_b_2_0),
                width = 10, alpha = 0.4) +
  annotate( "rect",xmin = as.POSIXct("2019-01-01", tz = "UTC"),
            xmax = as.POSIXct("2019-12-31 23:59:59", tz = "UTC"),ymin = -Inf,
            ymax = Inf,alpha = 0.08) +
  labs( title = "Temporal Variation of b-value",
        subtitle = "b-value estimated using events with M >= 2.0",x = "Time",y = "b-value" ) +
  theme_minimal()

print(p_b_temporal)

write.csv( temporal_mc_results,"Temporal_Mc_Rolling_Window.csv",row.names = FALSE)

ggsave("Temporal_Mc_Rolling_Window.png", p_temporal_mc,width = 10,height = 5,dpi = 300)

ggsave("Temporal_Mc_2019_Highlight.png", p_temporal_mc_2019, width = 10,height = 5, dpi = 300)

ggsave("Temporal_b_value_Mc_2_0.png",p_b_temporal,width = 10, height = 5,dpi = 300)
temporal_mc_results %>%
  filter( year(center_time) == 2019,mc_maxc >= 1.7 | mc_mbs >= 1.7) %>%
  select(start_time,center_time, end_time,mc_maxc, mc_mbs,n_above_1_7, b_at_1_7, n_above_2_0,
         b_at_2_0) %>%
  arrange(center_time)

#################################################################
##############################################################################
# sensitivity analysis
mc_main <- 2.0

mc_sensitivity <- c(2.0, 2.5)

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
# Main catalogue: Mc = 2.0
# ============================================================

# Create main catalogue for descriptive analysis
mc_main <- 2.0
catalog_main <- catalog %>%
  filter(mag >= mc_main) %>%
  mutate(date_only = as.Date(datetime),year = year(datetime),
         month = floor_date(datetime, "month"), year_month = format(datetime, "%Y-%m")) %>%
  arrange(datetime)

cat("Main catalogue Mc =", mc_main, "\n")
cat("Number of events in main catalogue:", nrow(catalog_main), "\n")
cat("Time range:", as.character(min(catalog_main$datetime)), "to",
    as.character(max(catalog_main$datetime)), "\n")
cat("Magnitude range:", min(catalog_main$mag), "to", max(catalog_main$mag), "\n")
cat("Depth range:", min(catalog_main$depth), "to", max(catalog_main$depth), "\n")

# Catalogue overview summary

catalog_summary <- data.frame(
  statistic = c("Mc threshold","Number of events","Start time","End time",
                "Minimum magnitude", "Maximum magnitude","Mean magnitude",
                "Median magnitude","Minimum depth","Maximum depth","Mean depth",
                "Median depth","Minimum latitude","Maximum latitude","Minimum longitude",
                "Maximum longitude"),
  value = c(mc_main,nrow(catalog_main),as.character(min(catalog_main$datetime, na.rm = TRUE)),
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
            round(max(catalog_main$lon, na.rm = TRUE), 3)))
print(catalog_summary)

write.csv(catalog_summary,"descriptive_catalog_summary_Mc_2.csv",row.names = FALSE)

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
  "annual_event_counts_Mc_2.csv",
  row.names = FALSE
)

p_annual <- ggplot(annual_counts, aes(x = year, y = n_events)) +
  geom_col() +
  labs(
    title = "Annual Number of Earthquakes",
    subtitle = "Main catalogue: M >= 2.0",
    x = "Year",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_annual)

ggsave(
  "annual_event_counts_Mc_2.png",
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
  "monthly_event_counts_Mc_2.csv",
  row.names = FALSE
)

p_monthly <- ggplot(monthly_counts, aes(x = month, y = n_events)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  labs(
    title = "Monthly Number of Earthquakes",
    subtitle = "Main catalogue: M >= 2",
    x = "Time",
    y = "Number of events"
  ) +
  theme_minimal()

print(p_monthly)

ggsave(
  "monthly_event_counts_Mc_2.png",
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
  "cumulative_event_count_Mc_2.png",
  p_cumulative,
  width = 9,
  height = 5,
  dpi = 300
)

# ============================================================
# MODEL COMPARISON 1
# Homogeneous Poisson baseline
# Main catalogue: Mc = 2.0
# ============================================================
catalog_model <- catalog_main %>%
  filter(!is.na(datetime), !is.na(mag)) %>%
  arrange(datetime)

# ------------------------------------------------------------
# Convert event times to days from catalogue start
# ------------------------------------------------------------

t0 <- min(catalog_model$datetime)
t_end <- max(catalog_model$datetime)

catalog_model <- catalog_model %>%
  mutate(
    time_days = as.numeric(
      difftime(datetime, t0, units = "days")
    )
  )

N <- nrow(catalog_model)

T_days <- as.numeric(
  difftime(t_end, t0, units = "days")
)

# ------------------------------------------------------------
# Homogeneous Poisson MLE
# ------------------------------------------------------------

lambda_hat <- N / T_days

logLik_poisson <-
  N * log(lambda_hat) -
  lambda_hat * T_days

# Number of free parameters
k_poisson <- 1

AIC_poisson <-
  -2 * logLik_poisson +
  2 * k_poisson

BIC_poisson <-
  -2 * logLik_poisson +
  log(N) * k_poisson

poisson_results <- data.frame(
  model = "Homogeneous Poisson",
  n_events = N,
  duration_days = T_days,
  lambda_per_day = lambda_hat,
  logLik = logLik_poisson,
  k = k_poisson,
  AIC = AIC_poisson,
  BIC = BIC_poisson
)

print(poisson_results)

write.csv(
  poisson_results,
  "model_comparison_poisson_Mc_2.csv",
  row.names = FALSE
)

# ============================================================
# MODEL COMPARISON 2
# Stationary temporal ETAS using SAPP
# ============================================================
# ------------------------------------------------------------
# Prepare event times
# SAPP needs event times in a numeric time scale
# ------------------------------------------------------------

etas_time <- catalog_model$time_days
etas_mag  <- catalog_model$mag

# Mc
mc_etas <- 2.0

# Magnitudes expressed relative to Mc
# depending on your exact ETAS parameterization,
# reference is normally set consistently with threshold
threshold_mag <- mc_etas
reference_mag <- mc_etas

# ------------------------------------------------------------
# Initial parameters
# IMPORTANT:
# These are starting values only.
# Replace later with initialization from your declustering scheme.
# ------------------------------------------------------------

param_initial <- c(
  0.5,    # mu/background parameter -- example only
  0.05,   # K
  0.01,   # c
  1.1,    # alpha or productivity-related parameter
  1.1     # p
)

# ------------------------------------------------------------
# Observation interval
# ------------------------------------------------------------

t_start <- min(etas_time)
t_end   <- max(etas_time)

# ------------------------------------------------------------
# Fit stationary ETAS
# ------------------------------------------------------------

fit_etas_stationary <- etasap(
  time = etas_time,
  mag = etas_mag,
  threshold = threshold_mag,
  reference = reference_mag,
  parami = param_initial,
  tstart = t_start,
  zte = t_end,
  approx = 2,
  plot = TRUE
)

print(fit_etas_stationary)
# ============================================================
# Extract stationary ETAS results
# ============================================================

# etasap returns NEGATIVE maximum log-likelihood
logLik_stationary <- -fit_etas_stationary$ngmle

# Five free parameters:
# mu, K, c, alpha, p
k_stationary <- 5

AIC_stationary <-
  -2 * logLik_stationary +
  2 * k_stationary

BIC_stationary <-
  -2 * logLik_stationary +
  log(N) * k_stationary

stationary_results <- data.frame(
  model = "Stationary ETAS",
  logLik = logLik_stationary,
  k = k_stationary,
  AIC = AIC_stationary,
  BIC = BIC_stationary
)

print(stationary_results)

print(fit_etas_stationary$param)
# ============================================================
# AIC / BIC model comparison
# ============================================================
model_comparison_basic <- data.frame(
  model = c(
    "Homogeneous Poisson",
    "Stationary ETAS"
  ),
  
  logLik = c(
    logLik_poisson,
    logLik_stationary
  ),
  
  k = c(
    k_poisson,
    k_stationary
  )
) %>%
  mutate(
    AIC = -2 * logLik + 2 * k,
    BIC = -2 * logLik + log(N) * k,
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)
  )

print(model_comparison_basic)


# ============================================================
# MODEL FORM SELECTION
# Stationary ETAS vs Non-stationary ETAS
# Main catalogue: Mc = 2.0
# ============================================================
catalog_model <- catalog_main %>%
  filter(
    !is.na(datetime),
    !is.na(mag)
  ) %>%
  arrange(datetime)

mc_etas <- 2.0

t0 <- min(catalog_model$datetime)

catalog_model <- catalog_model %>%
  mutate(
    time_days = as.numeric(
      difftime(datetime, t0, units = "days")
    )
  )

times <- catalog_model$time_days
mags  <- catalog_model$mag

N <- length(times)
T_end <- max(times)

cat("N =", N, "\n")
cat("Duration =", T_end, "days\n")

# ============================================================
# ETAS triggering integral
# ============================================================
trigger_integral <- function(
    times,
    mags,
    T_end,
    K,
    c,
    alpha,
    p,
    Mc
) {
  
  dt <- T_end - times
  
  productivity <-
    K * exp(alpha * (mags - Mc))
  
  if (abs(p - 1) < 1e-6) {
    
    integral <-
      productivity *
      log((dt + c) / c)
    
  } else {
    
    integral <-
      productivity *
      (
        (dt + c)^(1 - p) -
          c^(1 - p)
      ) /
      (1 - p)
  }
  
  sum(integral)
}
# ============================================================
# Stationary ETAS negative log-likelihood
# ============================================================
etas_nll_stationary <- function(
    par,
    times,
    mags,
    T_end,
    Mc
) {
  
  # Positive parameter transformations
  mu    <- exp(par[1])
  K     <- exp(par[2])
  c     <- exp(par[3])
  alpha <- exp(par[4])
  p     <- exp(par[5])
  
  N <- length(times)
  
  lambda <- numeric(N)
  
  for (j in seq_len(N)) {
    
    if (j == 1) {
      
      lambda[j] <- mu
      
    } else {
      
      dt <- times[j] - times[1:(j - 1)]
      
      triggering <-
        K *
        exp(
          alpha *
            (mags[1:(j - 1)] - Mc)
        ) /
        (dt + c)^p
      
      lambda[j] <-
        mu +
        sum(triggering)
    }
  }
  
  # Invalid intensity
  if (
    any(!is.finite(lambda)) ||
    any(lambda <= 0)
  ) {
    return(1e100)
  }
  
  # Integrated background rate
  integrated_background <-
    mu * T_end
  
  # Integrated triggering
  integrated_triggering <-
    trigger_integral(
      times = times,
      mags = mags,
      T_end = T_end,
      K = K,
      c = c,
      alpha = alpha,
      p = p,
      Mc = Mc
    )
  
  logLik <-
    sum(log(lambda)) -
    integrated_background -
    integrated_triggering
  
  return(-logLik)
}
# ============================================================
# Fit stationary ETAS
# ============================================================

initial_stationary <- log(
  c(
    mu = N / T_end * 0.2,
    K = 0.05,
    c = 0.01,
    alpha = 1.0,
    p = 1.1
  )
)

fit_stationary <- optim(
  par = initial_stationary,
  fn = etas_nll_stationary,
  times = times,
  mags = mags,
  T_end = T_end,
  Mc = mc_etas,
  method = "BFGS",
  control = list(
    maxit = 1000,
    reltol = 1e-8
  ),
  hessian = TRUE
)

fit_stationary$convergence
fit_stationary$value
stationary_parameters <-
  exp(fit_stationary$par)

names(stationary_parameters) <-
  c(
    "mu",
    "K",
    "c",
    "alpha",
    "p"
  )

print(stationary_parameters)
logLik_stationary_custom <-
  -fit_stationary$value

k_stationary_custom <- 5
# ============================================================
# Pre-defined equal-width time segments
# Example: approximately 5-year intervals
# ============================================================
segment_years <- 5

segment_days <-
  segment_years * 365.25

segment_breaks <- seq(
  from = 0,
  to = ceiling(T_end / segment_days) *
    segment_days,
  by = segment_days
)

# Ensure final boundary covers catalogue
if (
  tail(segment_breaks, 1) < T_end
) {
  segment_breaks <-
    c(segment_breaks, T_end)
}

J <- length(segment_breaks) - 1

cat(
  "Number of background-rate segments =",
  J,
  "\n"
)

print(segment_breaks)

segment_id <- findInterval(
  times,
  vec = segment_breaks,
  rightmost.closed = TRUE
)

segment_id[
  segment_id > J
] <- J

# ============================================================
# Non-stationary ETAS
# Piecewise-constant background rate
# ============================================================

etas_nll_nonstationary <- function(
    par,
    times,
    mags,
    T_end,
    Mc,
    segment_breaks,
    segment_id
) {
  
  J <- length(segment_breaks) - 1
  
  # Background rates for each segment
  mu_segments <-
    exp(par[1:J])
  
  # Common ETAS triggering parameters
  K <-
    exp(par[J + 1])
  
  c <-
    exp(par[J + 2])
  
  alpha <-
    exp(par[J + 3])
  
  p <-
    exp(par[J + 4])
  
  N <- length(times)
  
  lambda <- numeric(N)
  
  for (j in seq_len(N)) {
    
    mu_t <-
      mu_segments[
        segment_id[j]
      ]
    
    if (j == 1) {
      
      lambda[j] <- mu_t
      
    } else {
      
      dt <-
        times[j] -
        times[1:(j - 1)]
      
      triggering <-
        K *
        exp(
          alpha *
            (
              mags[1:(j - 1)] -
                Mc
            )
        ) /
        (dt + c)^p
      
      lambda[j] <-
        mu_t +
        sum(triggering)
    }
  }
  
  if (
    any(!is.finite(lambda)) ||
    any(lambda <= 0)
  ) {
    return(1e100)
  }
  
  # ----------------------------------------------------------
  # Integrated background intensity
  # ----------------------------------------------------------
  
  segment_lengths <-
    diff(segment_breaks)
  
  # Last segment may extend beyond catalogue end
  segment_lengths[J] <-
    T_end -
    segment_breaks[J]
  
  integrated_background <-
    sum(
      mu_segments *
        segment_lengths
    )
  
  # ----------------------------------------------------------
  # Integrated triggering intensity
  # ----------------------------------------------------------
  
  integrated_triggering <-
    trigger_integral(
      times = times,
      mags = mags,
      T_end = T_end,
      K = K,
      c = c,
      alpha = alpha,
      p = p,
      Mc = Mc
    )
  
  logLik <-
    sum(log(lambda)) -
    integrated_background -
    integrated_triggering
  
  return(-logLik)
}
# ============================================================
# Initial background rates
# ============================================================

segment_counts <- table(
  factor(
    segment_id,
    levels = 1:J
  )
)

segment_lengths <-
  diff(segment_breaks)

segment_lengths[J] <-
  T_end -
  segment_breaks[J]

raw_segment_rates <-
  as.numeric(segment_counts) /
  segment_lengths

initial_mu_segments <-
  pmax(
    raw_segment_rates * 0.2,
    1e-6
  )

initial_nonstationary <- c(
  log(initial_mu_segments),
  log(0.05),  # K
  log(0.01),  # c
  log(1.0),   # alpha
  log(1.1)    # p
)
# ============================================================
# Fit non-stationary ETAS
# ============================================================

fit_nonstationary <- optim(
  par = initial_nonstationary,
  fn = etas_nll_nonstationary,
  times = times,
  mags = mags,
  T_end = T_end,
  Mc = mc_etas,
  segment_breaks = segment_breaks,
  segment_id = segment_id,
  method = "BFGS",
  control = list(
    maxit = 1500,
    reltol = 1e-8
  ),
  hessian = TRUE
)

fit_nonstationary$convergence
fit_nonstationary$value

J <- length(segment_breaks) - 1

mu_nonstationary <-
  exp(
    fit_nonstationary$par[
      1:J
    ]
  )

trigger_par_nonstationary <-
  exp(
    fit_nonstationary$par[
      (J + 1):(J + 4)
    ]
  )

names(
  trigger_par_nonstationary
) <- c(
  "K",
  "c",
  "alpha",
  "p"
)

print(mu_nonstationary)
print(trigger_par_nonstationary)

logLik_nonstationary <-
  -fit_nonstationary$value

k_nonstationary <-
  J + 4

cat(
  "Stationary ETAS logLik =",
  logLik_stationary_custom,
  "\n"
)

cat(
  "Non-stationary ETAS logLik =",
  logLik_nonstationary,
  "\n"
)

cat(
  "Stationary k =",
  k_stationary_custom,
  "\n"
)

cat(
  "Non-stationary k =",
  k_nonstationary,
  "\n"
)

# ============================================================
# Stationary vs Non-stationary ETAS
# ============================================================

etas_model_comparison <- data.frame(
  model = c(
    "Stationary ETAS",
    "Non-stationary ETAS"
  ),
  
  logLik = c(
    logLik_stationary_custom,
    logLik_nonstationary
  ),
  
  k = c(
    k_stationary_custom,
    k_nonstationary
  )
)

etas_model_comparison <-
  etas_model_comparison %>%
  mutate(
    AIC =
      -2 * logLik +
      2 * k,
    
    BIC =
      -2 * logLik +
      log(N) * k,
    
    delta_AIC =
      AIC - min(AIC),
    
    delta_BIC =
      BIC - min(BIC)
  )

print(etas_model_comparison)

write.csv(
  etas_model_comparison,
  "stationary_vs_nonstationary_ETAS_Mc_2.csv",
  row.names = FALSE
)

background_rate_results <-
  data.frame(
    segment = 1:J,
    
    start_day =
      segment_breaks[
        1:J
      ],
    
    end_day =
      pmin(
        segment_breaks[
          2:(J + 1)
        ],
        T_end
      ),
    
    mu =
      mu_nonstationary
  )

background_rate_results <-
  background_rate_results %>%
  mutate(
    start_date =
      as.POSIXct(
        t0 +
          start_day *
          24 * 3600,
        origin = "1970-01-01"
      ),
    
    end_date =
      as.POSIXct(
        t0 +
          end_day *
          24 * 3600,
        origin = "1970-01-01"
      )
  )

print(background_rate_results)

# =0说明收敛正常
fit_stationary$convergence
fit_nonstationary$convergence


# ============================================================
# Spatial distribution map
# Point size = magnitude
# Point color = 2019 vs non-2019
# Main catalogue: Mc = 2.0
# ============================================================

# Add period group
catalog_main_spatial <- catalog_main %>%
  mutate(
    period_group = ifelse(
      year(datetime) == 2019,
      "2019",
      "non-2019"
    )
  )

# California map
usa_map <- map_data("state")

california_map <- usa_map %>%
  filter(region == "california")

# Plot
p_spatial_map <- ggplot() +
  
  geom_polygon(
    data = california_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "gray95",
    color = "gray60"
  ) +
  
  geom_point(
    data = catalog_main_spatial,
    aes(
      x = lon,
      y = lat,
      size = mag,
      color = period_group
    ),
    alpha = 0.45
  ) +
  
  coord_fixed(
    xlim = c(
      min(catalog_main_spatial$lon, na.rm = TRUE) - 0.2,
      max(catalog_main_spatial$lon, na.rm = TRUE) + 0.2
    ),
    ylim = c(
      min(catalog_main_spatial$lat, na.rm = TRUE) - 0.2,
      max(catalog_main_spatial$lat, na.rm = TRUE) + 0.2
    )
  ) +
  
  scale_size_continuous(
    range = c(0.5, 4)
  ) +
  
  labs(
    title = "Spatial Distribution of Earthquakes",
    subtitle = "Point size represents magnitude; color distinguishes 2019 from other years; M >= 2.0",
    x = "Longitude",
    y = "Latitude",
    size = "Magnitude",
    color = "Period"
  ) +
  
  theme_minimal()

print(p_spatial_map)

ggsave(
  "spatial_distribution_2019_vs_non2019_Mc_2.png",
  p_spatial_map,
  width = 8,
  height = 6,
  dpi = 300
)
catalog_non2019 <- catalog_main_spatial %>%
  filter(period_group == "non-2019")

catalog_2019 <- catalog_main_spatial %>%
  filter(period_group == "2019")

p_spatial_map <- ggplot() +
  
  geom_polygon(
    data = california_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "gray95",
    color = "gray60"
  ) +
  
  # Plot non-2019 first
  geom_point(
    data = catalog_non2019,
    aes(
      x = lon,
      y = lat,
      size = mag,
      color = period_group
    ),
    alpha = 0.25
  ) +
  
  # Plot 2019 on top
  geom_point(
    data = catalog_2019,
    aes(
      x = lon,
      y = lat,
      size = mag,
      color = period_group
    ),
    alpha = 0.65
  ) +
  
  coord_fixed(
    xlim = c(
      min(catalog_main_spatial$lon, na.rm = TRUE) - 0.2,
      max(catalog_main_spatial$lon, na.rm = TRUE) + 0.2
    ),
    ylim = c(
      min(catalog_main_spatial$lat, na.rm = TRUE) - 0.2,
      max(catalog_main_spatial$lat, na.rm = TRUE) + 0.2
    )
  ) +
  
  scale_size_continuous(
    range = c(0.5, 4)
  ) +
  
  labs(
    title = "Spatial Distribution of Earthquakes",
    subtitle = "Magnitude shown by point size; 2019 events highlighted; M >= 2.0",
    x = "Longitude",
    y = "Latitude",
    size = "Magnitude",
    color = "Period"
  ) +
  
  theme_minimal()

# Display plot
print(p_spatial_map)


# Save plot
ggsave(
  "spatial_distribution_2019_vs_non2019_Mc_2.png",
  p_spatial_map,
  width = 8,
  height = 6,
  dpi = 300
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
  "seismicity_around_largest_event_Mc_2.png",
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
  "daily_counts_around_largest_event_Mc_2.png",
  p_around_daily,
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
  "descriptive_analysis_text_summary_Mc_2.txt"
)
# ============================================================
# STEP 3: Declustering and ETAS initialization
# Main catalogue: Mc = 2.0
# ============================================================
# Methods
# A. Gardner-Knopoff window declustering
# B. Zaliapin-Ben-Zion-style nearest-neighbour cluster declustering
# C. Reasenberg-style interaction-link declustering
#
# Experimental principle:
# SAME input catalogue -> different declustering schemes -> harmonised
# background-like indicator -> comparable ETAS initialisation.
# ============================================================

`%||%` <- function(x, y) if (is.null(x) || length(x)==0 || all(is.na(x))) y else x

# ============================================================
# 0. LOCK THE COMMON INPUT CATALOGUE
# ============================================================
mc_main <- 2.0
bin_width <- 0.1
mc_tag <- gsub("\\.", "_", sprintf("%.1f", mc_main))

catalog_main <- catalog %>%
  filter(!is.na(datetime), !is.na(mag), !is.na(lat), !is.na(lon), mag >= mc_main) %>%
  arrange(datetime) %>%
  mutate(event_index = row_number(), original_order = row_number(), year = year(datetime))

stopifnot(all(catalog_main$mag >= mc_main))
cat("Main declustering catalogue: Mc =", mc_main, "; N =", nrow(catalog_main), "\n")

# ============================================================
# 1. SHARED HELPERS
# ============================================================
haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371.227
  lon1 <- lon1*pi/180; lat1 <- lat1*pi/180
  lon2 <- lon2*pi/180; lat2 <- lat2*pi/180
  dlon <- lon2-lon1; dlat <- lat2-lat1
  a <- sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2
  R * 2 * atan2(sqrt(a), sqrt(1-a))
}

decimal_year <- function(datetime) {
  tz <- attr(datetime, "tzone") %||% "UTC"
  y <- year(datetime)
  y0 <- as.POSIXct(paste0(y,"-01-01 00:00:00"), tz=tz)
  y1 <- as.POSIXct(paste0(y+1,"-01-01 00:00:00"), tz=tz)
  y + as.numeric(difftime(datetime,y0,units="secs")) /
    as.numeric(difftime(y1,y0,units="secs"))
}

estimate_b_value_local <- function(mags, mc, bin_width=0.1) {
  x <- mags[is.finite(mags) & mags >= mc]
  n <- length(x)
  if (n < 30) return(data.frame(mc=mc,n=n,b=NA,delta_b=NA))
  mean_mag <- mean(x)
  b <- log10(exp(1))/(mean_mag-(mc-bin_width/2))
  delta_b <- 2.3*b^2*sqrt(sum((x-mean_mag)^2)/(n*(n-1)))
  data.frame(mc=mc,n=n,b=b,delta_b=delta_b)
}

b_info_declustering <- estimate_b_value_local(catalog_main$mag, mc_main, bin_width)
b_catalogue <- b_info_declustering$b
cat("b-value estimated at Mc=2.0:", b_catalogue, "\n")

save_csv_mc <- function(x, stem) {
  f <- paste0(stem,"_Mc_",mc_tag,".csv")
  write.csv(x,f,row.names=FALSE); invisible(f)
}
save_plot_mc <- function(p, stem, width=9, height=5, dpi=300) {
  f <- paste0(stem,"_Mc_",mc_tag,".png")
  ggsave(f,p,width=width,height=height,dpi=dpi); invisible(f)
}

# ============================================================
# 2. A. GARDNER-KNOPOFF
# Literature/reference implementation style: OpenQuake/HMTK Type 1
# Baseline fs_time_prop = 1.0; 0.5 retained for sensitivity.
# ============================================================
gk_space_window_km <- function(M) 10^(0.1238*M + 0.983)
gk_time_window_days <- function(M) ifelse(M < 6.5,
                                          10^(0.5409*M - 0.547), 10^(0.032*M + 2.7389))

gardner_knopoff_declustering <- function(df, fs_time_prop=1.0, time_cutoff_days=NULL) {
  stopifnot(fs_time_prop>=0, fs_time_prop<=1)
  d0 <- df %>% arrange(datetime) %>%
    mutate(original_order=row_number(),
           gk_space_window_km=gk_space_window_km(mag),
           gk_time_window_days=gk_time_window_days(mag))
  if (!is.null(time_cutoff_days))
    d0$gk_time_window_days <- pmin(d0$gk_time_window_days,time_cutoff_days)
  
  d <- d0 %>% arrange(desc(mag),datetime)
  n <- nrow(d)
  d$gk_cluster_id <- 0L; d$gk_flag <- 0L
  d$gk_parent_evid <- NA; d$gk_parent_mag <- NA_real_
  cid <- 0L
  
  if (n>=2) for (i in seq_len(n-1)) {
    if (i%%1000==0) cat("GK",i,"/",n,"\n")
    if (d$gk_cluster_id[i]!=0L) next
    dt <- as.numeric(difftime(d$datetime,d$datetime[i],units="days"))
    keep_t <- d$gk_cluster_id==0L &
      dt >= -d$gk_time_window_days[i]*fs_time_prop &
      dt <= d$gk_time_window_days[i]
    idx <- which(keep_t); if (length(idx)<=1) next
    dist <- haversine_km(d$lon[idx],d$lat[idx],d$lon[i],d$lat[i])
    inside <- idx[dist <= d$gk_space_window_km[i]]
    if (length(setdiff(inside,i))==0) next
    cid <- cid+1L
    d$gk_cluster_id[inside] <- cid
    d$gk_flag[inside] <- 1L
    d$gk_flag[inside[dt[inside]<0]] <- -1L
    d$gk_flag[i] <- 0L
    nonmain <- setdiff(inside,i)
    d$gk_parent_evid[nonmain] <- d$evid[i]
    d$gk_parent_mag[nonmain] <- d$mag[i]
  }
  
  d %>% arrange(original_order) %>% mutate(
    gk_label=case_when(
      gk_cluster_id==0L ~ "background",
      gk_cluster_id!=0L & gk_flag==0L ~ "mainshock",
      gk_flag==1L ~ "aftershock",
      gk_flag==-1L ~ "foreshock",
      TRUE ~ "unknown"),
    gk_is_background_like=gk_label %in% c("background","mainshock"),
    gk_background_prob=as.numeric(gk_is_background_like))
}

gk_fs_baseline <- 1.0
declust_gk <- gardner_knopoff_declustering(catalog_main,gk_fs_baseline)
gk_summary <- declust_gk %>% count(gk_label,name="n") %>%
  mutate(percentage=100*n/sum(n),method="Gardner-Knopoff")
print(gk_summary)
save_csv_mc(declust_gk,"declustering_gardner_knopoff")
save_csv_mc(gk_summary,"declustering_gardner_knopoff_summary")

etas_init_gk <- declust_gk %>% select(
  evid,datetime,lat,lon,depth,mag,gk_cluster_id,gk_flag,gk_label,
  gk_is_background_like,gk_background_prob,gk_parent_evid,gk_parent_mag,
  gk_space_window_km,gk_time_window_days)
save_csv_mc(etas_init_gk,"etas_initialization_gardner_knopoff")

# ============================================================
# 3. B. ZALIAPIN-BEN-ZION-STYLE NEAREST NEIGHBOUR
# Literature-aligned deterministic cluster-analysis version:
# - eta combines time, distance, parent magnitude
# - d_f baseline = 1.6 for epicentral analysis
# - b estimated from this Mc=2.0 catalogue
# - 2-Gaussian mixture replaces the old k-means threshold
# NOTE: This is NOT the exact stochastic random-thinning algorithm of ZB2020.
# ============================================================
compute_nearest_neighbor_eta <- function(df,mc,b_value,d_f=1.6,use_depth=FALSE,
                                         max_lookback_years=Inf,max_previous_events=Inf) {
  d <- df %>% arrange(datetime) %>% mutate(
    event_index=row_number(),decimal_year=decimal_year(datetime),
    nn_parent_index=NA_integer_,nn_parent_evid=NA,
    nn_time_years=NA_real_,nn_distance_km=NA_real_,
    nn_parent_mag=NA_real_,nn_log_eta=NA_real_)
  n <- nrow(d)
  if (n>=2) for (j in 2:n) {
    if (j%%1000==0) cat("NN",j,"/",n,"\n")
    prev <- seq_len(j-1)
    dtall <- d$decimal_year[j]-d$decimal_year[prev]
    keep <- dtall>0
    if (is.finite(max_lookback_years)) keep <- keep & dtall<=max_lookback_years
    prev <- prev[keep]; if (!length(prev)) next
    if (is.finite(max_previous_events) && length(prev)>max_previous_events)
      prev <- tail(prev,as.integer(max_previous_events))
    dt <- d$decimal_year[j]-d$decimal_year[prev]
    dt[dt<=0] <- 1/(365.25*24*3600)
    epi <- haversine_km(d$lon[prev],d$lat[prev],d$lon[j],d$lat[j])
    epi[epi<=0] <- 0.001
    r <- if (use_depth) sqrt(epi^2+(d$depth[j]-d$depth[prev])^2) else epi
    r[r<=0] <- 0.001
    pm <- d$mag[prev]
    logeta <- log10(dt)+d_f*log10(r)-b_value*(pm-mc)
    q <- which.min(logeta); parent <- prev[q]
    d$nn_parent_index[j] <- parent; d$nn_parent_evid[j] <- d$evid[parent]
    d$nn_time_years[j] <- dt[q]; d$nn_distance_km[j] <- r[q]
    d$nn_parent_mag[j] <- pm[q]; d$nn_log_eta[j] <- logeta[q]
  }
  d
}

fit_nn_mixture_threshold <- function(log_eta,seed=123) {
  x <- log_eta[is.finite(log_eta)]
  if (length(x)<100) stop("Too few finite eta values")
  set.seed(seed)
  fit <- mclust::Mclust(x,G=2,modelNames=c("E","V"),verbose=FALSE)
  grid <- seq(quantile(x,.001),quantile(x,.999),length.out=4000)
  pr <- predict(fit,newdata=grid)
  thr <- if (!is.null(pr$z) && ncol(pr$z)==2)
    grid[which.min(abs(pr$z[,1]-pr$z[,2]))] else mean(fit$parameters$mean)
  list(fit=fit,threshold=thr,means=sort(as.numeric(fit$parameters$mean)))
}

classify_nn_clusters <- function(df,threshold) {
  d <- df %>% arrange(datetime) %>% mutate(
    nn_is_clustered_edge=ifelse(is.na(nn_log_eta),FALSE,nn_log_eta<threshold),
    nn_cluster_id=0L,nn_mainshock_evid=NA,nn_label="background")
  n <- nrow(d); adj <- vector("list",n); for(i in seq_len(n)) adj[[i]] <- integer(0)
  for(j in seq_len(n)) {
    p <- d$nn_parent_index[j]
    if (!is.na(p) && d$nn_is_clustered_edge[j]) {
      adj[[j]] <- unique(c(adj[[j]],p)); adj[[p]] <- unique(c(adj[[p]],j))
    }
  }
  visited <- rep(FALSE,n); cid <- 0L
  for(i in seq_len(n)) {
    if (visited[i]) next
    stack <- i; comp <- integer(0)
    while(length(stack)) {
      v <- stack[1]; stack <- stack[-1]
      if (visited[v]) next
      visited[v] <- TRUE; comp <- c(comp,v); stack <- unique(c(stack,adj[[v]]))
    }
    if (length(comp)==1) next
    cid <- cid+1L
    m <- comp[which.max(d$mag[comp])]; mt <- d$datetime[m]; me <- d$evid[m]
    d$nn_cluster_id[comp] <- cid; d$nn_mainshock_evid[comp] <- me
    for(k in comp) d$nn_label[k] <- if(k==m) "mainshock" else if(d$datetime[k]<mt) "foreshock" else "aftershock"
  }
  d %>% mutate(nn_is_background_like=nn_label %in% c("background","mainshock"),
               nn_background_prob=as.numeric(nn_is_background_like))
}

nn_df_baseline <- 1.6
nn_raw <- compute_nearest_neighbor_eta(catalog_main,mc_main,b_catalogue,nn_df_baseline,
                                       use_depth=FALSE,max_lookback_years=Inf,max_previous_events=Inf)
nn_mix <- fit_nn_mixture_threshold(nn_raw$nn_log_eta)
nn_eta_threshold <- nn_mix$threshold
cat("NN: d_f=",nn_df_baseline," b=",b_catalogue," threshold=",nn_eta_threshold,"\n")
declust_nn <- classify_nn_clusters(nn_raw,nn_eta_threshold)
nn_summary <- declust_nn %>% count(nn_label,name="n") %>%
  mutate(percentage=100*n/sum(n),method="Nearest-neighbour")
print(nn_summary)
save_csv_mc(declust_nn,"declustering_nearest_neighbor")
save_csv_mc(nn_summary,"declustering_nearest_neighbor_summary")

etas_init_nn <- declust_nn %>% select(
  evid,datetime,lat,lon,depth,mag,nn_parent_evid,nn_parent_index,
  nn_time_years,nn_distance_km,nn_parent_mag,nn_log_eta,nn_cluster_id,
  nn_mainshock_evid,nn_label,nn_is_background_like,nn_background_prob)
save_csv_mc(etas_init_nn,"etas_initialization_nearest_neighbor")

p_nn_eta <- ggplot(declust_nn %>% filter(is.finite(nn_log_eta)),aes(nn_log_eta))+
  geom_histogram(bins=80)+geom_vline(xintercept=nn_eta_threshold,linetype="dashed")+
  labs(title="Nearest-neighbour proximity distribution",
       subtitle="Two-component Gaussian-mixture threshold",
       x=expression(log[10](eta)),y="Number of events")+theme_minimal()
print(p_nn_eta); save_plot_mc(p_nn_eta,"nearest_neighbor_log_eta_distribution",8,5)

# ============================================================
# 4. C. REASENBERG-STYLE INTERACTION-LINK DECLUSTERING
# Standard literature baseline values:
# tau_min=1 day, tau_max=10 days, p1=.95, xk=.5, xmeff=4, rfact=10
# IMPORTANT: still label as "Reasenberg-style" unless the final thesis uses
# original CLUSTER2000/ZMAP output.
# ============================================================
reasenberg_crack_radius_km <- function(M) 10^(0.4*M-1.2)
reasenberg_interaction_radius_km <- function(M,rfact=10) rfact*reasenberg_crack_radius_km(M)

# Literature-inspired dynamic look-ahead, bounded by tau_min/tau_max.
# This makes p1 active; it is not a byte-for-byte CLUSTER2000 port.
reasenberg_tau <- function(cluster_max_mag,elapsed_days,tau_min=1,tau_max=10,
                           p1=.95,xmeff=4,xk=.5) {
  effective_cutoff <- xmeff + xk*pmax(cluster_max_mag-xmeff,0)
  delta_m <- pmax(cluster_max_mag-effective_cutoff,0)
  base <- pmax(elapsed_days,tau_min)
  tau <- -log(1-p1)*base/10^(2*(delta_m-1)/3)
  pmin(tau_max,pmax(tau_min,tau))
}

reasenberg_style_declustering <- function(df,tau_min=1,tau_max=10,p1=.95,
                                          xk=.5,xmeff=4,rfact=10) {
  d <- df %>% arrange(datetime) %>% mutate(
    reasenberg_cluster_id=0L,reasenberg_parent_index=NA_integer_,
    reasenberg_parent_evid=NA,reasenberg_distance_km=NA_real_,
    reasenberg_time_days=NA_real_)
  n <- nrow(d); uf <- seq_len(n)
  root <- function(x) { while(uf[x]!=x){uf[x]<<-uf[uf[x]];x<-uf[x]};x }
  unite <- function(a,b){ra<-root(a);rb<-root(b);if(ra!=rb)uf[rb]<<-ra}
  
  if(n>=2) for(j in 2:n) {
    if(j%%1000==0) cat("Reasenberg-style",j,"/",n,"\n")
    dtall <- as.numeric(difftime(d$datetime[j],d$datetime[seq_len(j-1)],units="days"))
    cand <- which(dtall>0 & dtall<=tau_max); if(!length(cand)) next
    linked <- integer(0); ldist <- numeric(0)
    for(i in cand) {
      ri <- root(i); members <- which(sapply(seq_len(j-1),root)==ri)
      Mmax <- max(d$mag[members]); tstart <- min(d$datetime[members])
      elapsed <- as.numeric(difftime(d$datetime[i],tstart,units="days"))
      tau_i <- reasenberg_tau(Mmax,elapsed,tau_min,tau_max,p1,xmeff,xk)
      dtij <- as.numeric(difftime(d$datetime[j],d$datetime[i],units="days"))
      if(dtij>tau_i) next
      dij <- haversine_km(d$lon[i],d$lat[i],d$lon[j],d$lat[j])
      if(dij<=reasenberg_interaction_radius_km(d$mag[i],rfact)) {
        linked <- c(linked,i); ldist <- c(ldist,dij)
      }
    }
    if(!length(linked)) next
    for(i in linked) unite(i,j)
    q <- which.min(ldist); p <- linked[q]
    d$reasenberg_parent_index[j] <- p; d$reasenberg_parent_evid[j] <- d$evid[p]
    d$reasenberg_distance_km[j] <- ldist[q]
    d$reasenberg_time_days[j] <- as.numeric(difftime(d$datetime[j],d$datetime[p],units="days"))
  }
  
  roots <- sapply(seq_len(n),root); tab <- table(roots)
  clustered_roots <- as.integer(names(tab[tab>1])); cmap <- setNames(seq_along(clustered_roots),clustered_roots)
  d$reasenberg_cluster_id <- ifelse(roots %in% clustered_roots,
                                    unname(cmap[as.character(roots)]),0L)
  d$reasenberg_label <- "background"; d$reasenberg_mainshock_evid <- NA
  for(cid in sort(unique(d$reasenberg_cluster_id[d$reasenberg_cluster_id>0]))) {
    idx <- which(d$reasenberg_cluster_id==cid); m <- idx[which.max(d$mag[idx])]
    mt <- d$datetime[m]; me <- d$evid[m]; d$reasenberg_mainshock_evid[idx] <- me
    d$reasenberg_label[idx] <- ifelse(idx==m,"mainshock",
                                      ifelse(d$datetime[idx]<mt,"foreshock","aftershock"))
  }
  d %>% mutate(reasenberg_is_background_like=reasenberg_label %in% c("background","mainshock"),
               reasenberg_background_prob=as.numeric(reasenberg_is_background_like))
}

reas_baseline <- list(tau_min=1,tau_max=10,p1=.95,xk=.5,xmeff=4,rfact=10)
declust_reasenberg <- reasenberg_style_declustering(
  catalog_main,reas_baseline$tau_min,reas_baseline$tau_max,reas_baseline$p1,
  reas_baseline$xk,reas_baseline$xmeff,reas_baseline$rfact)
reasenberg_summary <- declust_reasenberg %>% count(reasenberg_label,name="n") %>%
  mutate(percentage=100*n/sum(n),method="Reasenberg-style")
print(reasenberg_summary)
save_csv_mc(declust_reasenberg,"declustering_reasenberg_style")
save_csv_mc(reasenberg_summary,"declustering_reasenberg_style_summary")

etas_init_reasenberg <- declust_reasenberg %>% select(
  evid,datetime,lat,lon,depth,mag,reasenberg_cluster_id,reasenberg_parent_index,
  reasenberg_parent_evid,reasenberg_distance_km,reasenberg_time_days,
  reasenberg_mainshock_evid,reasenberg_label,reasenberg_is_background_like,
  reasenberg_background_prob)
save_csv_mc(etas_init_reasenberg,"etas_initialization_reasenberg_style")

# ============================================================
# 5. HARMONISE EVENT-LEVEL OUTPUTS
# ============================================================
event_level_compare <- catalog_main %>% select(evid,datetime,mag,lat,lon,year) %>%
  left_join(declust_gk %>% select(evid,gk_label,gk_is_background_like),by="evid") %>%
  left_join(declust_nn %>% select(evid,nn_label,nn_is_background_like),by="evid") %>%
  left_join(declust_reasenberg %>% select(evid,reasenberg_label,reasenberg_is_background_like),by="evid")
stopifnot(nrow(event_level_compare)==nrow(catalog_main))
stopifnot(!anyNA(event_level_compare[c("gk_is_background_like","nn_is_background_like","reasenberg_is_background_like")]))
save_csv_mc(event_level_compare,"declustering_event_level_comparison")

# ============================================================
# 6. DIAGNOSTIC 1: OVERALL BACKGROUND FRACTION
# ============================================================
overall_background_summary <- bind_rows(
  data.frame(method="Gardner-Knopoff",n_total=nrow(event_level_compare),
             n_background_like=sum(event_level_compare$gk_is_background_like)),
  data.frame(method="Nearest-neighbour",n_total=nrow(event_level_compare),
             n_background_like=sum(event_level_compare$nn_is_background_like)),
  data.frame(method="Reasenberg-style",n_total=nrow(event_level_compare),
             n_background_like=sum(event_level_compare$reasenberg_is_background_like))) %>%
  mutate(n_clustered=n_total-n_background_like,
         background_fraction=n_background_like/n_total,
         background_percentage=100*background_fraction)
print(overall_background_summary)
save_csv_mc(overall_background_summary,"diagnostic_overall_background_fraction")

# ============================================================
# 7. DIAGNOSTIC 2: EVENT-LEVEL AGREEMENT
# ============================================================
binary_jaccard <- function(x,y) { u<-sum(x|y); if(u==0) NA_real_ else sum(x&y)/u }
cohen_kappa_binary <- function(x,y) {
  x<-as.integer(x); y<-as.integer(y); po<-mean(x==y)
  pe<-mean(x==1)*mean(y==1)+mean(x==0)*mean(y==0)
  if(abs(1-pe)<1e-12) NA_real_ else (po-pe)/(1-pe)
}
agreement_pair <- function(x,y,m1,m2) data.frame(
  method_1=m1,method_2=m2,raw_agreement=mean(x==y),
  jaccard_background=binary_jaccard(x,y),kappa=cohen_kappa_binary(x,y),
  n_disagree=sum(x!=y),disagreement_percentage=100*mean(x!=y))

pairwise_agreement <- bind_rows(
  agreement_pair(event_level_compare$gk_is_background_like,event_level_compare$nn_is_background_like,
                 "Gardner-Knopoff","Nearest-neighbour"),
  agreement_pair(event_level_compare$gk_is_background_like,event_level_compare$reasenberg_is_background_like,
                 "Gardner-Knopoff","Reasenberg-style"),
  agreement_pair(event_level_compare$nn_is_background_like,event_level_compare$reasenberg_is_background_like,
                 "Nearest-neighbour","Reasenberg-style"))
print(pairwise_agreement)
save_csv_mc(pairwise_agreement,"diagnostic_pairwise_event_agreement")

event_level_disagreement <- event_level_compare %>% mutate(
  n_background_votes=as.integer(gk_is_background_like)+as.integer(nn_is_background_like)+as.integer(reasenberg_is_background_like),
  unanimous=n_background_votes %in% c(0,3),any_disagreement=!unanimous)
disagreement_summary <- event_level_disagreement %>% summarise(
  n_total=n(),n_unanimous=sum(unanimous),n_disagreement=sum(any_disagreement),
  disagreement_percentage=100*mean(any_disagreement))
print(disagreement_summary)
save_csv_mc(event_level_disagreement,"diagnostic_event_level_disagreement")
save_csv_mc(disagreement_summary,"diagnostic_event_level_disagreement_summary")

# ============================================================
# 8. DIAGNOSTIC 3: 2019 VS NON-2019
# ============================================================
period_summary <- event_level_compare %>% mutate(period_group=ifelse(year==2019,"2019","non-2019")) %>%
  group_by(period_group) %>% summarise(
    n_total=n(),
    gk_background=sum(gk_is_background_like),gk_background_fraction=mean(gk_is_background_like),
    nn_background=sum(nn_is_background_like),nn_background_fraction=mean(nn_is_background_like),
    reasenberg_background=sum(reasenberg_is_background_like),reasenberg_background_fraction=mean(reasenberg_is_background_like),
    all_three_background=sum(gk_is_background_like&nn_is_background_like&reasenberg_is_background_like),
    all_three_clustered=sum(!gk_is_background_like&!nn_is_background_like&!reasenberg_is_background_like),
    any_method_disagreement=sum((as.integer(gk_is_background_like)+as.integer(nn_is_background_like)+as.integer(reasenberg_is_background_like)) %in% c(1,2)),
    disagreement_fraction=any_method_disagreement/n_total,.groups="drop")
print(period_summary)
save_csv_mc(period_summary,"diagnostic_2019_vs_non2019_declustering")

agreement_by_period <- event_level_compare %>% mutate(period_group=ifelse(year==2019,"2019","non-2019")) %>%
  group_split(period_group) %>% map_dfr(function(z) bind_rows(
    agreement_pair(z$gk_is_background_like,z$nn_is_background_like,"Gardner-Knopoff","Nearest-neighbour"),
    agreement_pair(z$gk_is_background_like,z$reasenberg_is_background_like,"Gardner-Knopoff","Reasenberg-style"),
    agreement_pair(z$nn_is_background_like,z$reasenberg_is_background_like,"Nearest-neighbour","Reasenberg-style")) %>%
      mutate(period_group=unique(z$period_group)))
print(agreement_by_period)
save_csv_mc(agreement_by_period,"diagnostic_pairwise_agreement_2019_vs_non2019")

# ============================================================
# 9. MONTHLY ETAS INITIAL BACKGROUND RATE (events/day)
# ============================================================
monthly_background_rate <- function(df,background_col,method_name) {
  tmp <- df %>% mutate(month=floor_date(datetime,"month"),bg=as.numeric(.data[[background_col]])) %>%
    group_by(month) %>% summarise(total_events=n(),background_events=sum(bg),.groups="drop")
  full <- tibble(month=seq(floor_date(min(catalog_main$datetime),"month"),
                           floor_date(max(catalog_main$datetime),"month"),by="month"))
  full %>% left_join(tmp,by="month") %>% mutate(
    total_events=replace_na(total_events,0),background_events=replace_na(background_events,0),
    days_in_month=days_in_month(month),
    background_rate_per_day=background_events/days_in_month,
    background_count_per_month=background_events,method=method_name)
}

bg_rate_gk_monthly <- monthly_background_rate(declust_gk,"gk_is_background_like","Gardner-Knopoff")
bg_rate_nn_monthly <- monthly_background_rate(declust_nn,"nn_is_background_like","Nearest-neighbour")
bg_rate_reasenberg_monthly <- monthly_background_rate(declust_reasenberg,"reasenberg_is_background_like","Reasenberg-style")
bg_rate_three_methods <- bind_rows(bg_rate_gk_monthly,bg_rate_nn_monthly,bg_rate_reasenberg_monthly)
save_csv_mc(bg_rate_three_methods,"initial_background_rate_three_methods")

p_bg_three <- ggplot(bg_rate_three_methods,aes(month,background_rate_per_day,linetype=method))+
  geom_line(linewidth=.75)+labs(title="Comparison of Declustering-based Initial Background Rates",
                                subtitle=paste0("Common input catalogue: M >= ",mc_main),x="Time",y="Background-like events per day",
                                linetype="Declustering method")+theme_minimal()
print(p_bg_three); save_plot_mc(p_bg_three,"initial_background_rate_three_methods",10,5)

# ============================================================
# 10. DIAGNOSTIC 4: QUANTIFY INITIALIZATION DIFFERENCES
# ============================================================
bgwide <- bg_rate_three_methods %>% select(month,method,background_rate_per_day) %>%
  pivot_wider(names_from=method,values_from=background_rate_per_day)
init_metrics <- function(x,y,n1,n2) data.frame(method_1=n1,method_2=n2,
                                               MAE=mean(abs(x-y),na.rm=TRUE),RMSE=sqrt(mean((x-y)^2,na.rm=TRUE)),
                                               correlation=cor(x,y,use="complete.obs"))
init_difference_summary <- bind_rows(
  init_metrics(bgwide[["Gardner-Knopoff"]],bgwide[["Nearest-neighbour"]],"Gardner-Knopoff","Nearest-neighbour"),
  init_metrics(bgwide[["Gardner-Knopoff"]],bgwide[["Reasenberg-style"]],"Gardner-Knopoff","Reasenberg-style"),
  init_metrics(bgwide[["Nearest-neighbour"]],bgwide[["Reasenberg-style"]],"Nearest-neighbour","Reasenberg-style"))
print(init_difference_summary)
save_csv_mc(init_difference_summary,"diagnostic_initialisation_difference_metrics")

bg2019 <- bg_rate_three_methods %>% filter(year(month)==2019) %>%
  select(month,method,background_rate_per_day) %>% pivot_wider(names_from=method,values_from=background_rate_per_day)
init_difference_2019 <- bind_rows(
  init_metrics(bg2019[["Gardner-Knopoff"]],bg2019[["Nearest-neighbour"]],"Gardner-Knopoff","Nearest-neighbour"),
  init_metrics(bg2019[["Gardner-Knopoff"]],bg2019[["Reasenberg-style"]],"Gardner-Knopoff","Reasenberg-style"),
  init_metrics(bg2019[["Nearest-neighbour"]],bg2019[["Reasenberg-style"]],"Nearest-neighbour","Reasenberg-style")) %>%
  mutate(period="2019")
print(init_difference_2019)
save_csv_mc(init_difference_2019,"diagnostic_initialisation_difference_metrics_2019")

# ============================================================
# 11. OPTIONAL MINIMAL PARAMETER SENSITIVITY
# Do NOT tune for a desired result; this is a robustness analysis.
# ============================================================
RUN_SENSITIVITY <- FALSE
if (RUN_SENSITIVITY) {
  gk_sens <- map_dfr(c(.5,1.0), function(fs){
    z<-gardner_knopoff_declustering(catalog_main,fs)
    data.frame(method="Gardner-Knopoff",parameter="fs_time_prop",value=fs,
               background_fraction=mean(z$gk_is_background_like))})
  
  nn_sens <- map_dfr(c(1.4,1.6,1.8), function(dfv){
    r<-compute_nearest_neighbor_eta(catalog_main,mc_main,b_catalogue,dfv)
    mix<-fit_nn_mixture_threshold(r$nn_log_eta); z<-classify_nn_clusters(r,mix$threshold)
    data.frame(method="Nearest-neighbour",parameter="d_f",value=dfv,threshold=mix$threshold,
               background_fraction=mean(z$nn_is_background_like))})
  
  rg <- expand.grid(rfact=c(5,10,20),tau_max=c(5,10,15))
  reas_sens <- pmap_dfr(rg,function(rfact,tau_max){
    z<-reasenberg_style_declustering(catalog_main,1,tau_max,.95,.5,1.5,rfact)
    data.frame(method="Reasenberg-style",rfact=rfact,tau_max=tau_max,
               background_fraction=mean(z$reasenberg_is_background_like))})
  
  save_csv_mc(gk_sens,"sensitivity_gardner_knopoff")
  save_csv_mc(nn_sens,"sensitivity_nearest_neighbor")
  save_csv_mc(reas_sens,"sensitivity_reasenberg_style")
}

# ============================================================
# 12. FINAL SANITY CHECKS
# ============================================================
stopifnot(nrow(catalog_main)==nrow(declust_gk),
          nrow(catalog_main)==nrow(declust_nn),
          nrow(catalog_main)==nrow(declust_reasenberg))
cat("\nAll three methods used the same Mc=",mc_main," catalogue; N=",nrow(catalog_main),"\n",sep="")
cat("\nOverall background fractions:\n"); print(overall_background_summary)
cat("\nPairwise agreement:\n"); print(pairwise_agreement)
cat("\n2019 diagnostic:\n"); print(period_summary)
cat("\nInitialization differences:\n"); print(init_difference_summary)

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

mc_main <- 2

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
