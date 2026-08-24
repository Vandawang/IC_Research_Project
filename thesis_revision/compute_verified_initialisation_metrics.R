#!/usr/bin/env Rscript

# Reconstruct the actual spline background functions used as ETAS starting
# values and compare them with the final fitted functions on the same daily
# grid. This corrects the earlier propagation ratios, whose numerator used a
# daily grid while the denominator used unweighted calendar-month rates.

bg <- read.csv(
  "initial_background_rate_three_methods_Mc_2_0.csv",
  check.names = FALSE
)
final <- read.csv(
  "nonstationary_ETAS_final_background_rates_Mc_2_0.csv",
  check.names = FALSE
)
catalogue <- read.csv("catalog_for_etas_Mc_2.csv", check.names = FALSE)

event_datetime <- as.POSIXct(catalogue$datetime, tz = "UTC")
t_origin <- min(event_datetime)
t_end_datetime <- max(event_datetime)
t_event <- as.numeric(difftime(event_datetime, t_origin, units = "days"))
t_end <- max(t_event)

B_event <- splines::ns(
  t_event,
  df = 8,
  intercept = TRUE,
  Boundary.knots = c(0, t_end)
)

t_grid <- sort(unique(final$t_days))
B_grid <- predict(B_event, newx = t_grid)

make_initial_mu <- function(method_name) {
  method_data <- bg[bg$method == method_name, ]
  method_data <- method_data[order(method_data$month), ]

  month_start <- as.POSIXct(method_data$month, tz = "UTC")
  month_sequence <- seq(
    as.Date(format(min(month_start), "%Y-%m-01")),
    by = "month",
    length.out = nrow(method_data) + 1
  )
  month_end <- as.POSIXct(month_sequence[-1], tz = "UTC")

  month_start_clipped <- pmax(month_start, t_origin)
  month_end_clipped <- pmin(month_end, t_end_datetime)
  exposure_days <- as.numeric(difftime(
    month_end_clipped,
    month_start_clipped,
    units = "days"
  ))
  month_midpoint <- month_start_clipped +
    (month_end_clipped - month_start_clipped) / 2
  midpoint_days <- as.numeric(difftime(
    month_midpoint,
    t_origin,
    units = "days"
  ))

  log_initial_rate <- log(
    (method_data$background_events + 0.1) / pmax(exposure_days, 1)
  )
  B_midpoint <- predict(B_event, newx = midpoint_days)
  beta_initial <- solve(
    crossprod(B_midpoint) + 1e-4 * diag(ncol(B_midpoint)),
    crossprod(B_midpoint, log_initial_rate)
  )

  as.numeric(exp(B_grid %*% beta_initial))
}

methods <- c(
  "Gardner-Knopoff",
  "Nearest-neighbour",
  "Reasenberg-style"
)
initial_mu <- setNames(lapply(methods, make_initial_mu), methods)

final_mu <- lapply(methods, function(method_name) {
  z <- final[final$method == method_name, ]
  z <- z[order(z$t_days), ]
  z$mu
})
names(final_mu) <- methods

metric_row <- function(x, y, method_1, method_2) {
  data.frame(
    method_1 = method_1,
    method_2 = method_2,
    MAE = mean(abs(x - y)),
    RMSE = sqrt(mean((x - y)^2)),
    correlation = cor(x, y)
  )
}

method_pairs <- combn(methods, 2, simplify = FALSE)

initial_metrics <- do.call(rbind, lapply(method_pairs, function(pair) {
  metric_row(
    initial_mu[[pair[1]]],
    initial_mu[[pair[2]]],
    pair[1],
    pair[2]
  )
}))

final_metrics <- do.call(rbind, lapply(method_pairs, function(pair) {
  metric_row(
    final_mu[[pair[1]]],
    final_mu[[pair[2]]],
    pair[1],
    pair[2]
  )
}))

propagation <- merge(
  initial_metrics,
  final_metrics,
  by = c("method_1", "method_2"),
  suffixes = c("_initial", "_final"),
  sort = FALSE
)
propagation$MAE_ratio_final_to_initial <-
  propagation$MAE_final / propagation$MAE_initial
propagation$RMSE_ratio_final_to_initial <-
  propagation$RMSE_final / propagation$RMSE_initial
propagation$MAE_attenuation_percentage <-
  100 * (1 - propagation$MAE_ratio_final_to_initial)
propagation$RMSE_attenuation_percentage <-
  100 * (1 - propagation$RMSE_ratio_final_to_initial)

write.csv(
  initial_metrics,
  "thesis_revision/verified_daily_grid_initial_background_metrics_Mc_2_0.csv",
  row.names = FALSE
)
write.csv(
  propagation,
  "thesis_revision/verified_daily_grid_propagation_Mc_2_0.csv",
  row.names = FALSE
)

# Refit the saved nearest-neighbour proximity values to persist the exact
# mixture-model threshold that the main script prints but does not save.
suppressPackageStartupMessages(library(mclust))
nn <- read.csv(
  "declustering_nearest_neighbor_Mc_2_0.csv",
  check.names = FALSE
)
log_eta <- nn$nn_log_eta[is.finite(nn$nn_log_eta)]
set.seed(123)
mixture_fit <- Mclust(
  log_eta,
  G = 2,
  modelNames = c("E", "V"),
  verbose = FALSE
)
threshold_grid <- seq(
  quantile(log_eta, 0.001),
  quantile(log_eta, 0.999),
  length.out = 4000
)
posterior <- predict(mixture_fit, newdata = threshold_grid)$z
threshold <- threshold_grid[which.min(abs(posterior[, 1] - posterior[, 2]))]
ordered_components <- order(as.numeric(mixture_fit$parameters$mean))
mixture_diagnostic <- data.frame(
  model = mixture_fit$modelName,
  threshold_log10_eta = threshold,
  lower_component_mean = as.numeric(mixture_fit$parameters$mean)[ordered_components[1]],
  upper_component_mean = as.numeric(mixture_fit$parameters$mean)[ordered_components[2]],
  lower_component_proportion = mixture_fit$parameters$pro[ordered_components[1]],
  upper_component_proportion = mixture_fit$parameters$pro[ordered_components[2]]
)
write.csv(
  mixture_diagnostic,
  "thesis_revision/verified_nn_mixture_diagnostic_Mc_2_0.csv",
  row.names = FALSE
)

message("Verified same-grid initialisation and propagation metrics written.")
