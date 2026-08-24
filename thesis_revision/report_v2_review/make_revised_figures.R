#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

output_dir <- "thesis_revision/report_v2_review/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_thesis <- function() {
  theme_minimal(base_size = 11.5) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(colour = "#222222"),
      axis.text = element_text(colour = "#333333"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F1F3F5", colour = NA),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(8, 14, 8, 8)
    )
}

save_plot_pair <- function(plot, stem, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

# Figure 1: annual counts with the high-activity year made explicit.
annual <- read.csv("annual_event_counts.csv", check.names = FALSE)
annual$period <- ifelse(annual$year == 2019, "2019", "Other years")

p_annual <- ggplot(annual, aes(x = year, y = n_events, fill = period)) +
  geom_col(width = 0.82) +
  geom_text(
    data = annual[annual$year == 2019, ],
    aes(label = format(n_events, big.mark = ",", scientific = FALSE)),
    vjust = -0.45,
    size = 3.6,
    fontface = "bold",
    colour = "#9C3B00"
  ) +
  scale_fill_manual(values = c("2019" = "#D55E00", "Other years" = "#4477AA")) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 5), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = "Year",
    y = "Recorded earthquakes",
    caption = "The years 2000 and 2026 are partial observation years."
  ) +
  theme_thesis() +
  theme(legend.position = "none")

save_plot_pair(p_annual, "annual_event_counts_revised", 7.2, 4.1)

# Figure 2: rolling completeness estimates with directly labelled references.
rolling <- read.csv("Temporal_Mc_Rolling_Window.csv", check.names = FALSE)
rolling$center_time <- as.POSIXct(rolling$center_time, tz = "UTC")

p_mc <- ggplot(rolling, aes(x = center_time, y = mc_mbs)) +
  annotate(
    "rect",
    xmin = as.POSIXct("2019-01-01", tz = "UTC"),
    xmax = as.POSIXct("2019-12-31 23:59:59", tz = "UTC"),
    ymin = -Inf,
    ymax = Inf,
    fill = "#E69F00",
    alpha = 0.13
  ) +
  geom_line(linewidth = 0.42, colour = "#3B4C5A", alpha = 0.72) +
  geom_point(size = 0.72, colour = "#176B87", alpha = 0.72) +
  geom_hline(yintercept = 2.0, linetype = "dashed", linewidth = 0.55, colour = "#9C3B00") +
  geom_hline(yintercept = 1.7, linetype = "dotted", linewidth = 0.55, colour = "#555555") +
  annotate(
    "text",
    x = as.POSIXct("2025-09-01", tz = "UTC"),
    y = 2.06,
    label = "Operational threshold: 2.0",
    hjust = 1,
    vjust = 0,
    size = 3.15,
    colour = "#9C3B00"
  ) +
  annotate(
    "text",
    x = as.POSIXct("2025-09-01", tz = "UTC"),
    y = 1.64,
    label = "Catalogue-wide MBS-WW: 1.7",
    hjust = 1,
    vjust = 1,
    size = 3.15,
    colour = "#444444"
  ) +
  annotate(
    "text",
    x = as.POSIXct("2019-07-01", tz = "UTC"),
    y = 3.48,
    label = "2019",
    size = 3.2,
    fontface = "bold",
    colour = "#9C6500"
  ) +
  scale_x_datetime(date_breaks = "5 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(breaks = seq(0.5, 3.5, by = 0.5), limits = c(0.45, 3.58)) +
  labs(
    x = "Window midpoint",
    y = expression("MBS-WW estimate of " * M[c]),
    caption = "Each point uses 1,000 successive events; windows advance by 250 events."
  ) +
  theme_thesis()

save_plot_pair(p_mc, "rolling_mc_revised", 7.2, 4.4)

# Figure 3: reconstruct the three smooth starts on the final daily grid, then
# contrast them with the saved final background functions. This is a visual
# representation of the already verified same-grid comparison, not a new fit.
monthly <- read.csv(
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

make_initial_curve <- function(method_name) {
  z <- monthly[monthly$method == method_name, ]
  z <- z[order(z$month), ]
  month_start <- as.POSIXct(z$month, tz = "UTC")
  month_sequence <- seq(
    as.Date(format(min(month_start), "%Y-%m-01")),
    by = "month",
    length.out = nrow(z) + 1
  )
  month_end <- as.POSIXct(month_sequence[-1], tz = "UTC")
  observed_start <- pmax(month_start, t_origin)
  observed_end <- pmin(month_end, t_end_datetime)
  exposure_days <- as.numeric(difftime(observed_end, observed_start, units = "days"))
  midpoint <- observed_start + (observed_end - observed_start) / 2
  midpoint_days <- as.numeric(difftime(midpoint, t_origin, units = "days"))
  y <- log((z$background_events + 0.1) / pmax(exposure_days, 1))
  B_midpoint <- predict(B_event, newx = midpoint_days)
  beta <- solve(
    crossprod(B_midpoint) + 1e-4 * diag(ncol(B_midpoint)),
    crossprod(B_midpoint, y)
  )
  data.frame(
    datetime = t_origin + t_grid * 86400,
    rate = as.numeric(exp(B_grid %*% beta)),
    method = method_name,
    stage = "Smooth starting functions"
  )
}

methods <- c("Gardner-Knopoff", "Nearest-neighbour", "Reasenberg-style")
initial_curves <- do.call(rbind, lapply(methods, make_initial_curve))
final_curves <- data.frame(
  datetime = as.POSIXct(final$datetime, tz = "UTC"),
  rate = final$mu,
  method = final$method,
  stage = "Final fitted functions"
)
curves <- rbind(initial_curves, final_curves)
curves$stage <- factor(
  curves$stage,
  levels = c("Smooth starting functions", "Final fitted functions")
)

write.csv(
  curves,
  file.path(output_dir, "initial_and_final_background_curves_daily.csv"),
  row.names = FALSE
)

method_colours <- c(
  "Gardner-Knopoff" = "#0072B2",
  "Nearest-neighbour" = "#D55E00",
  "Reasenberg-style" = "#009E73"
)

p_curves <- ggplot(
  curves,
  aes(x = datetime, y = rate, colour = method, linetype = method)
) +
  annotate(
    "rect",
    xmin = as.POSIXct("2019-01-01", tz = "UTC"),
    xmax = as.POSIXct("2019-12-31 23:59:59", tz = "UTC"),
    ymin = -Inf,
    ymax = Inf,
    fill = "#E69F00",
    alpha = 0.08
  ) +
  geom_line(linewidth = 0.78, alpha = 0.88) +
  geom_label(
    data = data.frame(
      datetime = as.POSIXct("2022-06-01", tz = "UTC"),
      rate = 1.18,
      stage = factor(
        "Final fitted functions",
        levels = c("Smooth starting functions", "Final fitted functions")
      ),
      label = "Three fitted curves\nare almost coincident"
    ),
    aes(x = datetime, y = rate, label = label),
    inherit.aes = FALSE,
    size = 3.0,
    linewidth = 0.18,
    label.padding = grid::unit(0.15, "lines"),
    colour = "#333333",
    fill = "white"
  ) +
  facet_wrap(~stage, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = method_colours) +
  scale_linetype_manual(
    values = c(
      "Gardner-Knopoff" = "solid",
      "Nearest-neighbour" = "longdash",
      "Reasenberg-style" = "dotdash"
    )
  ) +
  scale_x_datetime(date_breaks = "5 years", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  labs(x = "Date", y = "Background rate (events/day)") +
  theme_thesis()

save_plot_pair(p_curves, "initial_vs_final_background_revised", 7.2, 6.4)

message("Revised thesis figures written to ", output_dir)
