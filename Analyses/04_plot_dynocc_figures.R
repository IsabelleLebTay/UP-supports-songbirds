# Create manuscript and supplement figures from the final dynamic occupancy model.
#
# Inputs:
#   Results/Model fits/Model_fit.rds
#   Data/full_data_for_stan.rds
#
# Outputs:
#   Results/Figures/*.jpeg

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)

figure_dir <- file.path("Results", "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path("Results", "Model fits", "Model_fit.rds"))
full_data <- readRDS(file.path("Data", "full_data_for_stan.rds"))

species_labels <- full_data$species_map$species
species_order <- c("TEWA", "WIWR", "HETH", "BRCR", "BTNW", "BBWA")
species_order <- species_order[species_order %in% species_labels]
treatment_labels <- c("CLH", "OF", "UP")
treatment_order <- c("OF", "UP", "CLH")
years <- full_data$year_base + seq_len(full_data$stan_data$n_years) - 1L
tsh_values <- seq_len(full_data$stan_data$n_tsh_grid) - 1L

treatment_colours <- c(
  "CLH" = "#D55E00",
  "OF" = "#009E73",
  "UP" = "#0072B2"
)

species_colours <- c(
  "BBWA" = "#D55E00",
  "BRCR" = "#7F3C8D",
  "BTNW" = "#E69F00",
  "HETH" = "#009E73",
  "TEWA" = "#0072B2",
  "WIWR" = "#CC79A7"
)

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

draw_var <- function(variable) {
  as.numeric(fit$draws(variable, format = "matrix")[, 1])
}

summarise_draws <- function(x, lower_prob = 0.10, upper_prob = 0.90) {
  tibble(
    median = median(x),
    lower = quantile(x, lower_prob),
    upper = quantile(x, upper_prob),
    lower95 = quantile(x, 0.025),
    upper95 = quantile(x, 0.975)
  )
}

# Figure 2: occupancy versus time since harvest ------------------------------
occ_tsh_summary <- crossing(
  species_id = seq_along(species_labels),
  treatment_id = seq_along(treatment_labels),
  tsh_id = seq_along(tsh_values)
) %>%
  mutate(summary = pmap(list(species_id, treatment_id, tsh_id), \(sp, trt, tsh_id) {
    vars <- paste0(
      "occ_trajectory[", sp, ",", trt, ",", tsh_id, ",",
      seq_len(full_data$stan_data$n_years), "]"
    )
    summarise_draws(rowMeans(fit$draws(vars, format = "matrix")))
  })) %>%
  unnest(summary) %>%
  mutate(
    species = factor(species_labels[species_id], levels = species_order),
    treatment = factor(treatment_labels[treatment_id], levels = treatment_order),
    tsh = tsh_values[tsh_id]
  )

figure_2 <- ggplot(occ_tsh_summary, aes(x = tsh, y = median, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_x_continuous(breaks = seq(0, max(tsh_values), by = 5)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Time since harvest (years)",
    y = "Year-averaged occupancy probability"
  ) +
  theme_pub()

ggsave(
  file.path(figure_dir, "Figure_2_occ_vs_tsh.jpeg"),
  figure_2,
  width = 8,
  height = 5,
  dpi = 300
)

# Figure 3: occupancy through study years by treatment ------------------------
project_occupancy_by_treatment <- function(sp, trt) {
  psi <- draw_var(paste0("psi_baseline[", sp, ",", trt, "]"))
  n_draws <- length(psi)
  psi_proj <- matrix(NA_real_, nrow = n_draws, ncol = length(years))
  psi_proj[, 1] <- psi

  for (yr in 2:length(years)) {
    phi <- draw_var(paste0("phi_by_trt_year[", sp, ",", trt, ",", yr, "]"))
    gam <- draw_var(paste0("gam_by_trt_year[", sp, ",", trt, ",", yr, "]"))
    psi_proj[, yr] <- psi_proj[, yr - 1] * phi + (1 - psi_proj[, yr - 1]) * gam
  }

  map_dfr(seq_along(years), function(yr) {
    summarise_draws(psi_proj[, yr]) %>%
      mutate(year = years[yr])
  })
}

occ_year_summary <- crossing(
  species_id = seq_along(species_labels),
  treatment_id = seq_along(treatment_labels)
) %>%
  mutate(summary = map2(species_id, treatment_id, project_occupancy_by_treatment)) %>%
  unnest(summary) %>%
  mutate(
    species = factor(species_labels[species_id], levels = species_order),
    treatment = factor(treatment_labels[treatment_id], levels = treatment_order)
  )

figure_3 <- ggplot(occ_year_summary, aes(x = year, y = median, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_x_continuous(breaks = seq(min(years), max(years), by = 2)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Year",
    y = "Occupancy probability"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(figure_dir, "Figure_3_occ_by_treatment.jpeg"),
  figure_3,
  width = 8,
  height = 5,
  dpi = 300
)

# Supplement: baseline psi, persistence, and colonization ---------------------
baseline_summary <- crossing(
  species_id = seq_along(species_labels),
  treatment_id = seq_along(treatment_labels),
  process = c("psi", "phi", "gam")
) %>%
  mutate(summary = pmap(list(species_id, treatment_id, process), \(sp, trt, process) {
    if (process == "psi") {
      summarise_draws(draw_var(paste0("psi_baseline[", sp, ",", trt, "]")), 0.025, 0.975)
    } else {
      vars <- paste0(process, "_by_trt_year[", sp, ",", trt, ",", seq_along(years), "]")
      summarise_draws(rowMeans(fit$draws(vars, format = "matrix")), 0.025, 0.975)
    }
  })) %>%
  unnest(summary) %>%
  mutate(
    species = factor(species_labels[species_id], levels = rev(species_order)),
    treatment = factor(treatment_labels[treatment_id], levels = treatment_order),
    process = factor(
      process,
      levels = c("psi", "phi", "gam"),
      labels = c("Initial occupancy", "Persistence", "Colonization")
    )
  )

figure_s1 <- ggplot(baseline_summary, aes(x = median, y = species, colour = treatment)) +
  geom_linerange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.65), alpha = 0.65) +
  geom_linerange(aes(xmin = lower95, xmax = upper95), position = position_dodge(width = 0.65), alpha = 0.35) +
  geom_point(position = position_dodge(width = 0.65), size = 1.8) +
  facet_wrap(~process, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = treatment_colours) +
  labs(x = "Probability", y = NULL) +
  theme_pub()

ggsave(
  file.path(figure_dir, "Figure_S1_baseline_params.jpeg"),
  figure_s1,
  width = 8,
  height = 4,
  dpi = 300
)

# Supplement: detection probability by time of day and day of year ------------
summarise_detection_curve <- function(covariate = c("tod", "doy"), n_points = 120) {
  covariate <- match.arg(covariate)
  raw_values <- full_data$recordings[[covariate]]
  value_seq <- seq(min(raw_values, na.rm = TRUE), max(raw_values, na.rm = TRUE), length.out = n_points)
  scaled_seq <- (value_seq - full_data$scaling[[paste0(covariate, "_mean")]]) /
    full_data$scaling[[paste0(covariate, "_sd")]]

  map_dfr(seq_along(species_labels), function(sp) {
    alpha_p <- draw_var(paste0("alpha_p[", sp, "]"))
    beta <- draw_var(paste0("beta_", covariate, "_p[", sp, "]"))
    probs <- sapply(scaled_seq, function(x) plogis(alpha_p + beta * x))

    tibble(
      species = species_labels[sp],
      value = value_seq,
      median = apply(probs, 2, median),
      lower = apply(probs, 2, quantile, 0.10),
      upper = apply(probs, 2, quantile, 0.90),
      covariate = covariate
    )
  })
}

detection_curves <- bind_rows(
  summarise_detection_curve("tod"),
  summarise_detection_curve("doy")
) %>%
  mutate(
    species = factor(species, levels = species_order),
    covariate = recode(covariate, tod = "Time of day", doy = "Day of year")
  )

figure_s2 <- ggplot(detection_curves, aes(x = value, y = median, colour = species, fill = species)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.75) +
  facet_wrap(~covariate, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = species_colours) +
  scale_fill_manual(values = species_colours) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Detection probability") +
  theme_pub()

ggsave(
  file.path(figure_dir, "Figure_S2_detection_tod_doy.jpeg"),
  figure_s2,
  width = 9,
  height = 4,
  dpi = 300
)

# Supplement: distance-to-edge marginal effect --------------------------------
dte_grid <- seq(
  min(full_data$stan_data$dte, na.rm = TRUE),
  max(full_data$stan_data$dte, na.rm = TRUE),
  length.out = 80
)

dte_summary <- crossing(
  species_id = seq_along(species_labels),
  treatment_id = seq_along(treatment_labels),
  dte = dte_grid
) %>%
  mutate(summary = pmap(list(species_id, treatment_id, dte), \(sp, trt, dte) {
    psi <- plogis(
      draw_var(paste0("alpha_psi[", sp, "]")) +
        draw_var(paste0("beta_trt_psi[", sp, ",", trt, "]")) +
        draw_var(paste0("beta_dte_out[", sp, ",", trt, "]")) * dte
    )
    summarise_draws(psi)
  })) %>%
  unnest(summary) %>%
  mutate(
    species = factor(species_labels[species_id], levels = species_order),
    treatment = factor(treatment_labels[treatment_id], levels = treatment_order)
  )

figure_s3 <- ggplot(dte_summary, aes(x = dte, y = median, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Distance to edge (standardized)", y = "Initial occupancy probability") +
  theme_pub()

ggsave(
  file.path(figure_dir, "Figure_S3_distance_to_edge.jpeg"),
  figure_s3,
  width = 8,
  height = 5,
  dpi = 300
)

# Supplement: understory height marginal effect for UP sites ------------------
up_uheight_scaled <- full_data$stan_data$uheight[full_data$stan_data$treatment == 3]
uheight_grid_scaled <- seq(min(up_uheight_scaled, na.rm = TRUE), max(up_uheight_scaled, na.rm = TRUE), length.out = 80)
uheight_grid_raw <- uheight_grid_scaled * full_data$scaling$uheight_sd + full_data$scaling$uheight_mean

uheight_summary <- crossing(
  species_id = seq_along(species_labels),
  grid_id = seq_along(uheight_grid_scaled)
) %>%
  mutate(summary = map2(species_id, grid_id, \(sp, grid_id) {
    psi <- plogis(
      draw_var(paste0("alpha_psi[", sp, "]")) +
        draw_var(paste0("beta_trt_psi[", sp, ",3]")) +
        draw_var(paste0("beta_uheight_psi_out[", sp, "]")) * uheight_grid_scaled[grid_id]
    )
    summarise_draws(psi)
  })) %>%
  unnest(summary) %>%
  mutate(
    species = factor(species_labels[species_id], levels = species_order),
    uheight = uheight_grid_raw[grid_id]
  )

figure_s4 <- ggplot(uheight_summary, aes(x = uheight, y = median)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = treatment_colours[["UP"]], alpha = 0.18) +
  geom_line(colour = treatment_colours[["UP"]], linewidth = 0.8) +
  facet_wrap(~species, ncol = 3, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Understory height (m)", y = "Initial occupancy probability") +
  theme_pub()

ggsave(
  file.path(figure_dir, "Figure_S4_understory_height.jpeg"),
  figure_s4,
  width = 8,
  height = 5,
  dpi = 300
)

# Supplement: TSH slope diagnostics ------------------------------------------
tsh_slope_summary <- crossing(
  species_id = seq_along(species_labels),
  treatment_id = c(1L, 3L),
  process = c("psi", "phi", "gam")
) %>%
  mutate(summary = pmap(list(species_id, treatment_id, process), \(sp, trt, process) {
    summarise_draws(draw_var(paste0("beta_tsh_", process, "_out[", sp, ",", trt, "]")), 0.025, 0.975)
  })) %>%
  unnest(summary) %>%
  mutate(
    species = factor(species_labels[species_id], levels = rev(species_order)),
    treatment = factor(treatment_labels[treatment_id], levels = c("UP", "CLH")),
    process = factor(
      process,
      levels = c("psi", "phi", "gam"),
      labels = c("Initial occupancy", "Persistence", "Colonization")
    )
  )

figure_s5 <- ggplot(tsh_slope_summary, aes(x = median, y = species, colour = treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_linerange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.6), alpha = 0.7) +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  facet_wrap(~process, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = treatment_colours) +
  labs(x = "TSH effect on log-odds scale", y = NULL) +
  theme_pub()

ggsave(
  file.path(figure_dir, "Figure_S5_tsh_slopes.jpeg"),
  figure_s5,
  width = 8,
  height = 4,
  dpi = 300
)
