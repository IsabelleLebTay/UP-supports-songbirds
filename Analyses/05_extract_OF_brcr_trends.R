# Extract old-forest baseline occupancy and BRCR background trend summaries.
#
# Inputs:
#   Results/Model fits/Model_fit.rds
#   Data/full_data_for_stan.rds
#
# Outputs:
#   Results/Projections/OF_*.csv
#   Results/Projections/BRCR_change_scenarios.csv

library(tidyverse)
library(cmdstanr)
library(posterior)

projection_dir <- file.path("Results", "Projections")
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path("Results", "Model fits", "Model_fit.rds"))
full_data <- readRDS(file.path("Data", "full_data_for_stan.rds"))

species_labels <- full_data$species_map$species
n_species <- length(species_labels)
n_years <- full_data$stan_data$n_years
years <- full_data$year_base + seq_len(n_years) - 1L
OF <- 2L

clamp <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

summarise_vector <- function(x) {
  tibble(
    mean = mean(x),
    median = median(x),
    q2.5 = quantile(x, 0.025),
    q10 = quantile(x, 0.10),
    q90 = quantile(x, 0.90),
    q97.5 = quantile(x, 0.975)
  )
}

# Baseline OF occupancy is inv_logit(alpha_psi), because OF is the reference
# treatment in the initial occupancy parameterization.
alpha_psi_draws <- fit$draws("alpha_psi", format = "matrix")
colnames(alpha_psi_draws) <- species_labels
of_psi_draws <- plogis(alpha_psi_draws)

of_psi_summary <- map_dfr(species_labels, function(sp) {
  summarise_vector(of_psi_draws[, sp]) %>%
    mutate(species = sp, .before = 1)
})

write_csv(as_tibble(of_psi_draws), file.path(projection_dir, "OF_baseline_psi_draws.csv"))

# Full OF occupancy trajectory by species. TSH grid choice is arbitrary for OF
# because TSH effects are zero for the OF treatment.
of_traj_draws <- map(seq_len(n_species), function(sp) {
  vars <- paste0("occ_trajectory[", sp, ",", OF, ",1,", seq_len(n_years), "]")
  x <- fit$draws(vars, format = "matrix")
  colnames(x) <- paste0("year_", years)
  x
})
names(of_traj_draws) <- species_labels

n_draws <- nrow(of_traj_draws[[1]])
year_index <- seq_len(n_years)

of_logit_slope_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_species)
of_logit_intercept_draws <- matrix(NA_real_, nrow = n_draws, ncol = n_species)
colnames(of_logit_slope_draws) <- species_labels
colnames(of_logit_intercept_draws) <- species_labels

for (sp in seq_len(n_species)) {
  traj <- of_traj_draws[[sp]]
  for (draw in seq_len(n_draws)) {
    logit_occ <- qlogis(clamp(traj[draw, ]))
    coefs <- .lm.fit(cbind(1, year_index), logit_occ)$coefficients
    of_logit_intercept_draws[draw, sp] <- coefs[1]
    of_logit_slope_draws[draw, sp] <- coefs[2]
  }
}

of_trend_summary <- map_dfr(species_labels, function(sp) {
  summarise_vector(of_logit_slope_draws[, sp]) %>%
    mutate(
      species = sp,
      prob_decline = mean(of_logit_slope_draws[, sp] < 0),
      .before = 1
    )
})

write_csv(as_tibble(of_logit_slope_draws), file.path(projection_dir, "OF_logit_slope_draws.csv"))

of_annual_diffs <- map_dfr(seq_along(species_labels), function(sp) {
  diffs <- of_traj_draws[[sp]][, 2:n_years] - of_traj_draws[[sp]][, 1:(n_years - 1)]
  colnames(diffs) <- paste0(years[1:(n_years - 1)], "_to_", years[2:n_years])
  as_tibble(diffs) %>%
    mutate(species = species_labels[sp], draw = row_number(), .before = 1)
})

write_csv(of_annual_diffs, file.path(projection_dir, "OF_annual_first_differences.csv"))

of_pct_change <- matrix(NA_real_, nrow = n_draws, ncol = n_species)
colnames(of_pct_change) <- species_labels

for (sp in seq_len(n_species)) {
  occ_first <- of_traj_draws[[sp]][, 1]
  occ_last <- of_traj_draws[[sp]][, n_years]
  of_pct_change[, sp] <- ((occ_last - occ_first) / occ_first) * 100
}

of_change_summary <- map_dfr(seq_along(species_labels), function(sp) {
  tibble(
    species = species_labels[sp],
    occ_first_median = median(of_traj_draws[[sp]][, 1]),
    occ_first_q10 = quantile(of_traj_draws[[sp]][, 1], 0.10),
    occ_first_q90 = quantile(of_traj_draws[[sp]][, 1], 0.90),
    occ_last_median = median(of_traj_draws[[sp]][, n_years]),
    occ_last_q10 = quantile(of_traj_draws[[sp]][, n_years], 0.10),
    occ_last_q90 = quantile(of_traj_draws[[sp]][, n_years], 0.90),
    pct_change_median = median(of_pct_change[, sp]),
    pct_change_q10 = quantile(of_pct_change[, sp], 0.10),
    pct_change_q90 = quantile(of_pct_change[, sp], 0.90),
    prob_decline = mean(of_pct_change[, sp] < 0)
  )
})

write_csv(of_change_summary, file.path(projection_dir, "OF_endpoint_change_summary.csv"))

of_combined <- of_psi_summary %>%
  rename(
    psi_mean = mean,
    psi_median = median,
    psi_q2.5 = q2.5,
    psi_q10 = q10,
    psi_q90 = q90,
    psi_q97.5 = q97.5
  ) %>%
  left_join(
    of_trend_summary %>%
      rename(
        logit_slope_mean = mean,
        logit_slope_median = median,
        logit_slope_q2.5 = q2.5,
        logit_slope_q10 = q10,
        logit_slope_q90 = q90,
        logit_slope_q97.5 = q97.5,
        logit_slope_prob_decline = prob_decline
      ),
    by = "species"
  )

write_csv(of_combined, file.path(projection_dir, "OF_baseline_summary.csv"))

brcr_idx <- match("BRCR", species_labels)

brcr_pct_change <- of_pct_change[, brcr_idx]

brcr_decline_scenarios <- tibble(
  scenario = c("optimistic", "median", "pessimistic"),
  pct_change = c(
    quantile(brcr_pct_change, 0.90),
    median(brcr_pct_change),
    quantile(brcr_pct_change, 0.10)
  ),
  probability_scale_multiplier = 1 + pct_change / 100
)

write_csv(brcr_decline_scenarios, file.path(projection_dir, "BRCR_change_scenarios.csv"))