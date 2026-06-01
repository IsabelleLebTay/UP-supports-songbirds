# Project BRCR habitat supply under alternative harvest scenarios.
#
# Inputs:
#   Results/Model fits/Model_fit.rds
#   Data/full_data_for_stan.rds
#
# Outputs:
#   Results/Projections/habitat_supply_projections.csv
#   Results/Projections/occ_avg_by_treatment_tsh.csv

library(tidyverse)
library(cmdstanr)
library(posterior)

options(width = as.numeric(Sys.getenv("COLUMNS", unset = "100")))

projection_dir <- file.path("Results", "Projections")
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path("Results", "Model fits", "Model_fit.rds"))
full_data <- readRDS(file.path("Data", "full_data_for_stan.rds"))

species_labels <- full_data$species_map$species
focal_species <- "BRCR"
focal_idx <- match(focal_species, species_labels)

n_years <- full_data$stan_data$n_years
tsh_mean <- full_data$scaling$tsh_mean
tsh_sd <- full_data$scaling$tsh_sd
dte_ref <- full_data$stan_data$dte_ref
uheight_ref <- full_data$stan_data$uheight_ref

CLH <- 1L
OF <- 2L
UP <- 3L
treatment_labels <- c("CLH", "OF", "UP")

# Landbase constants are derived from proprietary Al-Pac landbase summaries.
# Public data needed to recompute the model is in Data/; these constants define
# the management scenarios used for the habitat-supply projection.
landbase <- list(
  of_ha_initial = 535912,
  clh_ha_initial = 34211,
  up_ha_initial = 5252,
  of_annual_loss = 786,
  max_tsh = 22L
)

total_harvest_ha <- landbase$clh_ha_initial + landbase$up_ha_initial
clh_pct_current <- landbase$clh_ha_initial / total_harvest_ha
up_pct_current <- landbase$up_ha_initial / total_harvest_ha
harvest_annual <- total_harvest_ha / landbase$max_tsh

scenarios <- tribble(
  ~scenario, ~clh_pct_new, ~up_pct_new,
  "S1: Current conditions", clh_pct_current, up_pct_current,
  "S2: No UP", 1.00, 0.00,
  "S3: Maximum feasible UP", 0.80, 0.20
)

draw_var <- function(variable) {
  as.numeric(fit$draws(variable, format = "matrix")[, 1])
}

extract_species_treatment_array <- function(prefix, focal_idx) {
  n_draws <- nrow(fit$draws("alpha_psi", format = "matrix"))
  out <- array(NA_real_, dim = c(n_draws, 3))
  for (trt in 1:3) {
    out[, trt] <- draw_var(paste0(prefix, "[", focal_idx, ",", trt, "]"))
  }
  out
}

alpha_psi <- draw_var(paste0("alpha_psi[", focal_idx, "]"))
alpha_phi <- draw_var(paste0("alpha_phi[", focal_idx, "]"))
alpha_gam <- draw_var(paste0("alpha_gam[", focal_idx, "]"))
n_draws <- length(alpha_psi)

beta_trt_psi <- extract_species_treatment_array("beta_trt_psi", focal_idx)
beta_tsh_psi <- extract_species_treatment_array("beta_tsh_psi_out", focal_idx)
beta_tsh_phi <- extract_species_treatment_array("beta_tsh_phi_out", focal_idx)
beta_tsh_gam <- extract_species_treatment_array("beta_tsh_gam_out", focal_idx)
beta_dte <- extract_species_treatment_array("beta_dte_out", focal_idx)
beta_uheight_psi <- draw_var(paste0("beta_uheight_psi_out[", focal_idx, "]"))

tsh_ages <- 0:(landbase$max_tsh - 1L)
tsh_z <- outer(tsh_ages, seq_len(n_years), function(age, year_id) {
  (age + year_id - 1L - tsh_mean) / tsh_sd
})

occ_avg <- array(NA_real_, dim = c(n_draws, 3, length(tsh_ages)))

for (trt in 1:3) {
  for (age_id in seq_along(tsh_ages)) {
    logit_psi <- alpha_psi +
      beta_trt_psi[, trt] +
      beta_tsh_psi[, trt] * tsh_z[age_id, 1] +
      beta_dte[, trt] * dte_ref[trt]

    if (trt == UP) {
      logit_psi <- logit_psi + beta_uheight_psi * uheight_ref
    }

    occ_t <- plogis(logit_psi)
    occ_sum <- occ_t

    for (year_id in 2:n_years) {
      phi_t <- plogis(
        alpha_phi +
          beta_tsh_phi[, trt] * tsh_z[age_id, year_id] +
          beta_dte[, trt] * dte_ref[trt]
      )
      gam_t <- plogis(
        alpha_gam +
          beta_tsh_gam[, trt] * tsh_z[age_id, year_id] +
          beta_dte[, trt] * dte_ref[trt]
      )
      occ_t <- occ_t * phi_t + (1 - occ_t) * gam_t
      occ_sum <- occ_sum + occ_t
    }

    occ_avg[, trt, age_id] <- occ_sum / n_years
  }
}

compute_supply <- function(clh_cohorts, up_cohorts, of_ha, occ_avg) {
  harvest_supply <- rep(0, n_draws)

  for (age_id in seq_along(tsh_ages)) {
    harvest_supply <- harvest_supply +
      clh_cohorts[age_id] * occ_avg[, CLH, age_id] +
      up_cohorts[age_id] * occ_avg[, UP, age_id]
  }

  of_ha * occ_avg[, OF, 1] + harvest_supply
}

run_scenario <- function(clh_pct_new, up_pct_new) {
  clh_cohorts <- rep(landbase$clh_ha_initial / landbase$max_tsh, landbase$max_tsh)
  up_cohorts <- rep(landbase$up_ha_initial / landbase$max_tsh, landbase$max_tsh)
  of_ha <- landbase$of_ha_initial
  proj_years <- 0:landbase$max_tsh

  habitat_supply <- matrix(NA_real_, nrow = n_draws, ncol = length(proj_years))
  habitat_supply[, 1] <- compute_supply(clh_cohorts, up_cohorts, of_ha, occ_avg)

  for (year_id in 2:length(proj_years)) {
    of_ha <- of_ha - landbase$of_annual_loss
    clh_new <- harvest_annual * clh_pct_new
    up_new <- harvest_annual * up_pct_new

    clh_cohorts <- c(clh_new, clh_cohorts[1:(landbase$max_tsh - 1L)])
    up_cohorts <- c(up_new, up_cohorts[1:(landbase$max_tsh - 1L)])

    habitat_supply[, year_id] <- compute_supply(clh_cohorts, up_cohorts, of_ha, occ_avg)
  }

  habitat_supply
}

scenario_results <- pmap(
  list(scenarios$clh_pct_new, scenarios$up_pct_new),
  run_scenario
)
names(scenario_results) <- scenarios$scenario

summarise_scenario <- function(habitat_supply, scenario_name) {
  map_dfr(seq_len(ncol(habitat_supply)), function(year_id) {
    draws <- habitat_supply[, year_id]
    tibble(
      scenario = scenario_name,
      species = focal_species,
      proj_year = year_id - 1L,
      mean_ohe = mean(draws),
      median_ohe = median(draws),
      lower_ohe = quantile(draws, 0.025),
      upper_ohe = quantile(draws, 0.975),
      lower_80 = quantile(draws, 0.10),
      upper_80 = quantile(draws, 0.90)
    )
  })
}

all_summaries <- map2_dfr(scenario_results, names(scenario_results), summarise_scenario)
write_csv(all_summaries, file.path(projection_dir, "habitat_supply_projections.csv"))

occ_avg_summary <- crossing(
  species = focal_species,
  treatment = treatment_labels,
  tsh_age = tsh_ages
) %>%
  mutate(
    treatment_id = match(treatment, treatment_labels),
    age_id = match(tsh_age, tsh_ages),
    mean_occ = map2_dbl(treatment_id, age_id, \(trt, age) mean(occ_avg[, trt, age])),
    median_occ = map2_dbl(treatment_id, age_id, \(trt, age) median(occ_avg[, trt, age])),
    lower_occ = map2_dbl(treatment_id, age_id, \(trt, age) quantile(occ_avg[, trt, age], 0.025)),
    upper_occ = map2_dbl(treatment_id, age_id, \(trt, age) quantile(occ_avg[, trt, age], 0.975))
  ) %>%
  select(-treatment_id, -age_id)

write_csv(occ_avg_summary, file.path(projection_dir, "occ_avg_by_treatment_tsh.csv"))

summary_table <- all_summaries %>%
  filter(proj_year %in% c(0, landbase$max_tsh)) %>%
  select(scenario, species, proj_year, median_ohe) %>%
  pivot_wider(names_from = proj_year, values_from = median_ohe, names_prefix = "year_") %>%
  mutate(
    change_ohe = .data[[paste0("year_", landbase$max_tsh)]] - year_0,
    change_pct = 100 * change_ohe / year_0
  )

write_csv(summary_table, file.path(projection_dir, "habitat_supply_summary_table.csv"))