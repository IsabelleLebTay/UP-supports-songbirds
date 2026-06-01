# Extract posterior summaries from the final dynamic occupancy model.
#
# Inputs:
#   Results/Model fits/Model_fit.rds
#   Data/full_data_for_stan.rds
#
# Outputs:
#   Results/Summaries/*.csv

library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

summary_dir <- file.path("Results", "Summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(file.path("Results", "Model fits", "Model_fit.rds"))
full_data <- readRDS(file.path("Data", "full_data_for_stan.rds"))

species_labels <- full_data$species_map$species
treatment_labels <- c("CLH", "OF", "UP")
n_species <- length(species_labels)
n_treatments <- length(treatment_labels)
n_years <- full_data$stan_data$n_years
n_sites <- full_data$stan_data$n_sites
n_tsh_grid <- full_data$stan_data$n_tsh_grid
years <- full_data$year_base + seq_len(n_years) - 1L
tsh_values <- seq_len(n_tsh_grid) - 1L

ci_probs <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)

draw_var <- function(variable) {
  as.numeric(fit$draws(variable, format = "matrix")[, 1])
}

summarise_vector <- function(x) {
  qs <- quantile(x, probs = ci_probs, names = FALSE)
  tibble(
    mean = mean(x),
    sd = sd(x),
    q2.5 = qs[1],
    q5 = qs[2],
    q10 = qs[3],
    q25 = qs[4],
    median = qs[5],
    q75 = qs[6],
    q90 = qs[7],
    q95 = qs[8],
    q97.5 = qs[9],
    prob_positive = mean(x > 0),
    prob_negative = mean(x < 0)
  )
}

write_summary <- function(x, filename) {
  write_csv(x, file.path(summary_dir, filename))
}

# Convergence diagnostics -----------------------------------------------------
draw_summary <- posterior::summarise_draws(fit$draws())
sampler_diag <- fit$diagnostic_summary()

diagnostics <- tibble(
  metric = c(
    "max_rhat",
    "min_ess_bulk",
    "min_ess_tail",
    "n_rhat_above_1.01",
    "n_ess_bulk_below_400",
    "n_ess_tail_below_400",
    "n_divergent",
    "n_max_treedepth"
  ),
  value = c(
    max(draw_summary$rhat, na.rm = TRUE),
    min(draw_summary$ess_bulk, na.rm = TRUE),
    min(draw_summary$ess_tail, na.rm = TRUE),
    sum(draw_summary$rhat > 1.01, na.rm = TRUE),
    sum(draw_summary$ess_bulk < 400, na.rm = TRUE),
    sum(draw_summary$ess_tail < 400, na.rm = TRUE),
    sum(sampler_diag$num_divergent),
    sum(sampler_diag$num_max_treedepth)
  )
)

ebfmi <- tibble(
  chain = seq_along(sampler_diag$ebfmi),
  ebfmi = sampler_diag$ebfmi
)

write_summary(diagnostics, "convergence_diagnostics.csv")
write_summary(ebfmi, "ebfmi_per_chain.csv")

# LOO by species --------------------------------------------------------------
loo_summary <- map_dfr(seq_len(n_species), function(sp) {
  ll_vars <- paste0("log_lik[", sp, ",", seq_len(n_sites), "]")
  ll_array <- fit$draws(ll_vars, format = "draws_array")
  loo_sp <- loo(ll_array, r_eff = relative_eff(exp(ll_array)))
  pareto_k <- loo_sp$diagnostics$pareto_k

  tibble(
    species = species_labels[sp],
    elpd_loo = loo_sp$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = loo_sp$estimates["elpd_loo", "SE"],
    p_loo = loo_sp$estimates["p_loo", "Estimate"],
    n_pareto_k_high = sum(pareto_k > 0.7),
    n_pareto_k_very_high = sum(pareto_k > 1.0),
    max_pareto_k = max(pareto_k)
  )
})
write_summary(loo_summary, "loo_by_species.csv")

# Occupancy, treatment, and covariate effects --------------------------------
psi_summary <- map_dfr(seq_len(n_species), function(sp) {
  map_dfr(seq_len(n_treatments), function(trt) {
    summarise_vector(draw_var(paste0("psi_baseline[", sp, ",", trt, "]"))) %>%
      mutate(
        species = species_labels[sp],
        treatment = treatment_labels[trt],
        parameter = "psi_baseline",
        .before = 1
      )
  })
})

psi_contrasts <- crossing(
  species_id = seq_len(n_species),
  contrast = c("psi_UP_vs_CLH", "psi_UP_vs_OF", "psi_CLH_vs_OF")
) %>%
  mutate(summary = map2(species_id, contrast, \(sp, contrast) {
    summarise_vector(draw_var(paste0(contrast, "[", sp, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(species = species_labels[species_id], contrast, across(mean:prob_negative))

write_summary(psi_summary, "psi_baseline_by_species_treatment.csv")
write_summary(psi_contrasts, "psi_contrasts_by_species.csv")

tsh_slopes <- crossing(
  species_id = seq_len(n_species),
  treatment_id = seq_len(n_treatments),
  process = c("psi", "phi", "gam")
) %>%
  mutate(summary = pmap(list(species_id, treatment_id, process), \(sp, trt, proc) {
    summarise_vector(draw_var(paste0("beta_tsh_", proc, "_out[", sp, ",", trt, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(
    species = species_labels[species_id],
    treatment = treatment_labels[treatment_id],
    process,
    across(mean:prob_negative)
  )

tsh_contrasts <- crossing(
  species_id = seq_len(n_species),
  process = c("psi", "phi", "gam")
) %>%
  mutate(summary = map2(species_id, process, \(sp, proc) {
    summarise_vector(draw_var(paste0("beta_tsh_", proc, "_UP_vs_CLH[", sp, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(
    species = species_labels[species_id],
    process,
    contrast = "UP_vs_CLH",
    across(mean:prob_negative)
  )

tsh_community <- crossing(
  process = c("psi", "phi", "gam"),
  treatment = c("CLH", "UP"),
  statistic = c("mean", "sd")
) %>%
  mutate(
    variable = if_else(
      statistic == "mean",
      paste0("mu_beta_tsh_", process, "_", treatment),
      paste0("sigma_beta_tsh_", process, "_", treatment)
    ),
    summary = map(variable, \(x) summarise_vector(draw_var(x)))
  ) %>%
  unnest(summary)

write_summary(tsh_slopes, "tsh_slopes_by_species_treatment_process.csv")
write_summary(tsh_contrasts, "tsh_contrasts_UP_vs_CLH.csv")
write_summary(tsh_community, "tsh_community_hyperparams.csv")

dte_slopes <- crossing(
  species_id = seq_len(n_species),
  treatment_id = seq_len(n_treatments)
) %>%
  mutate(summary = map2(species_id, treatment_id, \(sp, trt) {
    summarise_vector(draw_var(paste0("beta_dte_out[", sp, ",", trt, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(
    species = species_labels[species_id],
    treatment = treatment_labels[treatment_id],
    across(mean:prob_negative)
  )

dte_contrasts <- crossing(
  species_id = seq_len(n_species),
  contrast = c("beta_dte_UP_vs_CLH", "beta_dte_UP_vs_OF", "beta_dte_CLH_vs_OF")
) %>%
  mutate(summary = map2(species_id, contrast, \(sp, contrast) {
    summarise_vector(draw_var(paste0(contrast, "[", sp, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(species = species_labels[species_id], contrast, across(mean:prob_negative))

uheight_slopes <- map_dfr(seq_len(n_species), function(sp) {
  summarise_vector(draw_var(paste0("beta_uheight_psi_out[", sp, "]"))) %>%
    mutate(species = species_labels[sp], parameter = "beta_uheight_psi", .before = 1)
})

uheight_community <- tibble(
  parameter = c("mu_beta_uheight_psi", "sigma_beta_uheight_psi")
) %>%
  mutate(summary = map(parameter, \(x) summarise_vector(draw_var(x)))) %>%
  unnest(summary)

write_summary(dte_slopes, "dte_slopes_by_species_treatment.csv")
write_summary(dte_contrasts, "dte_contrasts_by_species.csv")
write_summary(uheight_slopes, "uheight_slopes_by_species.csv")
write_summary(uheight_community, "uheight_community_hyperparams.csv")

# Detection summaries ---------------------------------------------------------
detection_summary <- map_dfr(seq_len(n_species), function(sp) {
  tibble(
    species = species_labels[sp],
    parameter = c("p", "beta_doy_p", "beta_tod_p"),
    summary = list(
      summarise_vector(draw_var(paste0("p[", sp, "]"))),
      summarise_vector(draw_var(paste0("beta_doy_p[", sp, "]"))),
      summarise_vector(draw_var(paste0("beta_tod_p[", sp, "]")))
    )
  ) %>%
    unnest(summary)
})

detection_community <- tibble(
  parameter = c(
    "mu_alpha_p", "sigma_alpha_p",
    "mu_beta_doy_p", "sigma_beta_doy_p",
    "mu_beta_tod_p", "sigma_beta_tod_p"
  )
) %>%
  mutate(summary = map(parameter, \(x) summarise_vector(draw_var(x)))) %>%
  unnest(summary)

write_summary(detection_summary, "detection_by_species.csv")
write_summary(detection_community, "detection_community_hyperparams.csv")

# Temporal and random-effect summaries ---------------------------------------
year_effects <- crossing(
  treatment_id = seq_len(n_treatments),
  year_id = seq_len(n_years),
  process = c("phi", "gam")
) %>%
  mutate(summary = pmap(list(treatment_id, year_id, process), \(trt, yr, proc) {
    summarise_vector(draw_var(paste0("year_", proc, "_shared[", trt, ",", yr, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(
    treatment = treatment_labels[treatment_id],
    year = years[year_id],
    process,
    across(mean:prob_negative)
  )

rw_variances <- crossing(
  treatment_id = seq_len(n_treatments),
  process = c("phi", "gam")
) %>%
  mutate(summary = map2(treatment_id, process, \(trt, proc) {
    summarise_vector(draw_var(paste0("sigma_year_", proc, "[", trt, "]")))
  })) %>%
  unnest(summary) %>%
  transmute(
    treatment = treatment_labels[treatment_id],
    process,
    across(mean:prob_negative)
  )

species_rw_variances <- tibble(
  parameter = c("sigma_species_year_phi", "sigma_species_year_gam")
) %>%
  mutate(summary = map(parameter, \(x) summarise_vector(draw_var(x)))) %>%
  unnest(summary)

cluster_summary <- map_dfr(seq_len(n_species), function(sp) {
  summarise_vector(draw_var(paste0("sigma_cluster_out[", sp, "]"))) %>%
    mutate(species = species_labels[sp], parameter = "sigma_cluster", .before = 1)
})

intercept_community <- tibble(
  parameter = c(
    "mu_alpha_psi", "sigma_alpha_psi",
    "mu_alpha_phi", "sigma_alpha_phi",
    "mu_alpha_gam", "sigma_alpha_gam"
  )
) %>%
  mutate(summary = map(parameter, \(x) summarise_vector(draw_var(x)))) %>%
  unnest(summary)

write_summary(year_effects, "random_walk_year_effects.csv")
write_summary(rw_variances, "random_walk_variances.csv")
write_summary(species_rw_variances, "species_random_walk_variances.csv")
write_summary(cluster_summary, "cluster_sd_by_species.csv")
write_summary(intercept_community, "intercept_community_hyperparams.csv")

# Observed detection summaries ------------------------------------------------
obs_detections <- full_data$recordings %>%
  group_by(species) %>%
  summarise(
    n_recordings = n(),
    n_detections = sum(detected),
    detection_rate = n_detections / n_recordings,
    .groups = "drop"
  )
write_summary(obs_detections, "observed_detection_summary.csv")

# Year-averaged occupancy trajectories ---------------------------------------
occ_trajectory_summary <- crossing(
  species_id = seq_len(n_species),
  treatment_id = seq_len(n_treatments),
  tsh_id = seq_len(n_tsh_grid)
) %>%
  mutate(summary = pmap(list(species_id, treatment_id, tsh_id), \(sp, trt, tsh_id) {
    variables <- paste0("occ_trajectory[", sp, ",", trt, ",", tsh_id, ",", seq_len(n_years), "]")
    draws <- fit$draws(variables, format = "matrix")
    summarise_vector(rowMeans(draws))
  })) %>%
  unnest(summary) %>%
  transmute(
    species = species_labels[species_id],
    treatment = treatment_labels[treatment_id],
    tsh = tsh_values[tsh_id],
    across(mean:prob_negative)
  )

write_summary(occ_trajectory_summary, "occ_trajectory_by_species_treatment_tsh.csv")