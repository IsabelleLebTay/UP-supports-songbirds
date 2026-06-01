# Fit the final multi-species dynamic occupancy model.
#
# Inputs:
#   Data/full_data_for_stan.rds
#   Code/Analyses/dynocc_multispecies.stan
#
# Outputs:
#   Results/Model fits/Multispp_height.rds
#   Results/Summaries/model_fit_diagnostics.csv

library(tidyverse)
library(cmdstanr)
library(posterior)

dotenv::load_dot_env(file = ".env")
cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN_PATH"))

fit_dir <- file.path("Results", "Model fits")
summary_dir <- file.path("Results", "Summaries")
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

full_data <- readRDS(file.path("Data", "full_data_for_stan.rds"))
stan_data <- full_data$stan_data

required_stan_fields <- c(
  "n_species", "n_sites", "n_visits", "n_obs", "n_treatments", "n_years",
  "n_clusters", "treatment", "cluster", "n_visits_site", "first_visit_idx",
  "tsh", "dte", "uheight", "visit_site", "visit_year", "n_recs_visit",
  "first_rec_idx", "y", "doy", "tod", "n_tsh_grid", "tsh_grid", "dte_ref",
  "uheight_ref"
)

init_fn <- function() {
  n_sp <- stan_data$n_species
  n_trt <- stan_data$n_treatments
  n_yr <- stan_data$n_years
  n_cl <- stan_data$n_clusters

  list(
    mu_alpha_psi = 0,
    mu_alpha_phi = 1.5,
    mu_alpha_gam = -1.5,
    mu_alpha_p = 0,

    sigma_alpha_psi = 0.5,
    sigma_alpha_phi = 0.5,
    sigma_alpha_gam = 0.5,
    sigma_alpha_p = 0.5,

    mu_beta_doy_p = 0,
    mu_beta_tod_p = 0,
    sigma_beta_doy_p = 0.25,
    sigma_beta_tod_p = 0.25,

    mu_beta_psi_CLH = 0,
    mu_beta_psi_UP = 0,
    sigma_beta_psi = 0.25,

    alpha_psi_raw = rep(0, n_sp),
    alpha_phi_raw = rep(0, n_sp),
    alpha_gam_raw = rep(0, n_sp),
    alpha_p_raw = rep(0, n_sp),
    beta_doy_p_raw = rep(0, n_sp),
    beta_tod_p_raw = rep(0, n_sp),
    beta_psi_CLH_raw = rep(0, n_sp),
    beta_psi_UP_raw = rep(0, n_sp),

    mu_beta_uheight_psi = 0,
    sigma_beta_uheight_psi = 0.25,
    beta_uheight_psi_raw = rep(0, n_sp),

    mu_beta_tsh_psi_CLH = 0,
    sigma_beta_tsh_psi_CLH = 0.25,
    beta_tsh_psi_CLH_raw = rep(0, n_sp),
    mu_beta_tsh_phi_CLH = 0,
    sigma_beta_tsh_phi_CLH = 0.25,
    beta_tsh_phi_CLH_raw = rep(0, n_sp),
    mu_beta_tsh_gam_CLH = 0,
    sigma_beta_tsh_gam_CLH = 0.25,
    beta_tsh_gam_CLH_raw = rep(0, n_sp),

    mu_beta_tsh_psi_UP = 0,
    sigma_beta_tsh_psi_UP = 0.25,
    beta_tsh_psi_UP_raw = rep(0, n_sp),
    mu_beta_tsh_phi_UP = 0,
    sigma_beta_tsh_phi_UP = 0.25,
    beta_tsh_phi_UP_raw = rep(0, n_sp),
    mu_beta_tsh_gam_UP = 0,
    sigma_beta_tsh_gam_UP = 0.25,
    beta_tsh_gam_UP_raw = rep(0, n_sp),

    mu_community_dte = 0,
    sigma_community_dte = 0.25,
    mu_beta_dte_raw = rep(0, n_sp),
    sigma_dte_trt_dev = 0.25,
    beta_dte_CLH_raw = rep(0, n_sp),
    beta_dte_UP_raw = rep(0, n_sp),

    mu_log_sigma_cluster = -1,
    sigma_log_sigma_cluster = 0.1,
    sigma_cluster_raw = rep(0, n_sp),
    cluster_effect_raw = matrix(0, nrow = n_sp, ncol = n_cl),

    sigma_year_phi = rep(0.1, n_trt),
    sigma_year_gam = rep(0.1, n_trt),
    year_init_phi = rep(0, n_trt),
    year_init_gam = rep(0, n_trt),
    year_raw_phi = matrix(0, nrow = n_trt, ncol = n_yr - 1),
    year_raw_gam = matrix(0, nrow = n_trt, ncol = n_yr - 1),

    sigma_species_year_phi = 0.1,
    sigma_species_year_gam = 0.1,
    species_year_phi_raw = array(0, dim = c(n_sp, n_trt, n_yr)),
    species_year_gam_raw = array(0, dim = c(n_sp, n_trt, n_yr))
  )
}

model <- cmdstan_model(file.path("Code", "Analyses", "dynocc_multispecies.stan"))


fit <- model$sample(
  data = stan_data,
  init = init_fn,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,

)

fit_path <- file.path(fit_dir, "Model_fit.rds")
saveRDS(fit, fit_path)

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

write_csv(diagnostics, file.path(summary_dir, "model_fit_diagnostics.csv"))
