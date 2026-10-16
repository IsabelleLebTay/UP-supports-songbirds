# Fit the multi-species dynamic occupancy model.
# Reads Data/model_data.rds; writes Results/model_fit.rds.

library(cmdstanr)

stan_data <- readRDS("Data/model_data.rds")$stan_data

n_sp <- stan_data$n_species
n_trt <- stan_data$n_treatments
n_yr <- stan_data$n_years
n_cl <- stan_data$n_clusters

# Non-centred deviates start at zero, location hyperparameters at their prior
# means, and standard deviations at small values
init_fn <- function() {
  zero <- rep(0, n_sp)
  list(
    mu_alpha_psi = 0, mu_alpha_phi = 1.5, mu_alpha_gam = -1.5, mu_alpha_p = 0,
    sigma_alpha_psi = 0.5, sigma_alpha_phi = 0.5, sigma_alpha_gam = 0.5, sigma_alpha_p = 0.5,
    alpha_psi_raw = zero, alpha_phi_raw = zero, alpha_gam_raw = zero, alpha_p_raw = zero,

    mu_beta_doy_p = 0, mu_beta_tod_p = 0,
    sigma_beta_doy_p = 0.25, sigma_beta_tod_p = 0.25,
    beta_doy_p_raw = zero, beta_tod_p_raw = zero,

    mu_beta_psi_CLH = 0, mu_beta_psi_UP = 0, sigma_beta_psi = 0.25,
    beta_psi_CLH_raw = zero, beta_psi_UP_raw = zero,

    mu_beta_uheight_psi = 0, sigma_beta_uheight_psi = 0.25, beta_uheight_psi_raw = zero,

    mu_beta_tsh_psi_CLH = 0, sigma_beta_tsh_psi_CLH = 0.25, beta_tsh_psi_CLH_raw = zero,
    mu_beta_tsh_phi_CLH = 0, sigma_beta_tsh_phi_CLH = 0.25, beta_tsh_phi_CLH_raw = zero,
    mu_beta_tsh_gam_CLH = 0, sigma_beta_tsh_gam_CLH = 0.25, beta_tsh_gam_CLH_raw = zero,
    mu_beta_tsh_psi_UP = 0, sigma_beta_tsh_psi_UP = 0.25, beta_tsh_psi_UP_raw = zero,
    mu_beta_tsh_phi_UP = 0, sigma_beta_tsh_phi_UP = 0.25, beta_tsh_phi_UP_raw = zero,
    mu_beta_tsh_gam_UP = 0, sigma_beta_tsh_gam_UP = 0.25, beta_tsh_gam_UP_raw = zero,

    mu_community_dte = 0, sigma_community_dte = 0.25, mu_beta_dte_raw = zero,
    sigma_dte_trt_dev = 0.25, beta_dte_CLH_raw = zero, beta_dte_UP_raw = zero,

    mu_log_sigma_cluster = -1, sigma_log_sigma_cluster = 0.1,
    sigma_cluster_raw = zero, cluster_effect_raw = matrix(0, n_sp, n_cl),

    sigma_year_phi = rep(0.1, n_trt), sigma_year_gam = rep(0.1, n_trt),
    year_init_phi = rep(0, n_trt), year_init_gam = rep(0, n_trt),
    year_raw_phi = matrix(0, n_trt, n_yr - 1), year_raw_gam = matrix(0, n_trt, n_yr - 1),
    sigma_species_year_phi = 0.1, sigma_species_year_gam = 0.1,
    species_year_phi_raw = array(0, c(n_sp, n_trt, n_yr)),
    species_year_gam_raw = array(0, c(n_sp, n_trt, n_yr))
  )
}

model <- cmdstan_model("Analyses/dynocc_multispecies.stan")

fit <- model$sample(
  data = stan_data,
  init = init_fn,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 2133273333
)

dir.create("Results", showWarnings = FALSE)
fit$save_object("Results/model_fit.rds")
