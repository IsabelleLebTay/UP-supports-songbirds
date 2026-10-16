# Posterior summaries reported in the manuscript and Appendix S1.
# Reads Results/model_fit.rds and Data/model_data.rds; writes CSVs to Results/.

library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

fit <- readRDS("Results/model_fit.rds")
md <- readRDS("Data/model_data.rds")

species <- md$species
treatments <- md$treatments
n_sp <- length(species)
n_yr <- length(md$years)
n_sites <- md$stan_data$n_sites

# Median, 80% and 95% intervals, and posterior probability of a positive value
# for each column of a draws matrix
summarise_cols <- function(m) {
  m <- as.matrix(m)
  q <- apply(m, 2, quantile, probs = c(0.5, 0.025, 0.1, 0.9, 0.975), names = FALSE)
  tibble(
    median = q[1, ], lower95 = q[2, ], lower80 = q[3, ], upper80 = q[4, ],
    upper95 = q[5, ], p_positive = colMeans(m > 0)
  )
}

draws <- as_draws_matrix(fit$draws(c(
  "alpha_psi", "alpha_p", "beta_trt_psi", "beta_doy_p", "beta_tod_p",
  "beta_tsh_psi", "beta_tsh_phi", "beta_tsh_gam", "beta_dte", "beta_uheight_psi",
  "p", "psi_baseline", "phi_by_trt_year", "gam_by_trt_year", "occ_trajectory"
)))
d <- function(...) as.numeric(draws[, sprintf(...)])

# Convergence --------------------------------------------------------------
diag_all <- summarise_draws(fit$draws(), "rhat", "ess_bulk", "ess_tail")
sampler <- fit$diagnostic_summary()

convergence <- tibble(
  max_rhat = max(diag_all$rhat, na.rm = TRUE),
  min_ess_bulk = min(diag_all$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(diag_all$ess_tail, na.rm = TRUE),
  divergences = sum(sampler$num_divergent),
  max_treedepth_hits = sum(sampler$num_max_treedepth),
  min_ebfmi = min(sampler$ebfmi)
)

# Leave-one-out cross-validation by species (Table S2) ---------------------
loo_table <- map_dfr(seq_len(n_sp), function(i) {
  ll <- fit$draws(sprintf("log_lik[%d,%d]", i, seq_len(n_sites)), format = "draws_array")
  res <- loo(ll, r_eff = relative_eff(exp(ll)))
  k <- res$diagnostics$pareto_k
  tibble(
    species = species[i],
    elpd_loo = res$estimates["elpd_loo", "Estimate"],
    se_elpd_loo = res$estimates["elpd_loo", "SE"],
    p_loo = res$estimates["p_loo", "Estimate"],
    n_k_high = sum(k >= 0.7 & k < 1),
    n_k_very_high = sum(k >= 1),
    max_k = max(k)
  )
})

# Detection probability ----------------------------------------------------
detection <- tibble(species = species) %>%
  bind_cols(summarise_cols(draws[, sprintf("p[%d]", seq_len(n_sp))]))

# Initial occupancy and year-averaged persistence and colonization (Fig S1)
process_grid <- expand_grid(sp = seq_len(n_sp), trt = 1:3, process = c("psi", "phi", "gam"))
process_draws <- pmap(process_grid, function(sp, trt, process) {
  if (process == "psi") return(d("psi_baseline[%d,%d]", sp, trt))
  rowMeans(sapply(2:n_yr, function(t) d("%s_by_trt_year[%d,%d,%d]", process, sp, trt, t)))
})
occupancy_processes <- process_grid %>%
  transmute(species = species[sp], treatment = treatments[trt], process) %>%
  bind_cols(summarise_cols(do.call(cbind, process_draws)))

# Covariate slopes (Fig 2 and Results) -------------------------------------
slope_vars <- bind_rows(
  expand_grid(sp = seq_len(n_sp), trt = c(1, 3), process = c("psi", "phi", "gam")) %>%
    mutate(covariate = "tsh", variable = sprintf("beta_tsh_%s[%d,%d]", process, sp, trt)),
  expand_grid(sp = seq_len(n_sp), trt = 1:3) %>%
    mutate(covariate = "dte", process = "psi, phi, gam",
           variable = sprintf("beta_dte[%d,%d]", sp, trt)),
  tibble(sp = seq_len(n_sp), trt = 3, covariate = "uheight", process = "psi",
         variable = sprintf("beta_uheight_psi[%d]", sp))
)
slopes <- slope_vars %>%
  transmute(covariate, process, species = species[sp], treatment = treatments[trt]) %>%
  bind_cols(summarise_cols(draws[, slope_vars$variable]))

# Occupancy with time since harvest, averaged over 2015-2025 (Fig 3) -------
traj_mean <- function(sp, trt, g) {
  rowMeans(sapply(seq_len(n_yr), function(t) d("occ_trajectory[%d,%d,%d,%d]", sp, trt, g, t)))
}
tsh_grid <- bind_rows(
  expand_grid(sp = seq_len(n_sp), trt = c(1, 3), g = seq_len(md$stan_data$n_tsh_grid)),
  expand_grid(sp = seq_len(n_sp), trt = 2, g = 1)
)
occupancy_tsh <- tsh_grid %>%
  transmute(species = species[sp], treatment = treatments[trt], tsh = g - 1) %>%
  bind_cols(summarise_cols(do.call(cbind, pmap(tsh_grid, traj_mean))))

# Occupancy through 2015-2025 from initial occupancy, persistence, and
# colonization at mean covariates (Fig S5) -----------------------------------
year_grid <- expand_grid(sp = seq_len(n_sp), trt = 1:3)
year_draws <- pmap(year_grid, function(sp, trt) {
  occ <- matrix(NA_real_, nrow(draws), n_yr)
  occ[, 1] <- d("psi_baseline[%d,%d]", sp, trt)
  for (t in 2:n_yr) {
    occ[, t] <- occ[, t - 1] * d("phi_by_trt_year[%d,%d,%d]", sp, trt, t) +
      (1 - occ[, t - 1]) * d("gam_by_trt_year[%d,%d,%d]", sp, trt, t)
  }
  occ
})
occupancy_year <- year_grid %>%
  transmute(species = species[sp], treatment = treatments[trt]) %>%
  slice(rep(seq_len(n()), each = n_yr)) %>%
  mutate(year = rep(md$years, nrow(year_grid))) %>%
  bind_cols(summarise_cols(do.call(cbind, year_draws)))

occupancy_change <- year_grid %>%
  transmute(species = species[sp], treatment = treatments[trt]) %>%
  mutate(
    first_year = map_dbl(year_draws, ~ median(.x[, 1])),
    last_year = map_dbl(year_draws, ~ median(.x[, n_yr])),
    p_decrease = map_dbl(year_draws, ~ mean(.x[, n_yr] < .x[, 1]))
  )

# Marginal effects (Figs S2-S4) --------------------------------------------
curve <- function(grid, f) {
  map_dfr(seq_len(n_sp), function(sp) {
    tibble(species = species[sp], x = grid) %>%
      bind_cols(summarise_cols(sapply(grid, function(x) f(sp, x))))
  })
}

detection_curves <- bind_rows(
  tod = curve(seq(md$tod_range[1], md$tod_range[2], length.out = 100), function(sp, x) {
    plogis(d("alpha_p[%d]", sp) + d("beta_tod_p[%d]", sp) *
             (x - md$scaling$tod_mean) / md$scaling$tod_sd)
  }),
  doy = curve(seq(md$doy_range[1], md$doy_range[2], length.out = 100), function(sp, x) {
    plogis(d("alpha_p[%d]", sp) + d("beta_doy_p[%d]", sp) *
             (x - md$scaling$doy_mean) / md$scaling$doy_sd)
  }),
  .id = "covariate"
)

# Understory height at UP sites, from 1.5 m to the tallest observed value
up_heights <- md$stan_data$uheight[md$stan_data$treatment == 3] *
  md$scaling$uheight_sd + md$scaling$uheight_mean
height_curve <- curve(seq(1.5, max(up_heights), length.out = 50), function(sp, x) {
  plogis(d("alpha_psi[%d]", sp) + d("beta_trt_psi[%d,3]", sp) + d("beta_uheight_psi[%d]", sp) *
           (x - md$scaling$uheight_mean) / md$scaling$uheight_sd)
})

edge_curve <- map_dfr(1:3, function(trt) {
  curve(seq(-2, 2, length.out = 50), function(sp, x) {
    plogis(d("alpha_psi[%d]", sp) + d("beta_trt_psi[%d,%d]", sp, trt) +
             d("beta_dte[%d,%d]", sp, trt) * x)
  }) %>%
    mutate(treatment = treatments[trt], .after = species)
})

# Write --------------------------------------------------------------------
outputs <- list(
  convergence = convergence, loo = loo_table, detection = detection,
  occupancy_processes = occupancy_processes, slopes = slopes,
  occupancy_tsh = occupancy_tsh, occupancy_year = occupancy_year,
  occupancy_change = occupancy_change, detection_curves = detection_curves,
  height_curve = height_curve, edge_curve = edge_curve
)
iwalk(outputs, ~ write_csv(.x, file.path("Results", paste0(.y, ".csv"))))
