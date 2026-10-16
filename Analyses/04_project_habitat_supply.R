# Regional projections of expected occupied area for Brown Creeper, Winter Wren,
# and Bay-breasted Warbler under three harvest allocation scenarios (Table 3).
# Reads Results/model_fit.rds, Data/model_data.rds, and Data/landbase.csv;
# writes Results/habitat_supply.csv and Results/table3.csv.

library(tidyverse)
library(cmdstanr)
library(posterior)

fit <- readRDS("Results/model_fit.rds")
md <- readRDS("Data/model_data.rds")
landbase <- deframe(select(read_csv("Data/landbase.csv", show_col_types = FALSE), component, value))

focal_species <- c("BRCR", "WIWR", "BBWA")
n_yr <- length(md$years)
horizon <- 22  # projection years; also the number of harvest cohorts (ages 0-21)
ages <- 0:(horizon - 1)

of_initial <- landbase[["reserved_OF"]] + landbase[["harvestable_OF"]]
harvested <- landbase[["CLH"]] + landbase[["UP"]]
of_loss <- landbase[["annual_fire_loss"]] + landbase[["annual_harvest_loss"]] -
  landbase[["annual_inflow"]]
new_harvest <- harvested / horizon

scenarios <- c(
  "S1: Current conditions" = landbase[["UP"]] / harvested,
  "S2: No UP" = 0,
  "S3: Maximum feasible UP" = 0.2
)

draws <- as_draws_matrix(fit$draws(c(
  "alpha_psi", "alpha_phi", "alpha_gam", "beta_trt_psi", "beta_tsh_psi",
  "beta_tsh_phi", "beta_tsh_gam", "beta_dte", "beta_uheight_psi",
  "species_year_effect_phi", "species_year_effect_gam"
)))
d <- function(...) as.numeric(draws[, sprintf(...)])

# Mean annual occupancy over an 11-year trajectory starting at each harvest age,
# at treatment-specific median distance to edge, median UP understory height,
# and zero cluster effect. Calendar-year variation is removed by holding each
# random walk at its 2016-2025 mean.
cohort_occupancy <- function(sp) {
  walk_level <- function(process, trt) {
    rowMeans(sapply(2:n_yr, function(t) d("species_year_effect_%s[%d,%d,%d]", process, sp, trt, t)))
  }

  occ <- array(NA_real_, c(nrow(draws), 3, horizon))
  for (trt in 1:3) {
    dte <- d("beta_dte[%d,%d]", sp, trt) * md$stan_data$dte_ref[trt]
    height <- if (trt == 3) d("beta_uheight_psi[%d]", sp) * md$stan_data$uheight_ref else 0
    phi0 <- d("alpha_phi[%d]", sp) + dte + walk_level("phi", trt)
    gam0 <- d("alpha_gam[%d]", sp) + dte + walk_level("gam", trt)

    for (a in seq_along(ages)) {
      tsh <- (ages[a] + seq_len(n_yr) - 1 - md$scaling$tsh_mean) / md$scaling$tsh_sd
      psi <- plogis(d("alpha_psi[%d]", sp) + d("beta_trt_psi[%d,%d]", sp, trt) +
                      d("beta_tsh_psi[%d,%d]", sp, trt) * tsh[1] + dte + height)
      total <- psi
      for (t in 2:n_yr) {
        phi <- plogis(phi0 + d("beta_tsh_phi[%d,%d]", sp, trt) * tsh[t])
        gam <- plogis(gam0 + d("beta_tsh_gam[%d,%d]", sp, trt) * tsh[t])
        psi <- psi * phi + (1 - psi) * gam
        total <- total + psi
      }
      occ[, trt, a] <- total / n_yr
    }
  }
  occ
}

# Harvested area is held in 22 equal annual cohorts. Each year the oldest
# cohort leaves, a new cohort is allocated to CLH or UP, and old forest
# declines by the net annual loss.
project <- function(occ, up_share) {
  clh <- rep(landbase[["CLH"]] / horizon, horizon)
  up <- rep(landbase[["UP"]] / horizon, horizon)
  of <- of_initial
  occupied <- function() of * occ[, 2, 1] + occ[, 1, ] %*% clh + occ[, 3, ] %*% up

  out <- matrix(NA_real_, dim(occ)[1], horizon + 1)
  out[, 1] <- occupied()
  for (year in seq_len(horizon)) {
    of <- of - of_loss
    clh <- c(new_harvest * (1 - up_share), clh[-horizon])
    up <- c(new_harvest * up_share, up[-horizon])
    out[, year + 1] <- occupied()
  }
  out
}

habitat_supply <- map_dfr(focal_species, function(code) {
  occ <- cohort_occupancy(match(code, md$species))
  imap_dfr(scenarios, function(up_share, scenario) {
    area <- project(occ, up_share)
    q <- apply(area, 2, quantile, probs = c(0.5, 0.025, 0.1, 0.9, 0.975), names = FALSE)
    tibble(
      species = code, scenario = scenario, year = 0:horizon,
      median = q[1, ], lower95 = q[2, ], lower80 = q[3, ], upper80 = q[4, ], upper95 = q[5, ]
    )
  })
}) %>%
  mutate(
    focal_area = of_initial + harvested - of_loss * year,
    percent_occupied = 100 * median / focal_area
  )

# Table 3 reports differences between posterior medians at years 0 and 22
table3 <- habitat_supply %>%
  filter(year %in% c(0, horizon)) %>%
  select(species, scenario, year, median) %>%
  pivot_wider(names_from = year, values_from = median, names_prefix = "year_") %>%
  mutate(
    change = .data[[paste0("year_", horizon)]] - year_0,
    percent_change = 100 * change / year_0
  )

write_csv(habitat_supply, "Results/habitat_supply.csv")
write_csv(table3, "Results/table3.csv")
