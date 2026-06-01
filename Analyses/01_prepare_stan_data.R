# Prepare data for the multi-species dynamic occupancy model.
#
# Inputs:
#   Data/response.csv   recording-level detection/non-detection data
#   Data/predictors.csv site-level covariates
#
# Output:
#   Data/full_data_for_stan.rds

library(tidyverse)

species_codes <- c("BBWA", "WIWR", "HETH", "TEWA", "BTNW", "BRCR")
year_base <- 2015L
n_years_total <- 11L
n_tsh_grid <- 36L

treatment_lookup <- c("CLH" = 1L, "OF" = 2L, "UP" = 3L)

add_site_covariates <- function(full_data, predictors) {
  site_covariates <- full_data$site_map %>%
    left_join(predictors, by = "location") %>%
    left_join(
      full_data$visits %>%
        distinct(site_id, treatment_code),
      by = "site_id"
    ) %>%
    arrange(site_id)

  full_data$stan_data$dte <- as.numeric(scale(site_covariates$distance_to_edge_m))

  dte_reference <- tibble(
    treatment = full_data$stan_data$treatment,
    dte = full_data$stan_data$dte
  ) %>%
    group_by(treatment) %>%
    summarise(median_dte = median(dte), .groups = "drop") %>%
    arrange(treatment)

  full_data$stan_data$dte_ref <- dte_reference$median_dte

  up_mean_uheight <- site_covariates %>%
    filter(treatment_code == treatment_lookup[["UP"]], !is.na(uheight)) %>%
    summarise(mean_uheight = mean(uheight)) %>%
    pull(mean_uheight)

  site_covariates <- site_covariates %>%
    mutate(
      uheight = if_else(
        treatment_code == treatment_lookup[["UP"]] & is.na(uheight),
        up_mean_uheight,
        uheight
      )
    )

  up_uheight <- site_covariates %>%
    filter(treatment_code == treatment_lookup[["UP"]]) %>%
    pull(uheight)

  uheight_mean <- mean(up_uheight, na.rm = TRUE)
  uheight_sd <- sd(up_uheight, na.rm = TRUE)

  site_covariates <- site_covariates %>%
    mutate(
      uheight_scaled = (uheight - uheight_mean) / uheight_sd,
      uheight_scaled = if_else(treatment_code == treatment_lookup[["UP"]], uheight_scaled, 0)
    )

  full_data$stan_data$uheight <- site_covariates$uheight_scaled
  full_data$stan_data$uheight_ref <- median(
    site_covariates$uheight_scaled[site_covariates$treatment_code == treatment_lookup[["UP"]]]
  )
  full_data$scaling$uheight_mean <- uheight_mean
  full_data$scaling$uheight_sd <- uheight_sd
  full_data$site_covariates <- site_covariates

  full_data
}

add_tsh_grid <- function(full_data) {
  tsh_grid <- matrix(NA_real_, nrow = n_tsh_grid, ncol = full_data$stan_data$n_years)

  for (g in seq_len(n_tsh_grid)) {
    for (t in seq_len(full_data$stan_data$n_years)) {
      tsh_grid[g, t] <- ((g - 1) + (t - 1) - full_data$scaling$tsh_mean) /
        full_data$scaling$tsh_sd
    }
  }

  full_data$stan_data$n_tsh_grid <- n_tsh_grid
  full_data$stan_data$tsh_grid <- tsh_grid
  full_data
}

prepare_dynocc_data <- function(response, predictors, species_codes) {
  required_response_cols <- c("location", "transect", "treatment", "year", "doy", "tod", species_codes)
  missing_response_cols <- setdiff(required_response_cols, names(response))

  required_predictor_cols <- c("location", "harvest_year", "distance_to_edge_m", "uheight")
  missing_predictor_cols <- setdiff(required_predictor_cols, names(predictors))

  sites <- response %>%
    distinct(location, transect, treatment) %>%
    mutate(
      treatment_code = unname(treatment_lookup[treatment]),
      cluster = transect
    )

  cluster_map <- sites %>%
    distinct(cluster) %>%
    arrange(cluster) %>%
    mutate(cluster_id = row_number())

  site_map <- sites %>%
    distinct(location) %>%
    arrange(location) %>%
    mutate(site_id = row_number())

  sites <- sites %>%
    left_join(cluster_map, by = "cluster") %>%
    left_join(site_map, by = "location")

  det <- response %>%
    left_join(
      sites %>% select(location, site_id, treatment_code, cluster_id),
      by = "location"
    ) %>%
    arrange(site_id, year, doy, tod) %>%
    mutate(recording_id = row_number())

  det_long <- det %>%
    pivot_longer(all_of(species_codes), names_to = "species", values_to = "detected")

  species_map <- tibble(species = species_codes, species_id = seq_along(species_codes))

  det_focal <- det_long %>%
    left_join(species_map, by = "species") %>%
    mutate(year_id = year - year_base + 1L)

  visits <- det_focal %>%
    distinct(location, site_id, year_id, year, treatment_code, cluster_id) %>%
    left_join(predictors %>% select(location, harvest_year), by = "location") %>%
    mutate(
      tsh = if_else(
        is.na(harvest_year) | harvest_year == 0 | treatment_code == treatment_lookup[["OF"]],
        0,
        year - harvest_year
      ),
      is_of = as.integer(treatment_code == treatment_lookup[["OF"]])
    ) %>%
    arrange(site_id, year_id) %>%
    mutate(visit_id = row_number())

  tsh_matrix_raw <- matrix(0, nrow = n_distinct(sites$site_id), ncol = n_years_total)
  site_harvest_info <- visits %>%
    distinct(site_id, treatment_code, harvest_year) %>%
    arrange(site_id)

  for (i in seq_len(nrow(site_harvest_info))) {
    site_id <- site_harvest_info$site_id[[i]]
    harvest_year <- site_harvest_info$harvest_year[[i]]
    treatment_code <- site_harvest_info$treatment_code[[i]]

    if (!is.na(harvest_year) && harvest_year != 0 && treatment_code != treatment_lookup[["OF"]]) {
      for (t in seq_len(n_years_total)) {
        tsh_matrix_raw[site_id, t] <- max(0, (year_base + t - 1) - harvest_year)
      }
    }
  }

  tsh_first_year <- tsh_matrix_raw[, 1]
  tsh_mean <- mean(tsh_first_year)
  tsh_sd <- sd(tsh_first_year)
  tsh_matrix_scaled <- (tsh_matrix_raw - tsh_mean) / tsh_sd

  site_visits <- visits %>%
    group_by(site_id) %>%
    summarise(n_visits = n(), first_visit_idx = min(visit_id), .groups = "drop") %>%
    arrange(site_id)

  doy_mean <- mean(det_focal$doy, na.rm = TRUE)
  doy_sd <- sd(det_focal$doy, na.rm = TRUE)
  tod_mean <- mean(det_focal$tod, na.rm = TRUE)
  tod_sd <- sd(det_focal$tod, na.rm = TRUE)

  det_focal <- det_focal %>%
    mutate(
      doy_scaled = (doy - doy_mean) / doy_sd,
      tod_scaled = (tod - tod_mean) / tod_sd
    )

  rec_id_map <- det_focal %>%
    distinct(site_id, year_id, recording_id, doy, tod, doy_scaled, tod_scaled) %>%
    left_join(visits %>% select(site_id, year_id, visit_id), by = c("site_id", "year_id")) %>%
    arrange(visit_id, recording_id) %>%
    mutate(rec_id = row_number())

  recordings <- det_focal %>%
    left_join(rec_id_map %>% select(recording_id, visit_id, rec_id), by = "recording_id") %>%
    arrange(rec_id, species_id)

  visit_recs <- rec_id_map %>%
    group_by(visit_id) %>%
    summarise(n_recs = n(), first_rec_idx = min(rec_id), .groups = "drop")

  visits <- visits %>%
    left_join(visit_recs, by = "visit_id")

  y_matrix <- matrix(0L, nrow = length(species_codes), ncol = nrow(rec_id_map))
  for (i in seq_len(nrow(recordings))) {
    y_matrix[recordings$species_id[[i]], recordings$rec_id[[i]]] <- as.integer(recordings$detected[[i]])
  }

  stan_data <- list(
    n_species = length(species_codes),
    n_sites = n_distinct(sites$site_id),
    n_years = n_years_total,
    n_treatments = 3L,
    n_clusters = n_distinct(sites$cluster_id),
    n_visits = nrow(visits),
    n_obs = nrow(rec_id_map),
    treatment = sites %>% arrange(site_id) %>% pull(treatment_code),
    cluster = sites %>% arrange(site_id) %>% pull(cluster_id),
    n_visits_site = site_visits$n_visits,
    first_visit_idx = site_visits$first_visit_idx,
    tsh = tsh_matrix_scaled,
    visit_site = visits$site_id,
    visit_year = visits$year_id,
    visit_tsh = visits$tsh,
    visit_is_of = visits$is_of,
    n_recs_visit = visits$n_recs,
    first_rec_idx = visits$first_rec_idx,
    y = y_matrix,
    doy = rec_id_map %>% arrange(rec_id) %>% pull(doy_scaled),
    tod = rec_id_map %>% arrange(rec_id) %>% pull(tod_scaled)
  )

  full_data <- list(
    stan_data = stan_data,
    species_map = species_map,
    site_map = site_map,
    cluster_map = cluster_map,
    visits = visits,
    recordings = recordings,
    scaling = list(
      doy_mean = doy_mean,
      doy_sd = doy_sd,
      tod_mean = tod_mean,
      tod_sd = tod_sd,
      tsh_mean = tsh_mean,
      tsh_sd = tsh_sd
    ),
    year_base = year_base
  )

  full_data %>%
    add_tsh_grid() %>%
    add_site_covariates(predictors)
}

response <- read_csv(file.path("Data", "response.csv"), show_col_types = FALSE)
predictors <- read_csv(file.path("Data", "predictors.csv"), show_col_types = FALSE)

full_data <- prepare_dynocc_data(response, predictors, species_codes)

saveRDS(full_data, file.path("Data", "full_data_for_stan.rds"))
