# Build the Stan input for the multi-species dynamic occupancy model.
# Run from the repository root. Reads Data/response.csv and Data/predictors.csv;
# writes Data/model_data.rds and Results/survey_summary.csv.

library(tidyverse)

species <- c("BBWA", "WIWR", "HETH", "TEWA", "BTNW", "BRCR")
treatments <- c("CLH", "OF", "UP")
years <- 2015:2025
n_tsh_grid <- 36  # trajectories start at 0-35 years since harvest

response <- read_csv("Data/response.csv", show_col_types = FALSE)
predictors <- read_csv("Data/predictors.csv", show_col_types = FALSE)

# Sites ----------------------------------------------------------------------
sites <- response %>%
  distinct(location, transect, treatment) %>%
  arrange(location) %>%
  mutate(
    site_id = row_number(),
    cluster_id = match(transect, sort(unique(transect))),
    trt = match(treatment, treatments)
  ) %>%
  left_join(predictors, by = "location")

# Time since harvest, scaled by its 2015 mean and SD across all sites (OF = 0)
tsh <- outer(sites$harvest_year, years, function(h, y) pmax(0, y - h))
tsh[sites$treatment == "OF", ] <- 0
tsh_mean <- mean(tsh[, 1])
tsh_sd <- sd(tsh[, 1])
tsh_grid <- outer(seq_len(n_tsh_grid) - 1, seq_along(years) - 1, `+`)

# Distance to edge, scaled across sites
dte <- as.numeric(scale(sites$distance_to_edge_m))

# Understory height at UP sites: missing values take the UP mean, then the
# covariate is scaled among UP sites and set to zero elsewhere
is_up <- sites$treatment == "UP"
uheight <- sites$uheight
uheight[is_up & is.na(uheight)] <- mean(uheight[is_up], na.rm = TRUE)
uheight_mean <- mean(uheight[is_up])
uheight_sd <- sd(uheight[is_up])
uheight_z <- ifelse(is_up, (uheight - uheight_mean) / uheight_sd, 0)

# Recordings and site-years ------------------------------------------------
recordings <- response %>%
  left_join(select(sites, location, site_id), by = "location") %>%
  arrange(site_id, year, doy, tod)

visits <- recordings %>%
  count(site_id, year, name = "n_recs") %>%
  mutate(first_rec_idx = cumsum(n_recs) - n_recs + 1L)

site_visits <- visits %>%
  mutate(visit_id = row_number()) %>%
  group_by(site_id) %>%
  summarise(n_visits = n(), first_visit_idx = min(visit_id))

# Detection covariates were scaled on the species-long table (one row per
# recording and species), as in the fitted model
scaling <- list(
  doy_mean = mean(recordings$doy),
  doy_sd = sd(rep(recordings$doy, length(species))),
  tod_mean = mean(recordings$tod),
  tod_sd = sd(rep(recordings$tod, length(species))),
  tsh_mean = tsh_mean,
  tsh_sd = tsh_sd,
  uheight_mean = uheight_mean,
  uheight_sd = uheight_sd
)

stan_data <- list(
  n_species = length(species),
  n_sites = nrow(sites),
  n_visits = nrow(visits),
  n_obs = nrow(recordings),
  n_treatments = length(treatments),
  n_years = length(years),
  n_clusters = max(sites$cluster_id),
  treatment = sites$trt,
  cluster = sites$cluster_id,
  n_visits_site = site_visits$n_visits,
  first_visit_idx = site_visits$first_visit_idx,
  tsh = (tsh - tsh_mean) / tsh_sd,
  dte = dte,
  uheight = uheight_z,
  visit_site = visits$site_id,
  visit_year = visits$year - min(years) + 1L,
  n_recs_visit = visits$n_recs,
  first_rec_idx = as.integer(visits$first_rec_idx),
  y = t(as.matrix(recordings[species])) * 1L,
  doy = (recordings$doy - scaling$doy_mean) / scaling$doy_sd,
  tod = (recordings$tod - scaling$tod_mean) / scaling$tod_sd,
  n_tsh_grid = n_tsh_grid,
  tsh_grid = (tsh_grid - tsh_mean) / tsh_sd,
  dte_ref = as.numeric(tapply(dte, sites$trt, median)),
  uheight_ref = median(uheight_z[is_up])
)

saveRDS(
  list(
    stan_data = stan_data,
    species = species,
    treatments = treatments,
    years = years,
    scaling = scaling,
    doy_range = range(recordings$doy),
    tod_range = range(recordings$tod)
  ),
  "Data/model_data.rds"
)

# Survey effort and naive site-year occurrence -------------------------------
dir.create("Results", showWarnings = FALSE)

site_years <- response %>%
  group_by(location, year) %>%
  summarise(n_recordings = n(), across(all_of(species), max), .groups = "drop")

survey_summary <- site_years %>%
  pivot_longer(all_of(species), names_to = "species", values_to = "detected") %>%
  group_by(species) %>%
  summarise(
    site_years_detected = sum(detected),
    site_years = n(),
    percent_detected = 100 * site_years_detected / site_years
  )

survey_effort <- bind_rows(
  count(site_years, value = n_recordings, name = "n") %>%
    mutate(measure = "recordings per site-year"),
  count(count(site_years, location), value = n, name = "n") %>%
    mutate(measure = "years surveyed per site")
) %>%
  select(measure, value, n)

write_csv(survey_summary, "Results/survey_summary.csv")
write_csv(survey_effort, "Results/survey_effort.csv")
