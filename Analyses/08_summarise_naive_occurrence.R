# Summarise naive detections in the public response data.
#
# Inputs:
#   Data/response.csv
#
# Outputs:
#   Results/Summaries/occurrence_summary_*.csv

library(tidyverse)

summary_dir <- file.path("Results", "Summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

species_cols <- c("BBWA", "BRCR", "BTNW", "HETH", "TEWA", "WIWR")

response <- read_csv(file.path("Data", "response.csv"), show_col_types = FALSE)

occ_long <- response %>%
  pivot_longer(all_of(species_cols), names_to = "species", values_to = "detected") %>%
  mutate(detected = as.integer(detected > 0))

site_summary <- occ_long %>%
  group_by(species, location) %>%
  summarise(ever_detected = as.integer(any(detected == 1)), .groups = "drop") %>%
  group_by(species) %>%
  summarise(
    n_sites_detected = sum(ever_detected),
    total_sites = n(),
    pct_sites = 100 * n_sites_detected / total_sites,
    .groups = "drop"
  )

site_treatment_summary <- occ_long %>%
  group_by(species, treatment, location) %>%
  summarise(ever_detected = as.integer(any(detected == 1)), .groups = "drop") %>%
  group_by(species, treatment) %>%
  summarise(
    n_sites_detected = sum(ever_detected),
    total_sites = n(),
    pct_sites = 100 * n_sites_detected / total_sites,
    .groups = "drop"
  )

year_summary <- occ_long %>%
  group_by(species, year) %>%
  summarise(
    n_recordings = n(),
    n_recordings_detected = sum(detected),
    n_sites_with_data = n_distinct(location),
    n_sites_detected = n_distinct(location[detected == 1]),
    recording_detection_rate = n_recordings_detected / n_recordings,
    site_detection_rate = n_sites_detected / n_sites_with_data,
    .groups = "drop"
  )

year_treatment_summary <- occ_long %>%
  group_by(species, year, treatment) %>%
  summarise(
    n_recordings = n(),
    n_recordings_detected = sum(detected),
    n_sites_with_data = n_distinct(location),
    n_sites_detected = n_distinct(location[detected == 1]),
    recording_detection_rate = n_recordings_detected / n_recordings,
    site_detection_rate = n_sites_detected / n_sites_with_data,
    .groups = "drop"
  )

site_year_matrix <- occ_long %>%
  group_by(species, location, treatment, year) %>%
  summarise(detected = as.integer(any(detected == 1)), .groups = "drop")

year_breadth <- year_summary %>%
  group_by(species) %>%
  summarise(
    total_years = n(),
    years_detected = sum(n_sites_detected > 0),
    pct_years = 100 * years_detected / total_years,
    years_present = paste(year[n_sites_detected > 0], collapse = ", "),
    .groups = "drop"
  )

overall <- site_summary %>%
  left_join(year_breadth, by = "species") %>%
  select(
    species,
    n_sites_detected,
    total_sites,
    pct_sites,
    years_detected,
    total_years,
    pct_years,
    years_present
  )

write_csv(site_summary, file.path(summary_dir, "occurrence_summary_by_site.csv"))
write_csv(site_treatment_summary, file.path(summary_dir, "occurrence_summary_by_site_treatment.csv"))
write_csv(year_summary, file.path(summary_dir, "occurrence_summary_by_year.csv"))
write_csv(year_treatment_summary, file.path(summary_dir, "occurrence_summary_by_year_treatment.csv"))
write_csv(site_year_matrix, file.path(summary_dir, "occurrence_summary_site_year_matrix.csv"))
write_csv(year_breadth, file.path(summary_dir, "occurrence_summary_year_breadth.csv"))
write_csv(overall, file.path(summary_dir, "occurrence_summary_overall.csv"))