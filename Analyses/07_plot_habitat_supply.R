# Optional visualization of the three-species projections from script 06.
# Endpoint tables are written by script 06, so plotting is not required.
library(tidyverse)

figure_dir <- file.path("Results", "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
projections <- read_csv(file.path("Results", "Projections", "habitat_supply_projections.csv"),
                        show_col_types = FALSE)
stopifnot(all(c("BRCR", "WIWR", "BBWA") %in% projections$species),
          all(projections$focal_landbase_ha > 0))
scenario_colours <- c("S1: Current conditions" = "#0072B2",
                      "S2: No UP" = "#D55E00",
                      "S3: Maximum feasible UP" = "#009E73")
species_labels <- c(BRCR = "Brown Creeper", WIWR = "Winter Wren",
                    BBWA = "Bay-breasted Warbler")
projections <- projections %>%
  mutate(species = factor(species, levels = names(species_labels)),
         scenario = factor(scenario, levels = names(scenario_colours)))

# The denominator declines with OF area; it is not the starting landbase.
figure_4 <- ggplot(projections, aes(proj_year, occupied_pct, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~species, scales = "free_y", ncol = 3,
             labeller = labeller(species = species_labels)) +
  scale_colour_manual(values = scenario_colours) +
  labs(x = "Projection year", y = "Remaining focal forest occupied (%)", colour = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "Figure_4_habitat_supply_proportion.jpg"),
       figure_4, width = 10, height = 4, dpi = 300)

# Companion uncertainty plot; no additional supplement figure number is implied.
occupied_area <- ggplot(projections, aes(proj_year, median_ohe, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~species, scales = "free_y", ncol = 3,
             labeller = labeller(species = species_labels)) +
  scale_colour_manual(values = scenario_colours) +
  scale_fill_manual(values = scenario_colours) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Projection year", y = "Expected occupied area (ha)", colour = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figure_dir, "habitat_supply_occupied_area.jpeg"),
       occupied_area, width = 10, height = 4, dpi = 300)
