# Plot BRCR habitat-supply projections.
#
# Inputs:
#   Results/Projections/habitat_supply_projections.csv
#
# Outputs:
#   Results/Figures/Figure_4_BRCR_habitat_supply_proportion.jpg
#   Results/Figures/Figure_S6_BRCR_habitat_supply_ohe.jpeg
#   Results/Projections/habitat_supply_summary_table.csv

library(tidyverse)

projection_dir <- file.path("Results", "Projections")
figure_dir <- file.path("Results", "Figures")
dir.create(projection_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

projections <- read_csv(
  file.path("Results", "Projections", "habitat_supply_projections.csv"),
  show_col_types = FALSE
)

scenario_colours <- c(
  "S1: Current conditions" = "#0072B2",
  "S2: No UP" = "#D55E00",
  "S3: Maximum feasible UP" = "#009E73"
)

scenario_labels <- c(
  "S1: Current conditions" = "S1: Current",
  "S2: No UP" = "S2: No UP",
  "S3: Maximum feasible UP" = "S3: Max UP"
)

landbase <- list(
  of_ha_initial = 535912,
  clh_ha_initial = 34211,
  up_ha_initial = 5252,
  of_annual_loss = 786,
  max_tsh = 22L
)

harvest_ha <- landbase$clh_ha_initial + landbase$up_ha_initial

brcr_proj <- projections %>%
  filter(species == "BRCR") %>%
  mutate(
    focal_landbase_ha = (landbase$of_ha_initial - landbase$of_annual_loss * proj_year) + harvest_ha,
    prop_occupied = median_ohe / focal_landbase_ha,
    scenario = factor(scenario, levels = names(scenario_colours))
  )

theme_projection <- theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black"),
    plot.margin = margin(5.5, 70, 5.5, 5.5)
  )

end_labels <- brcr_proj %>%
  filter(proj_year == max(proj_year)) %>%
  mutate(label = scenario_labels[as.character(scenario)])

figure_4 <- ggplot(brcr_proj, aes(x = proj_year, y = prop_occupied, colour = scenario)) +
  geom_line(linewidth = 1.2) +
  geom_text(
    data = end_labels,
    aes(x = proj_year + 0.6, label = label),
    hjust = 0,
    vjust = 0.5,
    size = 4.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = scenario_colours) +
  scale_x_continuous(
    breaks = seq(0, landbase$max_tsh, by = 5),
    limits = c(0, landbase$max_tsh + 8),
    expand = expansion(mult = c(0.01, 0))
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Projection year",
    y = "Focal habitat occupied (%)"
  ) +
  theme_projection

ggsave(
  file.path(figure_dir, "Figure_4_BRCR_habitat_supply_proportion.jpg"),
  figure_4,
  width = 6,
  height = 5,
  dpi = 300
)

figure_s6 <- ggplot(brcr_proj, aes(x = proj_year, y = median_ohe, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = scenario_colours, labels = scenario_labels) +
  scale_fill_manual(values = scenario_colours, labels = scenario_labels) +
  scale_x_continuous(breaks = seq(0, landbase$max_tsh, by = 5)) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Projection year",
    y = "Occupied habitat equivalents"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(figure_dir, "Figure_S6_BRCR_habitat_supply_ohe.jpeg"),
  figure_s6,
  width = 7,
  height = 5,
  dpi = 300
)

summary_table <- projections %>%
  filter(proj_year %in% c(0, landbase$max_tsh)) %>%
  select(scenario, species, proj_year, median_ohe) %>%
  pivot_wider(names_from = proj_year, values_from = median_ohe, names_prefix = "year_") %>%
  mutate(
    change_ohe = .data[[paste0("year_", landbase$max_tsh)]] - year_0,
    change_pct = 100 * change_ohe / year_0
  )

write_csv(summary_table, file.path(projection_dir, "habitat_supply_summary_table.csv"))