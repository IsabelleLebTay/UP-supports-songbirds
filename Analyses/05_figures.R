# Manuscript and Appendix S1 figures from the summaries written by scripts 03
# and 04. Writes image files to Results/figures/.

library(tidyverse)
library(patchwork)

dir.create("Results/figures", showWarnings = FALSE)
read <- function(name) read_csv(file.path("Results", paste0(name, ".csv")), show_col_types = FALSE)
save <- function(plot, name, width, height) {
  ggsave(file.path("Results/figures", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}

species_order <- c("TEWA", "WIWR", "HETH", "BRCR", "BTNW", "BBWA")
treatment_colours <- c(CLH = "#E69F00", OF = "#009E73", UP = "#0072B2")
process_labels <- c(psi = "Initial occupancy (ψ)", phi = "Persistence (φ)",
                    gam = "Colonization (γ)")

theme_set(
  theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "bottom"
    )
)

# Figure 2: time-since-harvest slopes ---------------------------------------
fig2 <- read("slopes") %>%
  filter(covariate == "tsh") %>%
  mutate(
    process = factor(process_labels[process], levels = process_labels),
    species = factor(species, levels = rev(species_order)),
    treatment = factor(treatment, levels = c("UP", "CLH")),
    excludes_zero = lower95 > 0 | upper95 < 0
  ) %>%
  ggplot(aes(median, species, colour = treatment)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_linerange(aes(xmin = lower95, xmax = upper95),
                 position = position_dodge(0.6), linewidth = 0.8) +
  geom_point(aes(shape = excludes_zero), position = position_dodge(0.6),
             size = 2, fill = "white", stroke = 0.8) +
  facet_wrap(~process) +
  scale_x_continuous(breaks = -2:2) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  scale_colour_manual(values = treatment_colours, breaks = c("CLH", "UP")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), breaks = c(TRUE, FALSE),
                     labels = c("95% CrI excludes 0", "95% CrI overlaps 0")) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2, override.aes = list(colour = "grey30"))) +
  labs(x = "Time since harvest effect (log-odds)", y = NULL, colour = NULL, shape = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_blank())

save(fig2, "Figure_2", 6, 4)

# Figure 3: occupancy with time since harvest -------------------------------
# UP curves stop at the oldest surveyed UP stand
md <- readRDS("Data/model_data.rds")
visit_tsh <- with(md$stan_data, tsh[cbind(visit_site, visit_year)]) * md$scaling$tsh_sd +
  md$scaling$tsh_mean
max_up_tsh <- max(visit_tsh[md$stan_data$treatment[md$stan_data$visit_site] == 3])

occ_tsh <- read("occupancy_tsh") %>%
  mutate(species = factor(species, levels = species_order))
harvested <- occ_tsh %>%
  filter(treatment != "OF", !(treatment == "UP" & tsh > max_up_tsh)) %>%
  mutate(treatment = factor(treatment, levels = c("OF", "UP", "CLH")))
old_forest <- occ_tsh %>%
  filter(treatment == "OF") %>%
  mutate(treatment = factor(treatment, levels = c("OF", "UP", "CLH")))
y_limits <- tibble(species = factor(species_order, levels = species_order),
                   ymax = rep(c(1, 0.5), each = 3))

fig3 <- ggplot() +
  geom_blank(data = y_limits, aes(x = 0, y = ymax)) +
  geom_blank(data = y_limits, aes(x = 0, y = 0)) +
  geom_hline(data = old_forest, aes(yintercept = median, colour = treatment,
                                    linetype = treatment), linewidth = 0.8) +
  geom_hline(data = old_forest, aes(yintercept = lower80), colour = treatment_colours[["OF"]],
             linetype = "dotted", linewidth = 0.4) +
  geom_hline(data = old_forest, aes(yintercept = upper80), colour = treatment_colours[["OF"]],
             linetype = "dotted", linewidth = 0.4) +
  geom_ribbon(data = harvested, aes(tsh, ymin = lower80, ymax = upper80, fill = treatment),
              alpha = 0.12) +
  geom_line(data = harvested, aes(tsh, median, colour = treatment, linetype = treatment),
            linewidth = 0.8) +
  facet_wrap(~species, scales = "free_y") +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours, guide = "none") +
  scale_linetype_manual(values = c(CLH = "dashed", OF = "solid", UP = "longdash")) +
  labs(x = "Time since harvest (years)", y = "Occupancy probability (year-averaged)",
       colour = "Treatment", linetype = "Treatment") +
  theme(panel.grid.major = element_blank())

save(fig3, "Figure_3", 8, 5)

# Figure 4: proportion of the focal forest occupied -------------------------
supply_species <- c(BRCR = "Brown Creeper", WIWR = "Winter Wren", BBWA = "Bay-breasted Warbler")
scenario_colours <- c("S1: Current conditions" = "#0072B2", "S2: No UP" = "#D55E00",
                      "S3: Maximum feasible UP" = "#009E73")
supply <- read("habitat_supply") %>%
  mutate(scenario = factor(scenario, levels = names(scenario_colours)))

# Panels share the same y-axis span so absolute changes are comparable
y_span <- supply %>%
  group_by(species) %>%
  summarise(span = diff(range(percent_occupied))) %>%
  pull(span) %>%
  max() %>%
  `+`(0.05) %>%
  `*`(10) %>%
  ceiling() / 10

supply_panel <- function(code) {
  d <- filter(supply, species == code)
  top <- max(d$percent_occupied) + 0.025
  ggplot(d, aes(year, percent_occupied, colour = scenario, linetype = scenario)) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(values = scenario_colours,
                        labels = c("S1: Current", "S2: No UP", "S3: Maximum feasible UP")) +
    scale_linetype_manual(values = c("solid", "longdash", "dotdash"),
                          labels = c("S1: Current", "S2: No UP", "S3: Maximum feasible UP")) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 22), expand = expansion(mult = c(0.02, 0.025))) +
    scale_y_continuous(breaks = scales::breaks_width(0.1),
                       labels = scales::label_number(accuracy = 0.1)) +
    coord_cartesian(ylim = c(top - y_span, top), expand = FALSE) +
    labs(x = "Projection year", y = if (code == "BRCR") "Focal landbase occupied (%)" else NULL,
         title = supply_species[[code]], colour = NULL, linetype = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
          legend.key.width = grid::unit(11, "mm"))
}

fig4 <- wrap_plots(map(names(supply_species), supply_panel), nrow = 1, guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Results/figures/Figure_4.png", fig4, width = 180, height = 88, units = "mm",
       dpi = 600, bg = "white")

# Figure S1: initial occupancy, persistence, and colonization ---------------
figS1 <- read("occupancy_processes") %>%
  mutate(
    process = factor(process_labels[process], levels = process_labels),
    species = factor(species, levels = rev(species_order)),
    treatment = factor(treatment, levels = c("OF", "UP", "CLH"))
  ) %>%
  ggplot(aes(median, species, colour = treatment)) +
  geom_linerange(aes(xmin = lower95, xmax = upper95), position = position_dodge(0.6),
                 linewidth = 0.5, alpha = 0.5) +
  geom_linerange(aes(xmin = lower80, xmax = upper80), position = position_dodge(0.6),
                 linewidth = 1.2) +
  geom_point(position = position_dodge(0.6), size = 2) +
  facet_wrap(~process) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_colour_manual(values = treatment_colours) +
  labs(x = "Probability", y = NULL, colour = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major.y = element_blank(),
        panel.spacing.x = unit(1, "lines"))

save(figS1, "Figure_S1", 6, 4)

# Figure S2: detection probability by time of day and day of year -----------
species_colours <- c(BBWA = "#d62728", WIWR = "#e377c2", HETH = "#2ca02c",
                     TEWA = "#1f77b4", BTNW = "#ff7f0e", BRCR = "#9467bd")
detection <- read("detection_curves") %>%
  mutate(species = factor(species, levels = species_order))
detection_panel <- function(cov, xlab) {
  ggplot(filter(detection, covariate == cov), aes(x, median)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = species), alpha = 0.2) +
    geom_line(aes(colour = species), linewidth = 1) +
    scale_colour_manual(values = species_colours) +
    scale_fill_manual(values = species_colours) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = xlab, y = "Detection probability", colour = "Species", fill = "Species") +
    theme(panel.grid = element_blank())
}
figS2 <- detection_panel("tod", "Time of day (minutes since midnight)") +
  detection_panel("doy", "Day of year") +
  plot_layout(guides = "collect") & theme(legend.position = "right")

save(figS2, "Figure_S2", 10, 4)

# Figure S3: understory height (UP sites) ----------------------------------
figS3 <- read("height_curve") %>%
  mutate(species = factor(species, levels = species_order)) %>%
  ggplot(aes(x, median)) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), fill = treatment_colours[["UP"]], alpha = 0.2) +
  geom_line(colour = treatment_colours[["UP"]], linewidth = 0.8) +
  facet_wrap(~species) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(1, 15), breaks = seq(1, 15, 2)) +
  labs(x = "Understory height (m)", y = "Initial occupancy probability") +
  theme(panel.grid = element_blank())

save(figS3, "Figure_S3", 9, 6)

# Figure S4: distance to edge ----------------------------------------------
figS4 <- read("edge_curve") %>%
  mutate(species = factor(species, levels = species_order),
         treatment = factor(treatment, levels = c("OF", "UP", "CLH"))) %>%
  ggplot(aes(x, median, colour = treatment, fill = treatment)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~species) +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Distance to edge (z-scored)", y = "Initial occupancy probability",
       colour = "Treatment", fill = "Treatment") +
  theme(panel.grid = element_blank())

save(figS4, "Figure_S4", 9, 6)

# Figure S5: occupancy through the study period ----------------------------
figS5 <- read("occupancy_year") %>%
  mutate(species = factor(species, levels = species_order),
         treatment = factor(treatment, levels = c("OF", "UP", "CLH"))) %>%
  ggplot(aes(year, median, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), alpha = 0.15, colour = NA,
              position = position_dodge(0.8)) +
  geom_line(aes(linetype = treatment), linewidth = 0.8, position = position_dodge(0.8)) +
  geom_point(size = 1.5, position = position_dodge(0.8)) +
  facet_wrap(~species) +
  scale_colour_manual(values = treatment_colours) +
  scale_fill_manual(values = treatment_colours) +
  scale_linetype_manual(values = c(CLH = "dashed", OF = "solid", UP = "solid")) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(labels = scales::label_number(accuracy = 1, big.mark = "")) +
  labs(x = "Year", y = "Occupancy probability", colour = "Treatment", fill = "Treatment",
       linetype = "Treatment") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

save(figS5, "Figure_S5", 10, 7)
