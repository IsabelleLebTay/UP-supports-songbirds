# Understory protection supports forest songbirds but does not offset regional fire-driven losses to old forest occupancy

Code for Lebeuf-Taylor, I., T. Hart, T. Habib, and E. Bayne. *Ecological Applications* (in review).

Data: [Borealis, doi:10.5683/SP3/OZUBY9](https://doi.org/10.5683/SP3/OZUBY9)
Contact: Isabelle Lebeuf-Taylor, University of Alberta

## Overview

Autonomous recording units surveyed 162 sites in 54 clusters, each with one conventional low-retention harvest (CLH), understory protection (UP), and unharvested old forest (OF) site, in northeastern Alberta in 2015, 2016, and 2020–2025. Recordings were processed with HawkEars v1.0.8 and thresholded at a modelled P(true positive) of 0.90 per species. A Bayesian multi-species dynamic occupancy model, fit in Stan, estimates initial occupancy, persistence, colonization, and detection for six species: Bay-breasted Warbler (BBWA), Winter Wren (WIWR), Hermit Thrush (HETH), Tennessee Warbler (TEWA), Black-throated Green Warbler (BTNW), and Brown Creeper (BRCR). The fitted model is then used to project expected occupied area of old spruce and mixedwood forest under three harvest allocation scenarios.

## Data

Download the files from Borealis into `Data/`.

| File | Contents |
| --- | --- |
| `response.csv` | One row per 3-minute recording (3,798 recordings from 480 site-years). |
| `predictors.csv` | One row per site (162). |
| `landbase.csv` | Areas and annual rates used in the habitat-supply projections. |
| `hawkears_thresholds.csv` | Appendix S1: Table S1. HawkEars score thresholds, validated true and false positives, ROC-AUC, precision and recall weighted by confidence-score and recording-origin strata (95% intervals), and detected site-years. Not read by the scripts. |

`response.csv`

| Column | Description |
| --- | --- |
| `location` | Site identifier |
| `transect` | Cluster identifier |
| `treatment` | `CLH`, `UP`, or `OF` |
| `year` | Survey year |
| `doy` | Day of year |
| `tod` | Recording start time, minutes after midnight |
| `BBWA` … `WIWR` | Detection (1) or non-detection (0) |

`predictors.csv`

| Column | Description |
| --- | --- |
| `location` | Site identifier |
| `distance_to_edge_m` | Distance to the nearest edge of the adjoining stand: CLH edge for UP and OF sites, OF edge for CLH sites (m) |
| `harvest_year` | Year of harvest (CLH and UP sites) |
| `uheight` | Height of retained understory spruce (m); used for UP sites only, with two missing values set to the UP mean |

The vegetation inventory, harvest, and wildfire layers behind `landbase.csv` belong to Alberta-Pacific Forest Industries Inc. and are not distributed.

## Analysis

Run the scripts from the repository root, in order.

| Script | Output |
| --- | --- |
| `Analyses/01_prepare_data.R` | Model input; survey effort and naive site-year occurrence |
| `Analyses/02_fit_model.R` | Model fit (4 chains, 1,000 warm-up and 1,000 sampling iterations each) |
| `Analyses/03_summarise_posterior.R` | Convergence, LOO (Table S2), and the posterior summaries reported in the Results |
| `Analyses/04_project_habitat_supply.R` | Habitat-supply projections and Table 3 |
| `Analyses/05_figures.R` | Figures 2–4 and S1–S5 |

The model is in `Analyses/dynocc_multispecies.stan`. Fitting takes roughly an hour.

## Software

R 4.5.3 with tidyverse, cmdstanr, posterior, loo, and patchwork. The published fit was run through cmdstanr with CmdStan 2.38.0 installed, using a model executable compiled under Stan 2.32.1, with seed 2133273333. Other versions or platforms will give slightly different draws.

## License

MIT; see [LICENSE](LICENSE).
