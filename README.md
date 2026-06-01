# Understory protection harvest supports forest songbirds over time but does not replace old forest for sensitive species

## General information
Code for: Understory protection harvest supports forest songbirds over time but does not replace old forest for sensitive species<br>
Dataset archivied on Borealis: https://doi.org/10.5683/SP3/OZUBY9

**Authors:**
- Isabelle Lebeuf-Taylor (corresponding author), Department of Biological Sciences, University of Alberta, Edmonton, Alberta, Canada. lebeufta@ualberta.ca. ORCID: 0000-0001-6809-7249
- Taylor Hart, Department of Biological Sciences, University of Alberta, Edmonton, Alberta, Canada
- Thomas Habib, Alberta-Pacific Forest Industries
- Erin Bayne, Department of Biological Sciences, University of Alberta, Edmonton, Alberta, Canada

**Date of data collection:** 2015-2025

**Geographic location of data collection:** Bird Conservation Region 6, Alberta, Canada

**Funding sources:** Research funding was provided by Alberta Pacific Forest Industries Inc. matched by Mitacs Accelerate (Number IT34177), the Forest Research Improvement Alliance of Alberta (FRIAA)

**Associated publication:** Lebeuf-Taylor I, Hart T, Habib T, Bayne E. Understory protection harvest supports forest songbirds over time but does not replace old forest for sensitive species. *Ecological Applications*. [under review]

**Code repository:** https://github.com/IsabelleLebTay/UP-supports-songbirds.git

## Dataset description
This dataset supports the analysis for “Understory protection harvest supports forest songbirds over time but does not replace old forest for sensitive species.” The study evaluates whether understory protection harvest, a partial harvest method that removes deciduous overstory while retaining younger white spruce, provides long-term conservation value for old-forest-associated boreal songbirds in northeastern Alberta, Canada. Autonomous recording units were deployed from 2015 to 2025 across 162 sites in the Alberta-Pacific Forest Industries Inc. Forest Management Agreement area. Sites were organized into 54 spatial clusters, each containing three treatments: understory protection harvest, conventional low-retention harvest, and unharvested old forest. Most sites were surveyed in three years, with four 3-minute recordings selected per site visit. Species detections were identified using the HawkEars acoustic classifier and species-specific confidence thresholds validated against manual review. The data support a Bayesian multi-species dynamic occupancy model for six focal songbird species: Tennessee Warbler, Hermit Thrush, Winter Wren, Brown Creeper, Bay-breasted Warbler, and Black-throated Green Warbler. The analysis estimates treatment-specific occupancy, colonization, and persistence while accounting for imperfect detection, time of day, day of year, time since harvest, distance to treatment edge, retained understory height, and temporal variation across the 11-year monitoring period. The dataset also supports landscape-level projections of expected occupied habitat for Brown Creeper under alternative understory protection implementation scenarios. The repository includes processed analysis data and derived data products needed to reproduce the analyses reported in the manuscript. Proprietary vegetation inventory and harvest spatial layers owned by Alberta-Pacific Forest Industries Inc. are not included. Analysis code is archived separately on GitHub. (2026-06-01)

## File structure and descriptions

### Data files
| File | Description |
|------|-------------|
| `Data/response.csv` | Recording-level detection/non-detection table for the six focal songbird species. Each row is one selected recording, with the site identifier, transect/cluster, treatment, survey year, day of year, time of day, and binary species detections used in the occupancy model. |
| `Data/predictors.csv` | Site-level covariate table joined to `response.csv` by `location`. Contains the distance-to-edge, harvest-year, and retained-understory-height variables used to prepare the Stan model inputs. |

### Code files
| File | Description |
|------|-------------|
| `Code/Analyses/01_prepare_stan_data.R` | Reads `Data/response.csv` and `Data/predictors.csv`; maps species, sites, clusters, treatments, visits, and recordings; scales day-of-year, time-of-day, time-since-harvest, distance-to-edge, and retained-understory-height covariates; and writes `Data/full_data_for_stan.rds`. |
| `Code/Analyses/02_fit_dynocc_model.R` | Fits the final multi-species dynamic occupancy model with `cmdstanr` and `Code/Analyses/dynocc_multispecies.stan`; writes `Results/Model fits/Model_fit.rds` and `Results/Summaries/model_fit_diagnostics.csv`. Requires `CMDSTAN_PATH` to point to a CmdStan installation. |
| `Code/Analyses/dynocc_multispecies.stan` | Stan model for multi-species dynamic occupancy, including initial occupancy, persistence, colonization, detection, treatment, time-since-harvest, distance-to-edge, retained-understory-height, cluster, and year effects, plus generated quantities used for summaries and contrasts. |
| `Code/Analyses/03_extract_posterior_summaries.R` | Extracts posterior summaries, convergence diagnostics, treatment contrasts, detection summaries, occupancy trajectories, time-since-harvest effects, distance-to-edge effects, and LOO inputs from `Model_fit.rds`; writes CSV outputs to `Results/Summaries`. |
| `Code/Analyses/04_plot_dynocc_figures.R` | Creates manuscript and supplement figures from posterior summaries and generated quantities; writes figure files to `Results/Figures`. |
| `Code/Analyses/05_extract_OF_brcr_trends.R` | Extracts old-forest baseline occupancy, old-forest annual trend summaries, and Brown Creeper background-change scenarios; writes CSV outputs to `Results/Projections`. |
| `Code/Analyses/06_project_habitat_supply.R` | Projects Brown Creeper occupied habitat equivalent under alternative harvest scenarios using the fitted model and landbase constants; writes projection tables to `Results/Projections`. |
| `Code/Analyses/07_plot_habitat_supply.R` | Plots Brown Creeper habitat-supply projections and writes the habitat-supply summary table. |
| `Code/Analyses/08_summarise_naive_occurrence.R` | Summarizes naive detections from `Data/response.csv` by site, treatment, year, and species; writes occurrence summary CSVs to `Results/Summaries`. |
| `Code/Analyses/optional_acoustic_thresholding.Rmd` | Optional provenance notebook for acoustic detector threshold selection and occurrence-table generation. It is not required to rerun the manuscript analysis pipeline from the public `response.csv` and `predictors.csv` files. |

## Variable definitions
### response.csv
| Variable | Description | Units/values |
|----------|-------------|--------------|
| `location` | Site/location identifier for the recording. This is the key used to join site-level covariates from `predictors.csv`. | Character string, e.g., `UP-1-1-CC`. |
| `transect` | Spatial transect or cluster identifier. Sites were organized in treatment triplets within these clusters. | Character string, e.g., `UP-1`. |
| `treatment` | Harvest treatment assigned to the site. | `CLH` = conventional low-retention harvest; `OF` = unharvested old forest; `UP` = understory protection harvest. |
| `year` | Calendar year in which the recording was collected. | Year; observed values in this file are 2015, 2016, and 2020-2025. |
| `doy` | Day of year on which the recording was collected. | Julian day; observed range 146-196. |
| `tod` | Time of day at the start of the recording. | Minutes after midnight; observed range 202-1286. |
| `BBWA` | Bay-breasted Warbler detection/non-detection for the recording after applying the species-specific acoustic confidence threshold. | Binary; `1` = detected, `0` = not detected. |
| `BRCR` | Brown Creeper detection/non-detection for the recording after applying the species-specific acoustic confidence threshold. | Binary; `1` = detected, `0` = not detected. |
| `BTNW` | Black-throated Green Warbler detection/non-detection for the recording after applying the species-specific acoustic confidence threshold. | Binary; `1` = detected, `0` = not detected. |
| `HETH` | Hermit Thrush detection/non-detection for the recording after applying the species-specific acoustic confidence threshold. | Binary; `1` = detected, `0` = not detected. |
| `TEWA` | Tennessee Warbler detection/non-detection for the recording after applying the species-specific acoustic confidence threshold. | Binary; `1` = detected, `0` = not detected. |
| `WIWR` | Winter Wren detection/non-detection for the recording after applying the species-specific acoustic confidence threshold. | Binary; `1` = detected, `0` = not detected. |

### predictors.csv
| Variable | Description | Units/values |
|----------|-------------|--------------|
| `location` | Site/location identifier. This key is used to join site-level covariates to the recording-level rows in `response.csv`. | Character string, e.g., `UP-10-67-UP`. |
| `distance_to_edge_m` | Distance from the recording site to the treatment or stand edge, used as a site-level occupancy covariate. | Metres; observed range 10-454. |
| `harvest_year` | Calendar year in which the site was harvested. Unharvested old-forest sites are coded as missing. This variable is used to calculate time since harvest for harvested treatments. | Year; observed harvested-site range 1985-2015; `NA` = no harvest year. |
| `uheight` | Retained understory height at the site, used by the final model as an understory-protection-only covariate on initial occupancy. | Metres; observed non-missing range 0-17; `NA` = missing value. |

## Software and dependencies

All analyses were performed in R (R Core Team 2024). Key packages:
- `tidyverse` — data manipulation, summaries, and plotting
- `cmdstanr` — Stan model compilation, sampling, and fitted-model access
- `posterior` — posterior draw summaries and diagnostics
- `loo` — leave-one-out cross-validation inputs and model-checking summaries
- `patchwork` — multi-panel figure assembly
- `dotenv` — loading the optional `CMDSTAN_PATH` environment variable from `.env`

The Stan model requires a working CmdStan installation available through `cmdstanr`.



## Reproduction

1. Clone the repository or download from Borealis
2. Place data files in a `Data/` directory relative to the code files
3. Run code files in order

## License

Code in this repository is released under the MIT License; see `LICENSE`.