# Understory protection harvest supports forest songbirds over time but does not offset regional fire-driven losses of old-forest songbird occupancy

Analysis code for the revised manuscript by Isabelle Lebeuf-Taylor, Taylor Hart, Thomas Habib, and Erin Bayne (*Ecological Applications*, under review).

- **Data:** [Borealis](https://doi.org/10.5683/SP3/OZUBY9)
- **Code:** [UP-supports-songbirds](https://github.com/IsabelleLebTay/UP-supports-songbirds)
- **Corresponding author:** Isabelle Lebeuf-Taylor, Department of Biological Sciences, University of Alberta; lebeufta@ualberta.ca; [ORCID](https://orcid.org/0000-0001-6809-7249)
- **Affiliations:** Lebeuf-Taylor, Hart, and Bayne: University of Alberta; Habib: Alberta-Pacific Forest Industries.
- **Funding:** Alberta-Pacific Forest Industries Inc., Mitacs Accelerate (IT34177), and FRIAA.

## Study and analysis

The study covers 2015–2025 in northeastern Alberta, Canada (Bird Conservation Region 6). It compares understory protection harvest (UP), conventional low-retention harvest (CLH), and unharvested old forest (OF), with one site per treatment in each of 54 spatial clusters.

The revised analysis uses 3,798 three-minute recordings from 162 sites and 480 site-years. There are 466 site-years with eight recordings and 14 with fewer (one with two, six with four, one with five, three with six, and three with seven). Eleven sites were surveyed in two years, 146 in three years, and five in four years. Survey years present are 2015, 2016, and 2020–2025; the model propagates occupancy through unsurveyed years.

The six focal species are Bay-breasted Warbler (BBWA), Brown Creeper (BRCR), Black-throated Green Warbler (BTNW), Hermit Thrush (HETH), Tennessee Warbler (TEWA), and Winter Wren (WIWR). Detections use species-specific HawkEars thresholds at modelled P(true positive) = 0.90. A recording is positive if any scored window reaches its species' threshold.

The Bayesian multi-species dynamic occupancy model estimates initial occupancy, persistence, colonization, and detection. Initial occupancy includes treatment, time since harvest (TSH), distance to edge (DTE), UP understory height, and cluster effects. Persistence and colonization include treatment-specific TSH and DTE slopes, cluster effects, and treatment-level temporal random walks with species deviations. Detection includes species-specific time-of-day and day-of-year effects. OF TSH slopes are fixed at zero. The likelihood integrates latent occupancy states and compounds annual transitions across gaps between surveys.

## Revised data package

Download the revised CSVs from Borealis and place them in **Data/REVISED/** at the repository root. The code deliberately uses this directory, so the original CSVs in Data/ can be retained. Earlier four-recording data releases cannot reproduce the revised analysis. Data and generated results are excluded from Git: GitHub distributes code, and Borealis distributes data.

| File in Data/REVISED | Purpose |
| --- | --- |
| response.csv | 3,798 recording-level rows of survey covariates and binary detections; required for scripts 01 and 08. |
| predictors.csv | 162 site-level rows of distance to edge, harvest year, and understory height; required for script 01. |
| landbase_accounting_2003_2025.csv | Aggregate areas and rounded annual rates used in the three-species projections; required for script 06. |
| annual_old_forest_rates_2003_2025.csv | Unrounded GIS-derived annual rates, retained as provenance. Script 06 uses the rounded accounting values to reproduce Table 3. |
| hawkears_thresholds_8recs.csv | Saved threshold-estimation output, validation counts, and detected site-year totals; provenance, not an input to the occupancy fit. |

These are processed analysis inputs. Proprietary Al-Pac vegetation inventory and harvest spatial layers, raw audio, detector-window tables, and development analyses are not distributed here. Aggregate landscape accounting values allow projections without those spatial layers.

### Response variables

| Column | Definition |
| --- | --- |
| location | Site identifier and join key. Identifiers retain historical CC/OG suffixes; use the treatment column for analysis labels. |
| transect | Spatial cluster identifier. |
| treatment | CLH, OF, or UP. |
| year | Calendar year of the recording. |
| doy | Day of year (146–197 in the revised data). |
| tod | Recording start time in minutes after midnight (286–460). |
| BBWA, BRCR, BTNW, HETH, TEWA, WIWR | Binary recording detections: 1 = detected, 0 = not detected. |

### Site predictors

| Column | Definition |
| --- | --- |
| location | Unique site identifier matching response.csv. |
| distance_to_edge_m | Distance to treatment edge, in metres. |
| harvest_year | Harvest year; NA for unharvested OF sites. |
| uheight | Retained understory height, in metres, used only at UP sites. Missing UP values are imputed with the observed UP mean before scaling; non-UP model inputs are zero. |

Script 01 saves scaling values and maps with the Stan input. To reproduce the fitted analysis, TSH is standardized using the 2015 values across all sites, including OF zeros; detection covariates use the species-long recording table. DTE is standardized across sites and height across UP sites. Structural trajectory predictions use treatment-specific median DTE, median UP height, zero cluster effect, and no calendar-year effects. These references differ from mean-covariate baseline summaries.

### Acoustic and landbase provenance

In the threshold table, threshold is the score at modelled P(true positive) = 0.90; threshold_lo90/hi90 are saved simulation limits. n_true/n_false count validated windows, roc_auc_sample is the unweighted validation-sample ROC-AUC, and n_above/false_above count validated windows above threshold. Precision and its bootstrap limits are weighted by confidence-score × recording-origin strata. detected_siteyears and n_siteyears give detected and total site-years. These are saved computational outputs; the table does not contain recall. The Brown Creeper empirical bootstrap precision interval is degenerate at 1 because no validated retained window was false; it does not imply certainty about population precision.

The landbase accounting table gives component, value_ha, and definition; annual components have units ha/year. The companion rates table gives annual inflow, harvest loss, fire loss, and signed net OF change, all in ha/year.

## Software

The fitted analysis used CmdStan 2.38.0 and cmdstanr 0.9.0. The revised workflow was checked with R 4.5.3. Required R packages are tidyverse, cmdstanr, posterior, and loo. The optional plots use patchwork; dotenv is needed only when loading a local .env file.

Install CmdStan and its C++ toolchain using cmdstanr. If cmdstanr already finds the installation, no environment file is needed. Otherwise copy .env.example to .env and set CMDSTAN_PATH to your installation directory. The local .env file is excluded from Git.

## Run the analysis

Run from the repository root, where Analyses/ and Data/ are located:

    source("Analyses/01_prepare_stan_data.R")
    source("Analyses/02_fit_dynocc_model.R")
    source("Analyses/03_extract_posterior_summaries.R")
    source("Analyses/06_project_habitat_supply.R")
    source("Analyses/08_summarise_naive_occurrence.R")

Script 02 runs four chains with 1,000 warm-up and 1,000 sampling iterations each, using the saved analysis seed (2133273333). Fitting requires substantial time and memory. Monte Carlo results can vary across software versions and platforms.

| Existing file in Analyses/ | Main output or purpose |
| --- | --- |
| 01_prepare_stan_data.R | Data/REVISED/full_data_for_stan.rds, recreated from revised CSVs; never overwrites the original cache in Data/. |
| 02_fit_dynocc_model.R | Self-contained Results/Model fits/Model_fit.rds, persistent chain CSVs, convergence diagnostics, E-BFMI, and session information. |
| dynocc_multispecies.stan | Final model, site-by-species log likelihoods, contrasts, and generated trajectories. |
| 03_extract_posterior_summaries.R | CSVs in Results/Summaries/: detection, covariate and treatment effects, structural trajectories, convergence, and species-level LOO diagnostics. |
| 06_project_habitat_supply.R | Three-species scenario trajectories, occupancy lookup summaries, and Table 3 endpoint changes in Results/Projections/. |
| 08_summarise_naive_occurrence.R | Observed occurrence summaries in Results/Summaries/. |

The existing 04_plot_dynocc_figures.R and 07_plot_habitat_supply.R are optional reference plotting scripts. They are not required to produce analysis tables. Script 07 requires script 06 first. The existing 05_extract_OF_brcr_trends.R is an optional structural-trajectory diagnostic: its year index is elapsed trajectory time and its BRCR change scenarios are **not** inputs to Table 3. Its filename is retained for compatibility. No additional scripts are needed.

## Projection assumptions

Projections preserve the Table 3 computation and evaluate BRCR, WIWR, and BBWA over 22 projection years. They start with 455,169 ha of OF (83,502 reserved and 371,667 harvestable), 22,336 ha of CLH, and 5,192 ha of UP.

All scenarios use net OF loss of 1,403 ha/year: fire loss of 4,803 plus harvest loss of 595 minus successional inflow of 3,995. New annual harvested area is 27,528 / 22 ha. Only new harvest allocation changes: S1 uses the current UP fraction (5,192 / 27,528), S2 uses 0% UP, and S3 uses 20% UP.

Initial harvested area is distributed equally across 22 age classes, indexed 0–21. Each year cohorts age, the oldest exits the accounting window, and a new age-zero cohort enters. Total harvested area therefore stays constant. This is an aggregate accounting assumption, not the observed distribution of mapped harvest cohorts. The OF pool changes separately; burned area is never assigned to a harvest treatment.

For each starting age, the occupancy lookup averages an 11-year dynamic trajectory in which stand age advances annually. Starting ages 0–21 therefore evaluate ages through 31, with year and cluster effects set to zero. These values are trajectory averages, not occupancy at a single age. The 22-year accounting window is not an empirical UP age limit: observed UP survey ages extend to 20 years. Predictions beyond observed ages involve extrapolation.

Expected occupied area sums each pool's area multiplied by its occupancy probability. Occupied percentage divides that value by remaining focal area in each projection year. Table 3 changes are differences between endpoint posterior medians, rather than posterior medians of draw-wise differences. No future temporal decline in OF occupancy is imposed in addition to the specified area loss.

## License

Code is released under the MIT License; see [LICENSE](LICENSE). Consult Borealis for dataset licensing and version information.
