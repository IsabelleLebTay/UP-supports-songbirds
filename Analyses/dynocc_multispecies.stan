// Final multi-species dynamic occupancy model.
//
// Includes treatment-specific time-since-harvest effects, distance-to-edge
// effects, cluster random effects, treatment-specific temporal random walks,
// and understory height as a UP-only covariate on initial occupancy.

data {
  int<lower=1> n_species;
  int<lower=1> n_sites;
  int<lower=1> n_visits;
  int<lower=1> n_obs;
  int<lower=1> n_treatments;
  int<lower=1> n_years;
  int<lower=1> n_clusters;

  array[n_sites] int<lower=1, upper=n_treatments> treatment;
  array[n_sites] int<lower=1, upper=n_clusters> cluster;
  array[n_sites] int<lower=1> n_visits_site;
  array[n_sites] int<lower=1> first_visit_idx;
  array[n_sites, n_years] real tsh;
  array[n_sites] real dte;

  // Understory height (standardized), meaningful for UP sites only.
  // CLH and OF sites should be set to 0 in the data so the uheight
  // term contributes nothing for those treatments.
  array[n_sites] real uheight;

  array[n_visits] int<lower=1, upper=n_sites> visit_site;
  array[n_visits] int<lower=1, upper=n_years> visit_year;
  array[n_visits] int<lower=1> n_recs_visit;
  array[n_visits] int<lower=1> first_rec_idx;

  array[n_species, n_obs] int<lower=0, upper=1> y;
  array[n_obs] real doy;
  array[n_obs] real tod;

  int<lower=1> n_tsh_grid;
  array[n_tsh_grid, n_years] real tsh_grid;
  array[n_treatments] real dte_ref;

  // Reference uheight value for UP trajectory projections (e.g. mean)
  real uheight_ref;
}

parameters {
  real mu_alpha_psi;
  real mu_alpha_phi;
  real mu_alpha_gam;
  real mu_alpha_p;

  real<lower=0> sigma_alpha_psi;
  real<lower=0> sigma_alpha_phi;
  real<lower=0> sigma_alpha_gam;
  real<lower=0> sigma_alpha_p;

  real mu_beta_doy_p;
  real mu_beta_tod_p;
  real<lower=0> sigma_beta_doy_p;
  real<lower=0> sigma_beta_tod_p;

  real mu_beta_psi_CLH;
  real mu_beta_psi_UP;
  real<lower=0> sigma_beta_psi;

  array[n_species] real alpha_psi_raw;
  array[n_species] real alpha_phi_raw;
  array[n_species] real alpha_gam_raw;
  array[n_species] real alpha_p_raw;

  array[n_species] real beta_doy_p_raw;
  array[n_species] real beta_tod_p_raw;

  array[n_species] real beta_psi_CLH_raw;
  array[n_species] real beta_psi_UP_raw;

  // Understory height effect on psi (UP only)
  real mu_beta_uheight_psi;
  real<lower=0> sigma_beta_uheight_psi;
  array[n_species] real beta_uheight_psi_raw;

  // ---- Treatment-specific TSH effects ----
  real mu_beta_tsh_psi_CLH;
  real<lower=0> sigma_beta_tsh_psi_CLH;
  array[n_species] real beta_tsh_psi_CLH_raw;

  real mu_beta_tsh_phi_CLH;
  real<lower=0> sigma_beta_tsh_phi_CLH;
  array[n_species] real beta_tsh_phi_CLH_raw;

  real mu_beta_tsh_gam_CLH;
  real<lower=0> sigma_beta_tsh_gam_CLH;
  array[n_species] real beta_tsh_gam_CLH_raw;

  real mu_beta_tsh_psi_UP;
  real<lower=0> sigma_beta_tsh_psi_UP;
  array[n_species] real beta_tsh_psi_UP_raw;

  real mu_beta_tsh_phi_UP;
  real<lower=0> sigma_beta_tsh_phi_UP;
  array[n_species] real beta_tsh_phi_UP_raw;

  real mu_beta_tsh_gam_UP;
  real<lower=0> sigma_beta_tsh_gam_UP;
  array[n_species] real beta_tsh_gam_UP_raw;

  // Distance-to-edge: treatment-specific
  real mu_community_dte;
  real<lower=0> sigma_community_dte;
  array[n_species] real mu_beta_dte_raw;
  real<lower=0> sigma_dte_trt_dev;
  array[n_species] real beta_dte_CLH_raw;
  array[n_species] real beta_dte_UP_raw;

  // Cluster random effects
  real mu_log_sigma_cluster;
  real<lower=0> sigma_log_sigma_cluster;
  array[n_species] real sigma_cluster_raw;
  array[n_species, n_clusters] real cluster_effect_raw;

  // Treatment-specific year random walks
  array[n_treatments] real<lower=0> sigma_year_phi;
  array[n_treatments] real<lower=0> sigma_year_gam;
  array[n_treatments] real year_init_phi;
  array[n_treatments] real year_init_gam;
  array[n_treatments, n_years - 1] real year_raw_phi;
  array[n_treatments, n_years - 1] real year_raw_gam;

  // Species deviations from treatment-specific year effects
  real<lower=0> sigma_species_year_phi;
  real<lower=0> sigma_species_year_gam;
  array[n_species, n_treatments, n_years] real species_year_phi_raw;
  array[n_species, n_treatments, n_years] real species_year_gam_raw;
}

transformed parameters {
  array[n_species] real alpha_psi;
  array[n_species] real alpha_phi;
  array[n_species] real alpha_gam;
  array[n_species] real alpha_p;
  array[n_species] real beta_doy_p;
  array[n_species] real beta_tod_p;

  for (i in 1:n_species) {
    alpha_psi[i] = mu_alpha_psi + sigma_alpha_psi * alpha_psi_raw[i];
    alpha_phi[i] = mu_alpha_phi + sigma_alpha_phi * alpha_phi_raw[i];
    alpha_gam[i] = mu_alpha_gam + sigma_alpha_gam * alpha_gam_raw[i];
    alpha_p[i] = mu_alpha_p + sigma_alpha_p * alpha_p_raw[i];
    beta_doy_p[i] = mu_beta_doy_p + sigma_beta_doy_p * beta_doy_p_raw[i];
    beta_tod_p[i] = mu_beta_tod_p + sigma_beta_tod_p * beta_tod_p_raw[i];
  }

  array[n_species, n_treatments] real beta_trt_psi;

  // Understory height effect on psi (species-specific, UP only)
  array[n_species] real beta_uheight_psi;

  for (i in 1:n_species) {
    // CLH (treatment 1)
    beta_trt_psi[i, 1] = mu_beta_psi_CLH + sigma_beta_psi * beta_psi_CLH_raw[i];

    // OF (treatment 2) - reference
    beta_trt_psi[i, 2] = 0;

    // UP (treatment 3)
    beta_trt_psi[i, 3] = mu_beta_psi_UP + sigma_beta_psi * beta_psi_UP_raw[i];

    // uheight: community hyperprior -> species
    beta_uheight_psi[i] = mu_beta_uheight_psi + sigma_beta_uheight_psi * beta_uheight_psi_raw[i];
  }

  // ---- Treatment-specific TSH slopes ----
  array[n_species, n_treatments] real beta_tsh_psi;
  array[n_species, n_treatments] real beta_tsh_phi;
  array[n_species, n_treatments] real beta_tsh_gam;

  for (i in 1:n_species) {
    // CLH (treatment 1)
    beta_tsh_psi[i, 1] = mu_beta_tsh_psi_CLH + sigma_beta_tsh_psi_CLH * beta_tsh_psi_CLH_raw[i];
    beta_tsh_phi[i, 1] = mu_beta_tsh_phi_CLH + sigma_beta_tsh_phi_CLH * beta_tsh_phi_CLH_raw[i];
    beta_tsh_gam[i, 1] = mu_beta_tsh_gam_CLH + sigma_beta_tsh_gam_CLH * beta_tsh_gam_CLH_raw[i];

    // OF (treatment 2) - no TSH effect
    beta_tsh_psi[i, 2] = 0;
    beta_tsh_phi[i, 2] = 0;
    beta_tsh_gam[i, 2] = 0;

    // UP (treatment 3)
    beta_tsh_psi[i, 3] = mu_beta_tsh_psi_UP + sigma_beta_tsh_psi_UP * beta_tsh_psi_UP_raw[i];
    beta_tsh_phi[i, 3] = mu_beta_tsh_phi_UP + sigma_beta_tsh_phi_UP * beta_tsh_phi_UP_raw[i];
    beta_tsh_gam[i, 3] = mu_beta_tsh_gam_UP + sigma_beta_tsh_gam_UP * beta_tsh_gam_UP_raw[i];
  }

  // DTE: species-level mean + treatment deviations (OF = reference)
  array[n_species] real mu_beta_dte;
  array[n_species, n_treatments] real beta_dte;
  for (i in 1:n_species) {
    mu_beta_dte[i] = mu_community_dte + sigma_community_dte * mu_beta_dte_raw[i];
    beta_dte[i, 1] = mu_beta_dte[i] + sigma_dte_trt_dev * beta_dte_CLH_raw[i];
    beta_dte[i, 2] = mu_beta_dte[i];
    beta_dte[i, 3] = mu_beta_dte[i] + sigma_dte_trt_dev * beta_dte_UP_raw[i];
  }

  array[n_species] real<lower=0> sigma_cluster;
  array[n_species, n_clusters] real cluster_effect;
  for (i in 1:n_species) {
    sigma_cluster[i] = exp(mu_log_sigma_cluster + sigma_log_sigma_cluster * sigma_cluster_raw[i]);
    for (c in 1:n_clusters) {
      cluster_effect[i, c] = sigma_cluster[i] * cluster_effect_raw[i, c];
    }
  }

  array[n_treatments, n_years] real year_effect_phi;
  array[n_treatments, n_years] real year_effect_gam;
  for (k in 1:n_treatments) {
    year_effect_phi[k, 1] = year_init_phi[k];
    year_effect_gam[k, 1] = year_init_gam[k];
    for (t in 2:n_years) {
      year_effect_phi[k, t] = year_effect_phi[k, t-1] + sigma_year_phi[k] * year_raw_phi[k, t-1];
      year_effect_gam[k, t] = year_effect_gam[k, t-1] + sigma_year_gam[k] * year_raw_gam[k, t-1];
    }
  }

  array[n_species, n_treatments, n_years] real species_year_effect_phi;
  array[n_species, n_treatments, n_years] real species_year_effect_gam;
  for (i in 1:n_species) {
    for (k in 1:n_treatments) {
      for (t in 1:n_years) {
        species_year_effect_phi[i, k, t] = year_effect_phi[k, t] +
          sigma_species_year_phi * species_year_phi_raw[i, k, t];
        species_year_effect_gam[i, k, t] = year_effect_gam[k, t] +
          sigma_species_year_gam * species_year_gam_raw[i, k, t];
      }
    }
  }
}

model {
  // === Priors: Community hyperparameters ===
  mu_alpha_psi ~ normal(0, 1.5);
  mu_alpha_phi ~ normal(1.5, 1);
  mu_alpha_gam ~ normal(-1.5, 1);
  mu_alpha_p ~ normal(0, 1.5);

  sigma_alpha_psi ~ normal(0, 1);
  sigma_alpha_phi ~ normal(0, 1);
  sigma_alpha_gam ~ normal(0, 1);
  sigma_alpha_p ~ normal(0, 1);

  mu_beta_doy_p ~ normal(0, 1);
  mu_beta_tod_p ~ normal(0, 1);
  sigma_beta_doy_p ~ normal(0, 0.5);
  sigma_beta_tod_p ~ normal(0, 0.5);

  mu_beta_psi_CLH ~ normal(0, 1);
  mu_beta_psi_UP ~ normal(0, 1);
  sigma_beta_psi ~ normal(0, 0.5);

  // Species-specific raw parameters
  alpha_psi_raw ~ std_normal();
  alpha_phi_raw ~ std_normal();
  alpha_gam_raw ~ std_normal();
  alpha_p_raw ~ std_normal();
  beta_doy_p_raw ~ std_normal();
  beta_tod_p_raw ~ std_normal();

  // Treatment-specific TSH effects - CLH
  mu_beta_tsh_psi_CLH ~ normal(0, 1);
  sigma_beta_tsh_psi_CLH ~ normal(0, 0.5);
  beta_tsh_psi_CLH_raw ~ std_normal();

  mu_beta_tsh_phi_CLH ~ normal(0, 1);
  sigma_beta_tsh_phi_CLH ~ normal(0, 0.5);
  beta_tsh_phi_CLH_raw ~ std_normal();

  mu_beta_tsh_gam_CLH ~ normal(0, 1);
  sigma_beta_tsh_gam_CLH ~ normal(0, 0.5);
  beta_tsh_gam_CLH_raw ~ std_normal();

  // Treatment-specific TSH effects - UP
  mu_beta_tsh_psi_UP ~ normal(0, 1);
  sigma_beta_tsh_psi_UP ~ normal(0, 0.5);
  beta_tsh_psi_UP_raw ~ std_normal();

  mu_beta_tsh_phi_UP ~ normal(0, 1);
  sigma_beta_tsh_phi_UP ~ normal(0, 0.5);
  beta_tsh_phi_UP_raw ~ std_normal();

  mu_beta_tsh_gam_UP ~ normal(0, 1);
  sigma_beta_tsh_gam_UP ~ normal(0, 0.5);
  beta_tsh_gam_UP_raw ~ std_normal();

  // Distance-to-edge (treatment-specific)
  mu_community_dte ~ normal(0, 1);
  sigma_community_dte ~ normal(0, 0.5);
  mu_beta_dte_raw ~ std_normal();
  sigma_dte_trt_dev ~ normal(0, 0.5);
  beta_dte_CLH_raw ~ std_normal();
  beta_dte_UP_raw ~ std_normal();

  // Treatment effects on psi
  beta_psi_CLH_raw ~ std_normal();
  beta_psi_UP_raw ~ std_normal();

  // Understory height on psi (UP only)
  mu_beta_uheight_psi ~ normal(0, 1);
  sigma_beta_uheight_psi ~ normal(0, 0.5);
  beta_uheight_psi_raw ~ std_normal();

  // Cluster random effects
  mu_log_sigma_cluster ~ normal(-1, 0.5);
  sigma_log_sigma_cluster ~ normal(0, 0.25);
  sigma_cluster_raw ~ std_normal();
  for (i in 1:n_species) {
    cluster_effect_raw[i] ~ std_normal();
  }

  // Year random walks
  sigma_year_phi ~ normal(0, 0.5);
  sigma_year_gam ~ normal(0, 0.5);
  year_init_phi ~ normal(0, 1);
  year_init_gam ~ normal(0, 1);
  for (k in 1:n_treatments) {
    year_raw_phi[k] ~ std_normal();
    year_raw_gam[k] ~ std_normal();
  }

  // Species deviations from year effects
  sigma_species_year_phi ~ normal(0, 0.3);
  sigma_species_year_gam ~ normal(0, 0.3);
  for (i in 1:n_species) {
    for (k in 1:n_treatments) {
      species_year_phi_raw[i, k] ~ std_normal();
      species_year_gam_raw[i, k] ~ std_normal();
    }
  }

  // === Likelihood ===
  for (sp in 1:n_species) {
    for (s in 1:n_sites) {
      int trt = treatment[s];
      int clust = cluster[s];
      int n_vis = n_visits_site[s];
      int first_vis = first_visit_idx[s];

      array[2] real log_forward;

      // First visit
      {
        int v = first_vis;
        int yr_first = visit_year[v];
        int first_rec = first_rec_idx[v];
        int n_recs = n_recs_visit[v];

        // uheight enters psi only; data zeroed for CLH/OF so term vanishes
        real psi_s = inv_logit(alpha_psi[sp] + beta_trt_psi[sp, trt] +
                                beta_tsh_psi[sp, trt] * tsh[s, yr_first] +
                                beta_dte[sp, trt] * dte[s] +
                                beta_uheight_psi[sp] * uheight[s] +
                                cluster_effect[sp, clust]);

        real log_p_data_z1 = 0;
        int any_det = 0;
        for (k in 1:n_recs) {
          int rec_idx = first_rec + k - 1;
          real p_k = inv_logit(alpha_p[sp] +
                                beta_doy_p[sp] * doy[rec_idx] +
                                beta_tod_p[sp] * tod[rec_idx]);
          if (y[sp, rec_idx] == 1) {
            log_p_data_z1 += log(p_k);
            any_det = 1;
          } else {
            log_p_data_z1 += log1m(p_k);
          }
        }

        real log_p_data_z0 = (any_det == 1) ? negative_infinity() : 0;
        log_forward[1] = log1m(psi_s) + log_p_data_z0;
        log_forward[2] = log(psi_s) + log_p_data_z1;
      }

      // Subsequent visits
      for (vi in 2:n_vis) {
        int v = first_vis + vi - 1;
        int v_prev = first_vis + vi - 2;
        int yr = visit_year[v];
        int yr_prev = visit_year[v_prev];
        int first_rec = first_rec_idx[v];
        int n_recs = n_recs_visit[v];

        array[2, 2] real log_trans_compound;
        log_trans_compound[1, 1] = 0;
        log_trans_compound[1, 2] = negative_infinity();
        log_trans_compound[2, 1] = negative_infinity();
        log_trans_compound[2, 2] = 0;

        for (t in (yr_prev + 1):yr) {
          real phi_t = inv_logit(alpha_phi[sp] +
                                  beta_tsh_phi[sp, trt] * tsh[s, t] +
                                  beta_dte[sp, trt] * dte[s] +
                                  species_year_effect_phi[sp, trt, t] +
                                  cluster_effect[sp, clust]);
          real gam_t = inv_logit(alpha_gam[sp] +
                                  beta_tsh_gam[sp, trt] * tsh[s, t] +
                                  beta_dte[sp, trt] * dte[s] +
                                  species_year_effect_gam[sp, trt, t] +
                                  cluster_effect[sp, clust]);

          array[2, 2] real log_trans_t;
          log_trans_t[1, 1] = log1m(gam_t);
          log_trans_t[1, 2] = log(gam_t);
          log_trans_t[2, 1] = log1m(phi_t);
          log_trans_t[2, 2] = log(phi_t);

          array[2, 2] real log_trans_new;
          for (i in 1:2) {
            for (j in 1:2) {
              log_trans_new[i, j] = log_sum_exp(
                log_trans_compound[i, 1] + log_trans_t[1, j],
                log_trans_compound[i, 2] + log_trans_t[2, j]
              );
            }
          }
          log_trans_compound = log_trans_new;
        }

        real log_p_data_z1 = 0;
        int any_det = 0;
        for (k in 1:n_recs) {
          int rec_idx = first_rec + k - 1;
          real p_k = inv_logit(alpha_p[sp] +
                                beta_doy_p[sp] * doy[rec_idx] +
                                beta_tod_p[sp] * tod[rec_idx]);
          if (y[sp, rec_idx] == 1) {
            log_p_data_z1 += log(p_k);
            any_det = 1;
          } else {
            log_p_data_z1 += log1m(p_k);
          }
        }

        real log_p_data_z0 = (any_det == 1) ? negative_infinity() : 0;

        array[2] real log_forward_new;
        log_forward_new[1] = log_sum_exp(
          log_forward[1] + log_trans_compound[1, 1],
          log_forward[2] + log_trans_compound[2, 1]
        ) + log_p_data_z0;
        log_forward_new[2] = log_sum_exp(
          log_forward[1] + log_trans_compound[1, 2],
          log_forward[2] + log_trans_compound[2, 2]
        ) + log_p_data_z1;

        log_forward = log_forward_new;
      }

      target += log_sum_exp(log_forward[1], log_forward[2]);
    }
  }
}

generated quantities {
  // === Detection probabilities ===
  array[n_species] real p;
  for (i in 1:n_species) {
    p[i] = inv_logit(alpha_p[i]);
  }

  // === phi/gam by treatment x year (at zero covariates, zero cluster) ===
  array[n_species, n_treatments, n_years] real phi_by_trt_year;
  array[n_species, n_treatments, n_years] real gam_by_trt_year;
  array[n_species, n_treatments, n_years] real eps_by_trt_year;

  for (i in 1:n_species) {
    for (k in 1:n_treatments) {
      for (t in 1:n_years) {
        phi_by_trt_year[i, k, t] = inv_logit(alpha_phi[i] +
                                              species_year_effect_phi[i, k, t]);
        gam_by_trt_year[i, k, t] = inv_logit(alpha_gam[i] +
                                              species_year_effect_gam[i, k, t]);
        eps_by_trt_year[i, k, t] = 1 - phi_by_trt_year[i, k, t];
      }
    }
  }

  // === Baseline psi by treatment (at mean TSH, zero uheight) ===
  array[n_species, n_treatments] real psi_baseline;
  for (i in 1:n_species) {
    for (k in 1:n_treatments) {
      psi_baseline[i, k] = inv_logit(alpha_psi[i] + beta_trt_psi[i, k]);
    }
  }

  // =========================================================================
  // OCCUPANCY TRAJECTORY PROJECTIONS BY TREATMENT x TSH
  // Evaluated at reference site: treatment-specific DTE (dte_ref[k]),
  // zero cluster. For UP, uheight evaluated at uheight_ref.
  // For CLH and OF, uheight term is zero (no retained understory).
  // =========================================================================
  array[n_species, n_treatments, n_tsh_grid, n_years] real occ_trajectory;

  for (i in 1:n_species) {
    for (k in 1:n_treatments) {
      for (g in 1:n_tsh_grid) {
        // uheight enters initial occupancy for UP (treatment 3) only
        real uheight_term = (k == 3) ? beta_uheight_psi[i] * uheight_ref : 0.0;

        occ_trajectory[i, k, g, 1] = inv_logit(alpha_psi[i] + beta_trt_psi[i, k] +
                                                 beta_tsh_psi[i, k] * tsh_grid[g, 1] +
                                                 beta_dte[i, k] * dte_ref[k] +
                                                 uheight_term);
        for (t in 2:n_years) {
          real phi_t = inv_logit(alpha_phi[i] +
                                  beta_tsh_phi[i, k] * tsh_grid[g, t] +
                                  beta_dte[i, k] * dte_ref[k] +
                                  species_year_effect_phi[i, k, t]);
          real gam_t = inv_logit(alpha_gam[i] +
                                  beta_tsh_gam[i, k] * tsh_grid[g, t] +
                                  beta_dte[i, k] * dte_ref[k] +
                                  species_year_effect_gam[i, k, t]);
          occ_trajectory[i, k, g, t] = occ_trajectory[i, k, g, t-1] * phi_t +
                                        (1 - occ_trajectory[i, k, g, t-1]) * gam_t;
        }
      }
    }
  }

  // === Community-level occupancy trajectories ===
  array[n_treatments, n_tsh_grid, n_years] real occ_trajectory_community;
  for (k in 1:n_treatments) {
    for (g in 1:n_tsh_grid) {
      for (t in 1:n_years) {
        occ_trajectory_community[k, g, t] = mean(occ_trajectory[, k, g, t]);
      }
    }
  }

  // === Treatment contrasts on occupancy trajectories ===
  array[n_species, n_tsh_grid, n_years] real occ_UP_vs_OF;
  array[n_species, n_tsh_grid, n_years] real occ_CLH_vs_OF;
  array[n_species, n_tsh_grid, n_years] real occ_UP_vs_CLH;

  for (i in 1:n_species) {
    for (g in 1:n_tsh_grid) {
      for (t in 1:n_years) {
        occ_UP_vs_OF[i, g, t] = occ_trajectory[i, 3, g, t] - occ_trajectory[i, 2, g, t];
        occ_CLH_vs_OF[i, g, t] = occ_trajectory[i, 1, g, t] - occ_trajectory[i, 2, g, t];
        occ_UP_vs_CLH[i, g, t] = occ_trajectory[i, 3, g, t] - occ_trajectory[i, 1, g, t];
      }
    }
  }

  // === Community-level trajectory contrasts ===
  array[n_tsh_grid, n_years] real occ_UP_vs_OF_community;
  array[n_tsh_grid, n_years] real occ_CLH_vs_OF_community;
  array[n_tsh_grid, n_years] real occ_UP_vs_CLH_community;

  for (g in 1:n_tsh_grid) {
    for (t in 1:n_years) {
      occ_UP_vs_OF_community[g, t] = mean(occ_UP_vs_OF[, g, t]);
      occ_CLH_vs_OF_community[g, t] = mean(occ_CLH_vs_OF[, g, t]);
      occ_UP_vs_CLH_community[g, t] = mean(occ_UP_vs_CLH[, g, t]);
    }
  }

  // === Shared year effects ===
  array[n_treatments, n_years] real year_phi_shared;
  array[n_treatments, n_years] real year_gam_shared;
  for (k in 1:n_treatments) {
    for (t in 1:n_years) {
      year_phi_shared[k, t] = year_effect_phi[k, t];
      year_gam_shared[k, t] = year_effect_gam[k, t];
    }
  }

  // === Year-specific treatment contrasts on phi/gam ===
  array[n_species, n_years] real phi_UP_vs_CLH_by_year;
  array[n_species, n_years] real phi_UP_vs_OF_by_year;
  array[n_species, n_years] real phi_CLH_vs_OF_by_year;
  array[n_species, n_years] real gam_UP_vs_CLH_by_year;
  array[n_species, n_years] real gam_UP_vs_OF_by_year;
  array[n_species, n_years] real gam_CLH_vs_OF_by_year;

  for (i in 1:n_species) {
    for (t in 1:n_years) {
      phi_UP_vs_CLH_by_year[i, t] = phi_by_trt_year[i, 3, t] - phi_by_trt_year[i, 1, t];
      phi_UP_vs_OF_by_year[i, t] = phi_by_trt_year[i, 3, t] - phi_by_trt_year[i, 2, t];
      phi_CLH_vs_OF_by_year[i, t] = phi_by_trt_year[i, 1, t] - phi_by_trt_year[i, 2, t];
      gam_UP_vs_CLH_by_year[i, t] = gam_by_trt_year[i, 3, t] - gam_by_trt_year[i, 1, t];
      gam_UP_vs_OF_by_year[i, t] = gam_by_trt_year[i, 3, t] - gam_by_trt_year[i, 2, t];
      gam_CLH_vs_OF_by_year[i, t] = gam_by_trt_year[i, 1, t] - gam_by_trt_year[i, 2, t];
    }
  }

  // === Psi contrasts ===
  array[n_species] real psi_UP_vs_CLH;
  array[n_species] real psi_UP_vs_OF;
  array[n_species] real psi_CLH_vs_OF;
  for (i in 1:n_species) {
    psi_UP_vs_CLH[i] = psi_baseline[i, 3] - psi_baseline[i, 1];
    psi_UP_vs_OF[i] = psi_baseline[i, 3] - psi_baseline[i, 2];
    psi_CLH_vs_OF[i] = psi_baseline[i, 1] - psi_baseline[i, 2];
  }

  // === TSH effect outputs (treatment-specific) ===
  array[n_species, n_treatments] real beta_tsh_psi_out;
  array[n_species, n_treatments] real beta_tsh_phi_out;
  array[n_species, n_treatments] real beta_tsh_gam_out;
  array[n_species, n_treatments] real beta_dte_out;
  array[n_species] real sigma_cluster_out;
  for (i in 1:n_species) {
    for (k in 1:n_treatments) {
      beta_tsh_psi_out[i, k] = beta_tsh_psi[i, k];
      beta_tsh_phi_out[i, k] = beta_tsh_phi[i, k];
      beta_tsh_gam_out[i, k] = beta_tsh_gam[i, k];
      beta_dte_out[i, k] = beta_dte[i, k];
    }
    sigma_cluster_out[i] = sigma_cluster[i];
  }

  // === TSH treatment contrasts ===
  array[n_species] real beta_tsh_psi_UP_vs_CLH;
  array[n_species] real beta_tsh_phi_UP_vs_CLH;
  array[n_species] real beta_tsh_gam_UP_vs_CLH;
  for (i in 1:n_species) {
    beta_tsh_psi_UP_vs_CLH[i] = beta_tsh_psi[i, 3] - beta_tsh_psi[i, 1];
    beta_tsh_phi_UP_vs_CLH[i] = beta_tsh_phi[i, 3] - beta_tsh_phi[i, 1];
    beta_tsh_gam_UP_vs_CLH[i] = beta_tsh_gam[i, 3] - beta_tsh_gam[i, 1];
  }

  // === Community-level year-specific summaries ===
  array[n_treatments, n_years] real phi_community_by_year;
  array[n_treatments, n_years] real gam_community_by_year;
  array[n_treatments, n_years] real eps_community_by_year;
  for (k in 1:n_treatments) {
    for (t in 1:n_years) {
      phi_community_by_year[k, t] = mean(phi_by_trt_year[, k, t]);
      gam_community_by_year[k, t] = mean(gam_by_trt_year[, k, t]);
      eps_community_by_year[k, t] = mean(eps_by_trt_year[, k, t]);
    }
  }

  // === Understory height effect outputs ===
  array[n_species] real beta_uheight_psi_out;
  for (i in 1:n_species) {
    beta_uheight_psi_out[i] = beta_uheight_psi[i];
  }

  // === DTE treatment contrasts ===
  array[n_species] real beta_dte_UP_vs_CLH;
  array[n_species] real beta_dte_UP_vs_OF;
  array[n_species] real beta_dte_CLH_vs_OF;
  for (i in 1:n_species) {
    beta_dte_UP_vs_CLH[i] = beta_dte[i, 3] - beta_dte[i, 1];
    beta_dte_UP_vs_OF[i] = beta_dte[i, 3] - beta_dte[i, 2];
    beta_dte_CLH_vs_OF[i] = beta_dte[i, 1] - beta_dte[i, 2];
  }

  // === Pointwise log-likelihood for LOO/WAIC ===
  array[n_species, n_sites] real log_lik;

  for (sp in 1:n_species) {
    for (s in 1:n_sites) {
      int trt = treatment[s];
      int clust = cluster[s];
      int n_vis = n_visits_site[s];
      int first_vis = first_visit_idx[s];

      array[2] real log_fwd;

      {
        int v = first_vis;
        int yr_first = visit_year[v];
        int first_rec = first_rec_idx[v];
        int n_recs = n_recs_visit[v];

        real psi_s = inv_logit(alpha_psi[sp] + beta_trt_psi[sp, trt] +
                                beta_tsh_psi[sp, trt] * tsh[s, yr_first] +
                                beta_dte[sp, trt] * dte[s] +
                                beta_uheight_psi[sp] * uheight[s] +
                                cluster_effect[sp, clust]);

        real log_p_data_z1 = 0;
        int any_det = 0;
        for (k in 1:n_recs) {
          int rec_idx = first_rec + k - 1;
          real p_k = inv_logit(alpha_p[sp] +
                                beta_doy_p[sp] * doy[rec_idx] +
                                beta_tod_p[sp] * tod[rec_idx]);
          if (y[sp, rec_idx] == 1) {
            log_p_data_z1 += log(p_k);
            any_det = 1;
          } else {
            log_p_data_z1 += log1m(p_k);
          }
        }

        real log_p_data_z0 = (any_det == 1) ? negative_infinity() : 0;
        log_fwd[1] = log1m(psi_s) + log_p_data_z0;
        log_fwd[2] = log(psi_s) + log_p_data_z1;
      }

      for (vi in 2:n_vis) {
        int v = first_vis + vi - 1;
        int v_prev = first_vis + vi - 2;
        int yr = visit_year[v];
        int yr_prev = visit_year[v_prev];
        int first_rec = first_rec_idx[v];
        int n_recs = n_recs_visit[v];

        array[2, 2] real log_trans_compound;
        log_trans_compound[1, 1] = 0;
        log_trans_compound[1, 2] = negative_infinity();
        log_trans_compound[2, 1] = negative_infinity();
        log_trans_compound[2, 2] = 0;

        for (t in (yr_prev + 1):yr) {
          real phi_t = inv_logit(alpha_phi[sp] +
                                  beta_tsh_phi[sp, trt] * tsh[s, t] +
                                  beta_dte[sp, trt] * dte[s] +
                                  species_year_effect_phi[sp, trt, t] +
                                  cluster_effect[sp, clust]);
          real gam_t = inv_logit(alpha_gam[sp] +
                                  beta_tsh_gam[sp, trt] * tsh[s, t] +
                                  beta_dte[sp, trt] * dte[s] +
                                  species_year_effect_gam[sp, trt, t] +
                                  cluster_effect[sp, clust]);

          array[2, 2] real log_trans_t;
          log_trans_t[1, 1] = log1m(gam_t);
          log_trans_t[1, 2] = log(gam_t);
          log_trans_t[2, 1] = log1m(phi_t);
          log_trans_t[2, 2] = log(phi_t);

          array[2, 2] real log_trans_new;
          for (i in 1:2) {
            for (j in 1:2) {
              log_trans_new[i, j] = log_sum_exp(
                log_trans_compound[i, 1] + log_trans_t[1, j],
                log_trans_compound[i, 2] + log_trans_t[2, j]
              );
            }
          }
          log_trans_compound = log_trans_new;
        }

        real log_p_data_z1 = 0;
        int any_det = 0;
        for (k in 1:n_recs) {
          int rec_idx = first_rec + k - 1;
          real p_k = inv_logit(alpha_p[sp] +
                                beta_doy_p[sp] * doy[rec_idx] +
                                beta_tod_p[sp] * tod[rec_idx]);
          if (y[sp, rec_idx] == 1) {
            log_p_data_z1 += log(p_k);
            any_det = 1;
          } else {
            log_p_data_z1 += log1m(p_k);
          }
        }

        real log_p_data_z0 = (any_det == 1) ? negative_infinity() : 0;

        array[2] real log_fwd_new;
        log_fwd_new[1] = log_sum_exp(
          log_fwd[1] + log_trans_compound[1, 1],
          log_fwd[2] + log_trans_compound[2, 1]
        ) + log_p_data_z0;
        log_fwd_new[2] = log_sum_exp(
          log_fwd[1] + log_trans_compound[1, 2],
          log_fwd[2] + log_trans_compound[2, 2]
        ) + log_p_data_z1;

        log_fwd = log_fwd_new;
      }

      log_lik[sp, s] = log_sum_exp(log_fwd[1], log_fwd[2]);
    }
  }
}
