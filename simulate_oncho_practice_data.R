## ==========================================================
## Onchocerciasis geostatistical simulations for inlabru workshop
## Topics:
## 1. Spatial + spatio-temporal, different outcome types
## 2. Joint modelling of multiple processes
## 3. Non-stationary spatial process
## 4. Hybrid ML + geostatistical model
## Country: Nigeria
## Note: Onchocerciasis prevalence is higher in the north
##       and strongly related to multiple predictors
##       (e.g. rivers, environment, population, etc.)
## ==========================================================

## ---------- Packages ----------
library(mvtnorm)   # for multivariate Gaussian simulation
library(dplyr)
library(tidyr)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(INLA)
library(inlabru)

set.seed(1234)

## ==========================================================
## Nigeria polygons (admin-1 and dissolved country)
## ==========================================================

# Nigeria admin-1
nga_admin1 <- ne_states(country = "Nigeria", returnclass = "sf")

# Country polygon (dissolved)
nga_poly <- st_union(nga_admin1)

# Work in UTM 32N (same as the rest of your code)
nga_poly_utm   <- st_transform(nga_poly, 32632)
nga_admin1_utm <- st_transform(nga_admin1, 32632)

## ==========================================================
## 0. Common helpers
## ==========================================================

## Approximate Nigeria bounding box (lon, lat)
sample_nigeria_coords <- function(n) {
  lon <- runif(n, min = 3, max = 14)   # ~ Nigeria longitudes
  lat <- runif(n, min = 4, max = 14)   # ~ Nigeria latitudes
  cbind(lon = lon, lat = lat)
}

sample_nigeria_utm_only <- function(n, polygon_utm = nga_poly_utm) {
  pts <- st_sample(polygon_utm, size = n, type = "random") |>
    st_as_sf()

  coords <- st_coordinates(pts)
  data.frame(
    utm_x = coords[, 1],
    utm_y = coords[, 2]
  )
}

## Exponential covariance function (isotropic)
exp_covariance <- function(coords, range = 700000, sigma2 = 0.5) {
  D <- as.matrix(dist(coords))
  sigma2 * exp(-D / range)
}

## Simulate a spatial Gaussian process on given coordinates
simulate_spatial_field <- function(coords, range = 300000, sigma2 = 1) {
  Sigma <- exp_covariance(coords, range = range, sigma2 = sigma2)
  mvtnorm::rmvnorm(1, sigma = Sigma)[1, ]
}

## Simple AR(1) temporal random effect
simulate_AR1 <- function(T, rho = 0.8, sigma2 = 0.5) {
  eps <- rnorm(T, mean = 0, sd = sqrt(sigma2))
  eta <- numeric(T)
  eta[1] <- eps[1] / sqrt(1 - rho^2)
  for (t in 2:T) {
    eta[t] <- rho * eta[t - 1] + eps[t]
  }
  eta
}

## Synthetic environmental and socio-demographic covariates
## designed for ONCHOCERCIASIS in Nigeria
## - prevalence higher in the north (north_index)
## - higher risk near "river bands"
## - multiple predictors available for modelling
simulate_covariates <- function(coords) {
  # coords is a data.frame with utm_x, utm_y
  x <- scale(coords[, "utm_x"])
  y <- scale(coords[, "utm_y"])
  n <- nrow(coords)

  # North index: larger values in the north
  north_index <- as.numeric(scale(y))  # roughly south -> low, north -> high

  # Elevation: gentle N–S gradient + noise
  elevation <- 300 + 60 * sin(pi * y) + rnorm(n, sd = 20)

  # NDVI: vegetation index, somewhat greener in the south
  ndvi <- plogis(0.5 - 0.8 * north_index + 0.3 * sin(pi * x))

  # Urban: more urbanisation towards the south and east
  urban_p <- plogis(-0.2 + 0.6 * x - 0.4 * north_index)
  urban   <- rbinom(n, 1, urban_p)

  # Rainfall: higher in the south, lower in the north
  rainfall <- 1500 - 350 * north_index + rnorm(n, sd = 80) # mm/year, rough

  # Temperature: higher in the north, lower in the south
  temperature <- 24 + 2.5 * north_index + rnorm(n, sd = 0.8) # °C, rough pattern

  # Population density: more people in some central/south urban belt
  pop_density <- exp(4 + 0.5 * x - 0.7 * north_index + rnorm(n, sd = 0.4))
  pop_density <- pmin(pop_density, quantile(pop_density, 0.99)) # trim extreme

  # River proximity: synthetic "river band" running roughly NE–SW
  river_axis <- as.numeric(scale(x + 0.2 * y))
  river_index <- exp(- (river_axis^2) / (2 * 0.6^2))  # high near river, ~0 far

  data.frame(
    elevation   = elevation,
    ndvi        = ndvi,
    urban       = urban,
    rainfall    = rainfall,
    temperature = temperature,
    pop_density = pop_density,
    river_index = river_index,
    north_index = north_index
  )
}

## ==========================================================
## 1. Spatial & spatio-temporal ONCHOCERCIASIS data
## ==========================================================

## ---------- 1A. Spatial Binomial prevalence ----------
## Think: community surveys of onchocerciasis in Nigeria
simulate_spatial_oncho_binomial <- function(n_sites = 250) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Latent spatial field (slowly varying baseline risk)
  S <- simulate_spatial_field(coords, range = 400000, sigma2 = 2.0)

  ## Linear predictor (logit-prevalence)
  ## Higher risk in the north (north_index), near rivers, in certain eco-zones
  beta0       <- -0.7   # baseline
  beta_elev   <- -0.0005
  beta_ndvi   <- 0.8
  beta_urban  <- -0.2
  beta_rain   <- 0.0002
  beta_temp   <- 0.12
  beta_river  <- 1.3
  beta_north  <- 0.7
  beta_popden <- 0.00015

  linpred <- beta0 +
    beta_elev   * covs$elevation +
    beta_ndvi   * covs$ndvi +
    beta_urban  * covs$urban +
    beta_rain   * covs$rainfall +
    beta_temp   * covs$temperature +
    beta_river  * covs$river_index +
    beta_north  * covs$north_index +
    beta_popden * covs$pop_density +
    S

  p <- plogis(linpred)

  N <- sample(50:200, n_sites, replace = TRUE)  # individuals examined
  Y <- rbinom(n_sites, size = N, prob = p)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_examined = N,
    n_pos = Y,
    prevalence_true = p,
    S = S
  ) %>%
    bind_cols(covs)
}

spatial_oncho_binom <- simulate_spatial_oncho_binomial()

## ---------- 1B. Spatial Poisson counts (oncho case counts / nodules) ----------
simulate_spatial_oncho_poisson <- function(n_sites = 250) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Latent spatial field
  S <- simulate_spatial_field(coords, range = 500000, sigma2 = 0.7)

  ## Population at risk (offset)
  pop <- round(runif(n_sites, 500, 5000))

  ## Log intensity of cases (e.g. nodules, positive skin snips)
  alpha0       <- -7.2
  alpha_elev   <- -0.0004
  alpha_ndvi   <- 0.6
  alpha_urban  <- -0.3
  alpha_rain   <- 0.0003
  alpha_temp   <- 0.15
  alpha_river  <- 1.5
  alpha_north  <- 0.8
  alpha_popden <- 0.00018

  log_rate <- alpha0 +
    alpha_elev   * covs$elevation +
    alpha_ndvi   * covs$ndvi +
    alpha_urban  * covs$urban +
    alpha_rain   * covs$rainfall +
    alpha_temp   * covs$temperature +
    alpha_river  * covs$river_index +
    alpha_north  * covs$north_index +
    alpha_popden * covs$pop_density +
    S

  lambda <- pop * exp(log_rate)
  cases  <- rpois(n_sites, lambda = lambda)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    population = pop,
    cases = cases,
    lambda_true = lambda,
    log_rate_true = log_rate,
    S = S
  ) %>%
    bind_cols(covs)
}

spatial_oncho_pois <- simulate_spatial_oncho_poisson()

## ---------- 1C. Spatio-temporal Binomial prevalence ----------
## e.g. repeated community surveys over several years
simulate_spatiotemporal_oncho_binomial <- function(n_sites = 200, T = 6) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Spatial field
  S <- simulate_spatial_field(coords, range = 400000, sigma2 = 0.8)

  ## Temporal AR(1) (e.g. 6 annual surveys)
  gamma_t <- simulate_AR1(T = T, rho = 0.7, sigma2 = 0.3)

  beta0       <- -0.8
  beta_elev   <- -0.0004
  beta_ndvi   <- 0.7
  beta_urban  <- -0.1
  beta_rain   <- 0.0002
  beta_temp   <- 0.1
  beta_river  <- 1.2
  beta_north  <- 0.6
  beta_popden <- 0.00012

  df <- expand.grid(
    site = 1:n_sites,
    time = 1:T
  ) %>%
    arrange(site, time)

  df$utm_x <- coords[df$site, "utm_x"]
  df$utm_y <- coords[df$site, "utm_y"]
  df$S   <- S[df$site]
  df$gamma_t <- gamma_t[df$time]

  df <- df %>%
    left_join(
      covs %>% mutate(site = 1:n_sites),
      by = "site"
    )

  linpred <- beta0 +
    beta_elev   * df$elevation +
    beta_ndvi   * df$ndvi +
    beta_urban  * df$urban +
    beta_rain   * df$rainfall +
    beta_temp   * df$temperature +
    beta_river  * df$river_index +
    beta_north  * df$north_index +
    beta_popden * df$pop_density +
    df$S +
    df$gamma_t

  p <- plogis(linpred)
  N <- sample(50:200, nrow(df), replace = TRUE)
  Y <- rbinom(nrow(df), size = N, prob = p)

  df$n_examined <- N
  df$n_pos <- Y
  df$prevalence_true <- p

  df
}

st_oncho_binom <- simulate_spatiotemporal_oncho_binomial()

## ==========================================================
## 2. Joint modelling of multiple onchocerciasis processes
##    e.g. two risk groups or two diagnostic methods
## ==========================================================

simulate_joint_oncho_processes <- function(n_sites = 250) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Shared spatial field
  S_shared <- simulate_spatial_field(coords, range = 500000, sigma2 = 0.7)

  ## Outcome-specific small-scale fields
  S_group1 <- simulate_spatial_field(coords, range = 200000, sigma2 = 0.25)
  S_group2 <- simulate_spatial_field(coords, range = 200000, sigma2 = 0.25)

  ## Group 1: e.g. adults
  beta0_g1 <- -0.9
  linpred_g1 <- beta0_g1 +
    0.8  * covs$ndvi -
    0.0004 * covs$elevation +
    0.9  * covs$river_index +
    0.6  * covs$north_index +
    0.00015 * covs$pop_density +
    S_shared + S_group1
  p_g1 <- plogis(linpred_g1)
  N_g1 <- sample(60:220, n_sites, replace = TRUE)
  Y_g1 <- rbinom(n_sites, N_g1, p_g1)

  ## Group 2: e.g. children or alternative diagnostic
  beta0_g2 <- -1.1
  linpred_g2 <- beta0_g2 +
    0.7  * covs$ndvi -
    0.0003 * covs$elevation +
    0.7  * covs$river_index +
    0.4  * covs$north_index +
    0.00012 * covs$pop_density +
    S_shared + S_group2
  p_g2 <- plogis(linpred_g2)
  N_g2 <- sample(40:160, n_sites, replace = TRUE)
  Y_g2 <- rbinom(n_sites, N_g2, p_g2)

  df_g1 <- data.frame(
    id = 1:n_sites,
    group = "adults",
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_examined = N_g1,
    n_pos = Y_g1,
    prevalence_true = p_g1,
    S_shared = S_shared,
    S_specific = S_group1
  ) %>%
    bind_cols(covs)

  df_g2 <- data.frame(
    id = 1:n_sites,
    group = "children_or_alt_diag",
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_examined = N_g2,
    n_pos = Y_g2,
    prevalence_true = p_g2,
    S_shared = S_shared,
    S_specific = S_group2
  ) %>%
    bind_cols(covs)

  bind_rows(df_g1, df_g2)
}

joint_oncho <- simulate_joint_oncho_processes()

## ==========================================================
## 3. Non-stationary spatial onchocerciasis process
##    e.g. smoother in arid north, rougher in southern belt
## ==========================================================

simulate_nonstationary_oncho <- function(n_sites = 300) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Two spatial fields with different smoothness / range
  S_smooth <- simulate_spatial_field(coords, range = 650000, sigma2 = 0.5)
  S_rough  <- simulate_spatial_field(coords, range = 180000, sigma2 = 0.5)

  ## Spatially varying mixing weight (based on latitude)
  ## More "smooth" in the north, rougher in the south
  utm_y <- coords[, "utm_y"]
  w <- (utm_y - min(utm_y)) / (max(utm_y) - min(utm_y))  # 0 in south, 1 in north

  S_ns <- sqrt(w) * S_smooth + sqrt(1 - w) * S_rough

  ## Onchocerciasis prevalence with non-stationary spatial effect
  beta0       <- -0.8
  beta_ndvi   <- 0.7
  beta_elev   <- -0.0003
  beta_river  <- 1.0
  beta_north  <- 0.8
  beta_popden <- 0.00012

  linpred <- beta0 +
    beta_ndvi   * covs$ndvi -
    beta_elev   * covs$elevation +
    beta_river  * covs$river_index +
    beta_north  * covs$north_index +
    beta_popden * covs$pop_density +
    S_ns

  p <- plogis(linpred)
  N <- sample(50:200, n_sites, replace = TRUE)
  Y <- rbinom(n_sites, N, p)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_examined = N,
    n_pos = Y,
    prevalence_true = p,
    S_nonstationary = S_ns,
    weight_north = w,
    S_smooth = S_smooth,
    S_rough = S_rough
  ) %>%
    bind_cols(covs)
}

nonstat_oncho <- simulate_nonstationary_oncho()

## ==========================================================
## 4. Hybrid ML + geostatistical onchocerciasis data
##    (complex covariate effects + residual spatial field)
## ==========================================================

simulate_hybrid_ml_geo_oncho <- function(n_sites = 400) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## "True" complex ML-type signal (non-linear, interactions)
  ## approximating what RF / XGBoost might learn
  ml_signal <- with(covs, {
    f1 <- -2 * (ndvi - 0.5)^2                # peak risk around NDVI ~ 0.5
    f2 <- -0.00004 * (elevation - 250)^2     # higher risk at low–moderate elevation
    f3 <-  1.2 * river_index                 # strong river effect
    f4 <-  0.8 * north_index                 # higher risk in the north
    f5 <-  0.00015 * pop_density             # mildly higher in denser areas
    f6 <-  0.3 * sin(coords[, "utm_x"] / 150000) # east–west oscillation
    f1 + f2 + f3 + f4 + f5 + f6
  })

  ## Residual spatial field (to be modelled by inlabru)
  S_res <- simulate_spatial_field(coords, range = 400000, sigma2 = 0.3)

  beta0 <- -1.0
  linpred <- beta0 + ml_signal + S_res
  p <- plogis(linpred)

  N <- sample(50:200, n_sites, replace = TRUE)
  Y <- rbinom(n_sites, N, p)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_examined = N,
    n_pos = Y,
    prevalence_true = p,
    ml_signal_true = ml_signal,
    S_residual_true = S_res
  ) %>%
    bind_cols(covs)
}

hybrid_ml_geo_oncho <- simulate_hybrid_ml_geo_oncho()

## ==========================================================
## Prediction grid for Nigeria (UTM + covariates)
## ==========================================================

make_prediction_grid_utm <- function(nx = 200, ny = 200,
                                     polygon_ll  = st_transform(nga_poly_utm, 4326),
                                     crs_utm = 32632) {

  # 1. Bounding box from polygon (in lon/lat)
  bb   <- st_bbox(polygon_ll)
  lon_seq <- seq(bb["xmin"], bb["xmax"], length.out = nx)
  lat_seq <- seq(bb["ymin"], bb["ymax"], length.out = ny)

  grid_ll <- expand.grid(lon = lon_seq, lat = lat_seq)
  pts_ll  <- st_as_sf(grid_ll, coords = c("lon", "lat"), crs = 4326)

  # 2. Keep only points inside Nigeria
  inside  <- st_within(pts_ll, polygon_ll, sparse = FALSE)[, 1]
  pts_ll  <- pts_ll[inside, ]
  grid_ll <- grid_ll[inside, ]

  # 3. Transform to UTM
  pts_utm <- st_transform(pts_ll, crs_utm)
  utm_coords <- st_coordinates(pts_utm)

  df <- data.frame(
    lon   = grid_ll$lon,
    lat   = grid_ll$lat,
    utm_x = utm_coords[, 1],
    utm_y = utm_coords[, 2]
  )

  # 4. Covariates for prediction grid
  covs <- simulate_covariates(df[, c("utm_x", "utm_y")])

  dplyr::bind_cols(df, covs)
}

pred_grid_oncho <- make_prediction_grid_utm(nx = 200, ny = 200)

## ==========================================================
## Save data for participants
## ==========================================================

out_dir <- "oncho_data"
if (!dir.exists(out_dir)) dir.create(out_dir)

write.csv(spatial_oncho_binom,
          file = file.path(out_dir, "spatial_oncho_binomial_data.csv"),
          row.names = FALSE)

write.csv(spatial_oncho_pois,
          file = file.path(out_dir, "spatial_oncho_poisson_data.csv"),
          row.names = FALSE)

write.csv(st_oncho_binom,
          file = file.path(out_dir, "spatiotemporal_oncho_binomial_data.csv"),
          row.names = FALSE)

write.csv(joint_oncho,
          file = file.path(out_dir, "joint_oncho_processes.csv"),
          row.names = FALSE)

write.csv(nonstat_oncho,
          file = file.path(out_dir, "nonstationary_oncho_data.csv"),
          row.names = FALSE)

write.csv(hybrid_ml_geo_oncho,
          file = file.path(out_dir, "hybrid_ml_geostatistical_oncho_data.csv"),
          row.names = FALSE)

write.csv(pred_grid_oncho,
          file = file.path(out_dir, "oncho_prediction_grid_utm.csv"),
          row.names = FALSE)

## Save Nigeria shapefile (admin-1) for mapping
nga_admin1_utm |>
  st_write(file.path(out_dir, "nga_shapefile.shp"),
           delete_layer = TRUE)

cat("Onchocerciasis datasets saved in folder:", out_dir, "\n")
cat("Files:\n",
    " - spatial_oncho_binomial_data.csv\n",
    " - spatial_oncho_poisson_data.csv\n",
    " - spatiotemporal_oncho_binomial_data.csv\n",
    " - joint_oncho_processes.csv\n",
    " - nonstationary_oncho_data.csv\n",
    " - hybrid_ml_geostatistical_oncho_data.csv\n",
    " - oncho_prediction_grid_utm.csv\n",
    " - nga_shapefile.* (Nigeria admin-1)\n")
