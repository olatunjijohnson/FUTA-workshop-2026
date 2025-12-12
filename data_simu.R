## ==========================================================
## Malaria geostatistical simulations for inlabru workshop
## Topics:
## 1. Spatial + spatio-temporal, different outcome types
## 2. Joint modelling of multiple processes
## 3. Non-stationary spatial process
## 4. Hybrid ML + geostatistical model
## ==========================================================

## ---------- Packages ----------
library(mvtnorm)   # for multivariate Gaussian simulation
library(dplyr)
library(tidyr)
library(sf)
library(tidyverse)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)
# remotes::install_github("inlabru-org/inlabru", ref = "devel")
library(inlabru)

set.seed(1234)


library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

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

## Some synthetic environmental covariates
simulate_covariates <- function(coords) {
  # coords is a data.frame with utm_x, utm_y
  x <- scale(coords[, "utm_x"])
  y <- scale(coords[, "utm_y"])
  n <- nrow(coords)

  elevation <- 300 + 60 * sin(pi * y) + rnorm(n, sd = 20)   # gentle N–S gradient
  ndvi      <- plogis(-0.5 + 0.8 * y + 0.4 * sin(pi * x))  # smoother north–south pattern
  urban_p   <- plogis(-0.4 + 0.5 * x)                      # more urban east, say
  urban     <- rbinom(n, 1, urban_p)

  data.frame(
    elevation = elevation,
    ndvi      = ndvi,
    urban     = urban
  )
}


## ==========================================================
## 1. Spatial & spatio-temporal malaria data (different outcomes)
## ==========================================================

## ---------- 1A. Spatial Binomial prevalence ----------
simulate_spatial_binomial <- function(n_sites = 250) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Latent spatial field
  S <- simulate_spatial_field(coords, range = 400000, sigma2 = 2.7)

  ## Linear predictor (logit-prevalence)
  beta0 <- -1.2   # baseline prevalence ~ 23%
  beta_elev <- -0.001
  beta_ndvi <- 1.0
  beta_urban <- 0.4

  linpred <- beta0 +
    beta_elev * covs$elevation +
    beta_ndvi * covs$ndvi +
    beta_urban * covs$urban +
    S

  p <- plogis(linpred)

  N <- sample(50:200, n_sites, replace = TRUE)  # children tested
  Y <- rbinom(n_sites, size = N, prob = p)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_tested = N,
    n_pos = Y,
    prevalence_true = p,
    S = S
  ) %>%
    bind_cols(covs)
}

spatial_binom <- simulate_spatial_binomial()


## ---------- 1B. Spatial Poisson counts (case counts) ----------
simulate_spatial_poisson <- function(n_sites = 250) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Latent spatial field shared with prevalence or independent
  S <- simulate_spatial_field(coords, range = 500000, sigma2 = 0.5)

  ## Population at risk (offset)
  pop <- round(runif(n_sites, 200, 3000))

  ## Log incidence rate
  alpha0 <- -6.0           # baseline ~ few cases per year per person
  alpha_elev <- -0.0005
  alpha_ndvi <- 0.8
  alpha_urban <- 0.6

  log_rate <- alpha0 +
    alpha_elev * covs$elevation +
    alpha_ndvi * covs$ndvi +
    alpha_urban * covs$urban +
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

spatial_pois <- simulate_spatial_poisson()
head(spatial_pois)

## ---------- 1C. Spatio-temporal Binomial prevalence ----------
simulate_spatiotemporal_binomial <- function(n_sites = 200, T = 6) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Spatial field
  S <- simulate_spatial_field(coords, range = 400000, sigma2 = 0.6)

  ## Temporal AR(1) (e.g. 6 annual surveys, or 6 seasons)
  gamma_t <- simulate_AR1(T = T, rho = 0.7, sigma2 = 0.3)

  beta0 <- -1.2
  beta_elev <- -0.001
  beta_ndvi <- 0.9
  beta_urban <- 0.5

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
    beta_elev * df$elevation +
    beta_ndvi * df$ndvi +
    beta_urban * df$urban +
    df$S +
    df$gamma_t

  p <- plogis(linpred)
  N <- sample(50:200, nrow(df), replace = TRUE)
  Y <- rbinom(nrow(df), size = N, prob = p)

  df$n_tested <- N
  df$n_pos <- Y
  df$prevalence_true <- p

  df
}

st_binom <- simulate_spatiotemporal_binomial()

## ==========================================================
## 2. Joint modelling of multiple malaria processes
##    (e.g. two risk groups: children vs pregnant women)
## ==========================================================

simulate_joint_malaria_processes <- function(n_sites = 250) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Shared spatial field
  S_shared <- simulate_spatial_field(coords, range = 500000, sigma2 = 0.7)

  ## Outcome-specific small-scale fields
  S_child  <- simulate_spatial_field(coords, range = 200000, sigma2 = 0.2)
  S_preg   <- simulate_spatial_field(coords, range = 200000, sigma2 = 0.2)

  ## Children <5
  beta0_c <- -1.0
  linpred_c <- beta0_c +
    0.9 * covs$ndvi -
    0.001 * covs$elevation +
    0.4 * covs$urban +
    S_shared + S_child
  p_c <- plogis(linpred_c)
  N_c <- sample(50:200, n_sites, replace = TRUE)
  Y_c <- rbinom(n_sites, N_c, p_c)

  ## Pregnant women
  beta0_p <- -0.5     # typically different baseline
  linpred_p <- beta0_p +
    0.7 * covs$ndvi -
    0.0008 * covs$elevation +
    0.2 * covs$urban +
    S_shared + S_preg
  p_p <- plogis(linpred_p)
  N_p <- sample(20:100, n_sites, replace = TRUE)
  Y_p <- rbinom(n_sites, N_p, p_p)

  df_children <- data.frame(
    id = 1:n_sites,
    group = "children_u5",
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_tested = N_c,
    n_pos = Y_c,
    prevalence_true = p_c,
    S_shared = S_shared,
    S_specific = S_child
  ) %>%
    bind_cols(covs)

  df_preg <- data.frame(
    id = 1:n_sites,
    group = "pregnant_women",
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_tested = N_p,
    n_pos = Y_p,
    prevalence_true = p_p,
    S_shared = S_shared,
    S_specific = S_preg
  ) %>%
    bind_cols(covs)

  bind_rows(df_children, df_preg)
}

joint_malaria <- simulate_joint_malaria_processes()

## ==========================================================
## 3. Non-stationary spatial process
##    (e.g. smoother field in savannah, rougher in forest belt)
## ==========================================================

simulate_nonstationary_malaria <- function(n_sites = 300) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## Two spatial fields with different smoothness / range
  S_smooth <- simulate_spatial_field(coords, range = 600000, sigma2 = 0.5)
  S_rough  <- simulate_spatial_field(coords, range = 200000, sigma2 = 0.5)

  ## Spatially varying mixing weight (e.g. based on latitude)
  ## More "rough" in the south, smoother in the north
  utm_y <- coords[, "utm_y"]
  w <- (utm_y - min(utm_y)) / (max(utm_y) - min(utm_y))  # 0 in south, 1 in north

  S_ns <- sqrt(w) * S_smooth + sqrt(1 - w) * S_rough

  ## Malaria prevalence with non-stationary spatial effect
  beta0 <- -1.0
  linpred <- beta0 +
    0.8 * covs$ndvi -
    0.001 * covs$elevation +
    0.3 * covs$urban +
    S_ns

  p <- plogis(linpred)
  N <- sample(50:200, n_sites, replace = TRUE)
  Y <- rbinom(n_sites, N, p)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_tested = N,
    n_pos = Y,
    prevalence_true = p,
    S_nonstationary = S_ns,
    weight_north = w,
    S_smooth = S_smooth,
    S_rough = S_rough
  ) %>%
    bind_cols(covs)
}

nonstat_malaria <- simulate_nonstationary_malaria()

## ==========================================================
## 4. Hybrid ML + geostatistical malaria data
##    (complex covariate effects + residual spatial field)
## ==========================================================

simulate_hybrid_ml_geo <- function(n_sites = 400) {
  coords <- sample_nigeria_utm_only(n_sites)
  covs   <- simulate_covariates(coords)

  ## "True" complex ML-type signal (non-linear, interactions)
  ## Think of this as what a random forest / XGBoost might learn.
  ml_signal <- with(covs, {
    f1 <- 2 * (ndvi - 0.6)^2 * (-1)   # peak risk around NDVI ~ 0.6
    f2 <- 0.5 * (elevation / 100)^2 * (-1) # lower elevation -> higher risk
    f3 <- ifelse(urban == 1, 0.5, -0.2)
    f4 <- 0.7 * sin(coords[, "utm_x"] / 2)  # east-west oscillation
    f1 + f2 + f3 + f4
  })

  ## Residual spatial field (to be modelled by inlabru)
  S_res <- simulate_spatial_field(coords, range = 400000, sigma2 = 0.3)

  beta0 <- -1.2
  linpred <- beta0 + ml_signal + S_res
  p <- plogis(linpred)

  N <- sample(50:200, n_sites, replace = TRUE)
  Y <- rbinom(n_sites, N, p)

  data.frame(
    id = 1:n_sites,
    utm_x = coords[, "utm_x"],
    utm_y = coords[, "utm_y"],
    n_tested = N,
    n_pos = Y,
    prevalence_true = p,
    ml_signal_true = ml_signal,
    S_residual_true = S_res
  ) %>%
    bind_cols(covs)
}

hybrid_ml_geo <- simulate_hybrid_ml_geo()

## ==========================================================
## Objects now available:
##   spatial_binom   : Topic 1 (spatial binomial)
##   spatial_pois    : Topic 1 (spatial Poisson)
##   st_binom        : Topic 1 (spatio-temporal binomial)
##   joint_malaria   : Topic 2 (joint processes)
##   nonstat_malaria : Topic 3 (non-stationary)
##   hybrid_ml_geo   : Topic 4 (ML + geostatistical)
## ==========================================================



###################### Predictions #######################


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

  # 4. Smooth covariates
  covs <- simulate_covariates(df[, c("utm_x", "utm_y")])

  dplyr::bind_cols(df, covs)
}

pred_grid <- make_prediction_grid_utm(nx = 200, ny = 200)



head(pred_grid)

########################### Save the data ###################

# -----------------------------------------------------------
# Create an output folder for the workshop data
# -----------------------------------------------------------
out_dir <- "data"
if (!dir.exists(out_dir)) dir.create(out_dir)

# -----------------------------------------------------------
# Save each dataset as CSV
# -----------------------------------------------------------

write.csv(spatial_binom,
          file = file.path(out_dir, "spatial_binomial_data.csv"),
          row.names = FALSE)

write.csv(spatial_pois,
          file = file.path(out_dir, "spatial_poisson_data.csv"),
          row.names = FALSE)

write.csv(st_binom,
          file = file.path(out_dir, "spatiotemporal_binomial_data.csv"),
          row.names = FALSE)

write.csv(joint_malaria,
          file = file.path(out_dir, "joint_malaria_processes.csv"),
          row.names = FALSE)

write.csv(nonstat_malaria,
          file = file.path(out_dir, "nonstationary_malaria_data.csv"),
          row.names = FALSE)

write.csv(hybrid_ml_geo,
          file = file.path(out_dir, "hybrid_ml_geostatistical_data.csv"),
          row.names = FALSE)

# -----------------------------------------------------------
# Save prediction grid (UTM + covariates)
# -----------------------------------------------------------
write.csv(pred_grid,
          file = file.path(out_dir, "prediction_grid_utm.csv"),
          row.names = FALSE)

# -----------------------------------------------------------
# Confirmation message
# -----------------------------------------------------------
cat("All datasets saved in folder:", out_dir, "\n")
