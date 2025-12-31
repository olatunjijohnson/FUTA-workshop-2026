# FUTA inlabru Workshop – Onchocerciasis Practice Data

This folder contains **synthetic datasets** for practicing geostatistical modelling with **inlabru** in the context of **onchocerciasis in Nigeria**.

The data are **simulated**, not real, but they are designed to mimic plausible patterns:

- Higher onchocerciasis risk in the **north of Nigeria**  
- Strong association with **rivers**, **environmental factors**, and **population density**  
- Multiple datasets tailored to different workshop topics:
  - Spatial and spatio-temporal modelling  
  - Joint modelling of multiple processes  
  - Non-stationary spatial fields  
  - Hybrid machine learning + geostatistics  

You only need these files to do the exercises; the simulation code is **not required**.

---

## 1. Coordinate Systems and Geography

Most files use **projected coordinates**:

- `utm_x`, `utm_y`:  
  - Projected coordinates in **UTM Zone 32N** (EPSG:32632).  
  - Units are **meters**.  
  - These should be used for **SPDE meshes** and spatial modelling.

The prediction grid also includes:

- `lon`, `lat`:  
  - Geographic coordinates in **WGS84** (EPSG:4326).  
  - Useful for plotting on maps or exporting to GIS.

A simplified **Nigeria admin-1 shapefile** is provided for mapping:

- `nga_shapefile.shp` + associated files (`.dbf`, `.shx`, `.prj` etc.).
- Coordinates: **UTM Zone 32N** (EPSG:32632).

---

## 2. Common Covariates (Predictors)

Many datasets share the same set of predictors:

- `elevation`  
  Approximate elevation (meters). Slight north–south gradient plus noise.

- `ndvi`  
  Synthetic **Normalized Difference Vegetation Index** (0–1, via logistic transform).  
  Roughly: greener in the south, less green in the north.

- `urban`  
  Binary indicator (0/1) for **urban location**.  
  Higher probability in some southern/eastern regions.

- `rainfall`  
  Approximate annual rainfall (mm).  
  Generally higher in the south, lower in the north.

- `temperature`  
  Approximate mean temperature (°C).  
  Generally higher in the north, lower in the south.

- `pop_density`  
  Synthetic **population density** (people per km², on a relative scale).  
  Higher in some central/southern belt; values have been trimmed to avoid extreme outliers.

- `river_index`  
  A **proximity-to-river score** (0–something), higher near a synthetic river band.  
  Onchocerciasis risk is simulated to be **strongly associated** with this variable.

- `north_index`  
  Continuous indicator of **north–south location**  
  (low in the south, high in the north).  
  Onchocerciasis risk is higher for larger `north_index`.

These covariates are designed so that:

- Onchocerciasis prevalence is **higher in northern Nigeria**,  
- And **strongly related** to rivers and environmental conditions.

---

## 3. Data Files

### 3.1 `spatial_oncho_binomial_data.csv`

**Purpose:**  
Spatial **binomial prevalence** data (e.g., community surveys at locations across Nigeria).  
Suitable for **Day 1 spatial binomial** practice.

**Rows:**  
One row per survey location.

**Key columns:**

- `id`  
  Survey location ID.

- `utm_x`, `utm_y`  
  Coordinates in UTM Zone 32N (meters).

- `n_examined`  
  Number of individuals examined at the location.

- `n_pos`  
  Number of **onchocerciasis-positive** individuals.

- `prevalence_true`  
  The simulated true underlying prevalence (probability).  
  > Use this mainly for **validation / checking**. In typical modelling, you would **not** use this column.

- `S`  
  Latent spatial Gaussian field value at the location.  
  > Again, provided for **diagnostics or method evaluation**, not for fitting.

- Covariates:  
  - `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
    `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Fit a **spatial binomial model**:  
  `n_pos` / `n_examined` ~ covariates + SPDE spatial field.  
- Explore **mesh design**, **SPDE priors**, and **prediction surfaces**.

---

### 3.2 `spatial_oncho_poisson_data.csv`

**Purpose:**  
Spatial **Poisson count** data (e.g., counts of nodules or positive skin snips).  
Suitable for **Day 1 Poisson** practice.

**Rows:**  
One row per location.

**Key columns:**

- `id`  
  Location ID.

- `utm_x`, `utm_y`  
  UTM coordinates (meters).

- `population`  
  Local population at risk (for use as an **offset**).

- `cases`  
  Number of onchocerciasis cases / nodules / positive skin snips.

- `lambda_true`  
  Underlying expected count at the location (population × true rate).  
  > For validation / checking.

- `log_rate_true`  
  Underlying log incidence rate (before multiplying by `population`).

- `S`  
  Latent spatial field for the log-rate.

- Covariates:  
  `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
  `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Fit a **Poisson spatial model**:  
  `cases` ~ covariates + SPDE field + offset(log(population)).

---

### 3.3 `spatiotemporal_oncho_binomial_data.csv`

**Purpose:**  
Spatio-temporal **binomial prevalence** data (e.g., repeated surveys over years).  
Suitable for **Day 1 spatio-temporal** practice.

**Rows:**  
One row per (site, time) combination.

**Key columns:**

- `site`  
  Site ID (spatial location).

- `time`  
  Time index (e.g., 1–6).  
  Interpretable as years, seasons, or survey rounds.

- `utm_x`, `utm_y`  
  Coordinates of the site in UTM.

- `S`  
  Spatial field value at the site (constant over time).

- `gamma_t`  
  Temporal random effect at time `t` (common to all sites at that time).

- `n_examined`  
  Number of individuals examined at that site and time.

- `n_pos`  
  Number of positives at that site and time.

- `prevalence_true`  
  True underlying probability for that site and time.

- Covariates:  
  `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
  `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Fit **separable** or **inseparable** spatio-temporal models using:
  - Spatial SPDE field,
  - Temporal random effect (e.g., RW1/RW2),
  - Or structured space–time interactions.

---

### 3.4 `joint_oncho_processes.csv`

**Purpose:**  
Data for **joint modelling of multiple onchocerciasis processes**, e.g.:

- Different **risk groups** (adults vs children), or  
- Different **diagnostic methods**.

Suitable for **Day 2 joint modelling** practice.

**Rows:**  
One row per (location, group).

**Key columns:**

- `id`  
  Location ID.

- `group`  
  Group label (e.g., `"adults"`, `"children_or_alt_diag"`).

- `utm_x`, `utm_y`  
  UTM coordinates.

- `n_examined`  
  Number examined in that group at the location.

- `n_pos`  
  Number positive in that group.

- `prevalence_true`  
  True underlying prevalence for that group at that location.

- `S_shared`  
  Shared latent spatial field, common to both groups.

- `S_specific`  
  Group-specific latent spatial field.

- Covariates:  
  `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
  `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Fit **joint binomial models** with:
  - Shared spatial field,
  - Optional group-specific spatial components,
  - Group-specific fixed effects.

---

### 3.5 `nonstationary_oncho_data.csv`

**Purpose:**  
Data where the spatial effect is **non-stationary**, e.g.:

- Smoother in the **north**, rougher in the **south**.  

Suitable for **Day 2 non-stationarity** practice.

**Rows:**  
One row per location.

**Key columns:**

- `id`  
  Location ID.

- `utm_x`, `utm_y`  
  UTM coordinates.

- `n_examined`  
  Number examined.

- `n_pos`  
  Number positive.

- `prevalence_true`  
  True underlying prevalence.

- `S_nonstationary`  
  Combined non-stationary spatial field used in the simulation.

- `weight_north`  
  Weight (0–1) indicating how “north” a location is (low in south, high in north).  
  Used internally to mix two spatial fields.

- `S_smooth`  
  Smooth field (longer range).

- `S_rough`  
  Rough field (shorter range).

- Covariates:  
  `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
  `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Compare **stationary** vs **non-stationary** spatial models.  
- Explore models with:
  - Spatially varying parameters, or  
  - Mixtures of fields in different regions.

---

### 3.6 `hybrid_ml_geostatistical_oncho_data.csv`

**Purpose:**  
Data with a complex **machine learning-style signal** plus a residual spatial field.

- The outcome depends on **non-linear combinations** of covariates,
  interactions, and a residual spatial effect.
- Designed for **Day 3 hybrid ML + geostatistics** practice.

**Rows:**  
One row per survey location.

**Key columns:**

- `id`  
  Location ID.

- `utm_x`, `utm_y`  
  UTM coordinates.

- `n_examined`  
  Number examined.

- `n_pos`  
  Number positive.

- `prevalence_true`  
  True underlying probability.

- `ml_signal_true`  
  Deterministic signal from a complex (non-linear) function of the covariates.  
  Think: “what a powerful ML model might capture”.

- `S_residual_true`  
  Residual spatial field not explained by the ML signal.

- Covariates:  
  `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
  `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Fit **ML models** (e.g., Random Forest, XGBoost) using covariates only.  
- Examine **spatial autocorrelation of residuals**.  
- Build **hybrid models**:
  - ML part for complex covariate effects,
  - inlabru spatial field for residual spatial structure.

---

### 3.7 `oncho_prediction_grid_utm.csv`

**Purpose:**  
Prediction grid over Nigeria with covariates, for mapping model-based predictions.

**Rows:**  
Prediction locations inside the Nigeria polygon, typically on a regular grid.

**Key columns:**

- `lon`, `lat`  
  Longitude/latitude (WGS84).

- `utm_x`, `utm_y`  
  UTM coordinates (EPSG:32632).

- Covariates:  
  `elevation`, `ndvi`, `urban`, `rainfall`, `temperature`,  
  `pop_density`, `river_index`, `north_index`.

**Typical use:**

- Build a mesh over `utm_x`, `utm_y`.  
- Use as the `pred` / `integration` locations for:
  - Predicted prevalence,  
  - Exceedance probabilities,  
  - Mapping spatial fields.

---

### 3.8 `nga_shapefile.*`

**Purpose:**  
Simple **Nigeria admin-1** boundaries for maps.

**Files:**

- `nga_shapefile.shp`  
- `nga_shapefile.dbf`  
- `nga_shapefile.shx`  
- `nga_shapefile.prj`  
- ... (other standard shapefile components)

**Coordinate reference:**  
UTM Zone 32N (EPSG:32632).

**Typical use:**

- Plot boundaries under your prediction maps.  
- Clip or mask predictions to the country boundary.

---

## 4. Example: Loading Data in R

```r
# Set working directory to the folder containing the data/
setwd("path/to/project")

library(readr)
library(dplyr)

# 1. Spatial binomial onchocerciasis data
spatial_oncho <- read_csv("oncho_data/spatial_oncho_binomial_data.csv")

# 2. Spatio-temporal data
st_oncho <- read_csv("oncho_data/spatiotemporal_oncho_binomial_data.csv")

# 3. Prediction grid
pred_grid <- read_csv("oncho_data/oncho_prediction_grid_utm.csv")
