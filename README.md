# A repository for:

Kase, A et al. The impact of local and long-distance colonization in a fragmented landscape

## Links to different parts of the readme file

1. [What is in this repository?](#what-is-in-this-repository)
2. [The working directory](#the-working-directory)
3. [The data folder (`./data`)](#the-data-folder-data)
4. [The R folder (`./R`)](#the-r-folder-R)
5. [The nimble folder (`./nimble`)](#the-nimble-folder-nimble)
6. [The fuzzy plots folder (`./fuzzy_plots`)](#the-fuzzy-plots-folder-fuzzy_plots)
7. [The plots folder (`./plots`)](#the-plots-folder-plots)
8. [The results folder (`./skunk_rds`)](#the-results-folder-skunk_rds)

## What is in this repository?
This repository stores all of the data and code used to fit a dynamic occupancy model with explicit colonization terms, compare model outputs using Brier scores, and forecast model estimates across space and time. The folder organization separates the data (`./data`), figures from the manuscript (`./figures`), nimble models (`./nimble`), and the R code (`./R`).

This document describes all of the files present in the repository.

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

## The working directory
Aside from the aforementioned folders, the working directory here stores the `.gitignore` file for this repository, this README file (`README.md`) and the `.Rproj` file (for if you are using RStudio, `skunks.Rproj`).

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

## The data folder (`./data`)
This folder has 6 files and 1 sub-folder (`./data/raw_sim_covariates`). 

**| `./data/complete_data.csv` |** The camera trap data used in our analysis. This csv file has 3,075 rows and 9 columns.

| Column  | Data Type | Explanation                                                                                                                                                                                                  |
| ------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Species | `factor`    | The species associated to this data point                                                                                                                                                                    |
| Season  | `factor`    | A seasonal code for the season the data comes from. It combines the first two letters of the season and the last two digits of the year. Seasonal codes are SP = Spring, SU = Summer, FA = Fall, WI = Winter |
| Site    | `factor`    | The site associated to this data point                                                                                                                                                                       |
| City    | `factor`    | The city abbreviation associated to this data point                                                                                                                                                          |
| Long    | `numeric`   | The longitude of the site associated to this data point (WGS 84)                                                                                                                                             |
| Lat     | `numeric`   | The latitude of the site associated to this data point (WGS 84)                                                                                                                                              |
| Crs     | `integer`   | The coordinate reference system code for the site coordinates                                                                                                                                                |
| Y       | `integer`   | The number of days the species was detected during a given season                                                                                                                                            |
| J       | `integer`   | The number of days the camera trap was operational during a given season. If `J == 0` then no sampling occurred                                                                                                |


**| `./data/site_covariates.csv` |** The spatial covariate (urbanization score, distance to water, and proportion of managed open space) values for each site. This csv file is generated from the [`creating_spatial_points.R`, `dist_water_cov.R`, `managed_lawn_cov.R`, and `scaled_covariates.R`](#group 2---preparing-spatial-covariates) scripts located in the `./R` folder. This csv file has 107 rows and 7 columns.

| Column     | Data Type | Explanation                                                                                                                                                                     |
| ---------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Site       | `factor`  | The site associated to this data point                                                                                                                                          |
| HU10       | `numeric` | The mean housing density per meter<sup>2</sup> within a 1 kilometer radius of a site                                                                                            |
| tree       | `numeric` | The mean proportion of tree canopy cover within a 1 kilometer radius of a site                                                                                                  |
| imperv     | `numeric` | The mean proportion of impervious cover within a 1 kilometer radius of a site                                                                                                   |
| urb        | `numeric` | Urban intensity metric created from a principal component analysis using mean housing density, mean tree cover, and mean impervious cover within a 1 kilometer radius of a site |
| water_dist | `numeric` | The shortest Euclidian distance between a site and a permanent body of water in meters                                                                                          |
| open_dev   | `numeric` | The mean proportion of developed open space (areas in which impervious surfaces account for less than 20 percent of total cover) within a 1 kilometer radius of a site          |


**| `./data/original_z_sim_values.csv` |** The simulated initial occupancy values for each site across the Chicagoland study area used to forecast striped skunk occupancy predictions through space and time. This csv file has 4067 rows and 13 columns. Each row is a site, and each column is a timestep.

**| `./data/scaled_simulation_covariates.csv` |** The spatial covariate data for each site across the Chicagoland study area used to forecast striped skunk occupancy predictions through space and time. This file is generated by the [`./R/skunk_simulation_covariate_prep.R`](#group-6---simulations) located in the `./R` folder. This csv file contains 4067 rows and 7 columns.

| Column     | Data Type | Explanation                                                                                                                                                                     |
| ---------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Site       | `factor`  | The simulated site associated to this data point                                                                                                                                |
| x          | `integer` | Horizontal/Latitudinal location of simulated site in landscape grid                                                                                                             |
| y          | `integer` | Vertical/Longitudinal location of simulated site in landscape grid                                                                                                              |
| crs        | `integer` | The coordinate reference system code for the site coordinates                                                                                                                   |
| urb        | `numeric` | Urban intensity metric created from a principal component analysis using mean housing density, mean tree cover, and mean impervious cover within a 1 kilometer radius of a site |
| water_dist | `numeric` | The shortest Euclidian distance between a site and a permanent body of water in meters                                                                                          |
| open_dev   | `numeric` | The mean proportion of developed open space (areas in which impervious surfaces account for less than 20 percent of total cover) within a 1 kilometer radius of a site          |


**| `./data/phi_gamma_delta_by_season.RDS` |** The output of a model to simulate occupancy and its parameters across Chicagoland area and through time created by [`./R/simulation_model.R`](#group-6---simulations) in the `./R` folder.

**| `./data/simulation_neighbors.RDS` |** The output of a simulated model used to calculate the number of occupied neighbors in the previous timestep for each simulated site created by [`./R/skunk_simulation_covariate_prep.R`](#group-6---simulations) in the `./R` folder.

The `./data/raw_sim_covariates` sub-folder contains 4 csv files. 

**| `./data/raw_sim_covariates/dist2water.csv` |** The raw shortest Euclidian distance between a site and a permanent body of water for each simulated site created by [`./R/skunk_simulation_query_covars.R`](#group-6---simulations) in the `./R` folder. This file contains 4067 rows and 2 columns.

| Column     | Data Type | Explanation                                                                                |
| ---------- | --------- | ------------------------------------------------------------------------------------------ |
| Site       | `factor`  | The simulated site associated to this data point                                           |
| water_dist | `numeric` | The raw shortest Euclidian distance between a site and a permanent body of water in meters |


**| `./data/raw_sim_covariates/point_locs.csv` |** The locations of simulated sites in landscape grid created by [`./R/skunk_simulation_query_covars.R`](#group-6---simulations) in the `./R` folder. This file contains 4067 rows and 4 columns.

| Column | Data Type | Explanation                                                         |
| ------ | --------- | ------------------------------------------------------------------- |
| Site   | `factor`  | The simulated site associated to this data point                    |
| x      | `integer` | Horizontal/Latitudinal location of simulated site in landscape grid |
| y      | `integer` | Vertical/Longitudinal location of simulated site in landscape grid  |
| crs    | `integer` | The coordinate reference system code for the site coordinates       |


**| `./data/raw_sim_covariates/urb_covars.csv` |** The 

**| `./data/raw_sim_covariates/urban_openspace.csv` |**


[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

## The R folder (`./R`)

This folder contains 35 R scripts. We have grouped these scripts based upon the function they perform in the overall scheme of the analysis and the relative order in which they would be run if someone was interested in recreating this analysis.

### **Group 1 - Functions**  
These scripts contain various functions that are sourced in other scripts to complete the analysis


| File                          | Description                                                                  | Packages Required     |
| ----------------------------- | ---------------------------------------------------------------------------- | --------------------- |
| **`./R/init_functions.R`**    | Functions to assign initial values to parameters during model run            | Only Base R used      |
| **`./R.mcmc_functions.R`**    | Functions to extract the MCMC posterior from raw model RDS outputs           | Only Base R used      |
| **`./R/raster_extraction.R`** | Functions to extract spatial data from raster files and calculate proportion | `cli`, `raster`, `sf` |


### **Group 2 - Preparing Spatial Covariates**  
These scripts source and extract spatial data, calculate covariate values for each site, and compile data into the `./data/site_covariates.csv`


| File                                | Description                                                                                                                                                                                                                     | Packages Required                           |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| **`./R/creating_spatial_points.R`** | Creates spatial points for each site and projects them into UTM                                                                                                                                                                 | `dplyr`, `sf`                               |
| **`./R/dist_water_cov.R`**          | Generates the shortest Euclidian distances between each site and a permanent body of water                                                                                                                                      | `dplyr`, `sf`                               |
| **`./R/managed_lawn_cov.R`**        | Generates the mean proportion of developed open space (areas in which impervious surfaces account for less than 20 percent of total cover) within a 1 kilometer radius of a site                                                | `dplyr`, `sf`                               |
| **`./R/scaled_covariates.R`**       | Extracts the mean housing density, proportion of tree canopy cover, and the proportion of developed open space, and generates the urban intensity metric from a principal component analysis using the aforementioned variables | `dplyr`, `raster`, `sf`, `uwinspatialtools` |

**Note:** `uwinspatialtools` is an R package deleveloped by Dr. Mason Fidino, which essentially has some wrapper functions for `sf` and `raster`. It can be found at www.github/com/mfidino/uwinspatialtools.


### **Group 3 - Data Prep**

These scripts prepare the detection data and appropriate spatial and temporal covariates to be input into the model. The list below follows the same order as our hypotheses in Table 1 of the manuscript. 

**`./R/intercept_data_prep.R`** 

**`./R/spatial_covariates_data_prep.R`** 

**`./R/fall_term_data_prep.R`** 

**`./R/spatial_covariates_fall_data_prep.R`** 

**`./R/spatial_covariates_inx_data_prep.R`** 

**`./R/spatial_covariates_fall_urbless_data_prep.R`** 


These scripts require the R packages `dplyr`, and `sf`.


### **Group 4 - Model Run**

These scripts run the models by sourcing the appropriate data prep, initial values, and nimble scripts (from the `./nimble` folder), and saving the outputs as RDS files into the `./skunk_rds` folder. The end of these scripts also include a visual check of model convergence by plotting the MCMC chains and saving the plots into the `./fuzzy_plots` folder. The list below follows the same order as our hypotheses in Table 1 of the manuscript.

**`./R/intercept_model_run.R`** 

**`./R/spatial_covariates_model_run.R`** 

**`./R/fall_term_model_run.R`** 

**`./R/spatial_covariates_fall_model_run.R`** 

**`./R/spatial_covariates_inx_model_run.R`** 

**`./R/spatial_covariates_fall_urbless_model_run.R`** 


These scipts require the R packages `dplyr`, `MCMCvis`, `nimble`, `parallel`, and `scales`.


### **Group 5 - Brier Scores**

These scripts calculate the Brier score to be used in model selection for each model from the model output RDS files and out of sample detection data.  

**Note:** The `calculate_all_briers.R` is the only script that needs to be run as it sources the individual model Brier score calculation scripts listed below in the same order as our hypotheses in Table 1 of the manuscript.

**`./R/intercept_brier_score.R`** 

**`./R/spatial_covariates_brier_score.R`** 

**`./R/fall_term_brier_score.R`** 

**`./R/spatial_covariates_fall_brier_score.R`** 

**`./R/spatial_covariates_urbless_brier_score.R`** 

**`./R/spatial_covariates_fall_urbless_brier_score.R`** 


These scripts require the R packages `dplyr`, `MCMCvis`, and `nimble`.


### **Group 6 - Simulations**

These scripts are used to simulated striped skunk occupancy across the Chicagoland study area from the output of the best predicting model, spatial covariates and species life history (Hypothesis 3 from Table 1 in manuscript).


| File                                        | Description                                                                                                                                                                                                                                                                                                                                                                                                 | Required Packages                           |
| ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| **`./R/skunk_simulation_query_covars.R`**   | Creates a grid covering the Chicagoland area and extracts spatial covariate data from original data sources fore each grid cell rather than individual study sites. Grid cells and spatial covariate data are written into the following csv files and saved into the `raw_sim_covariates` sub-folder in the `./data` folder: `point_locs.csv`, `urb_covars.csv`, `dist2water.csv`,  `urban_openspace.csv`. | `dplyr`, `raster`, `sf`, `uwinspatialtools` |
| **`./R/skunk_simulation_covariate_prep.R`** | Prepares spatial covariate data and number of neighboring sites data for each site to be used in the simulated model. Spatial covariate data are written into a csv file (\`./data/scaled_simulation_covariates.csv\`), and neighbor data are written into an RDS file (\`./data/simulation_neighbors.RDS\`).                                                                                               | `dplyr`, `sf`                               |
| **`./R/simulation_model.R`**                | Runs a model to simulate occupancy and its parameters across Chicagoland area and through time. The output is saved into the `./data` folder as and RDS file (\`phi_gamma_delta_by_season.RDS\`).                                                                                                                                                                                                           | `dplyr`                                     |
| **`./R/occupancy_simulation_example.R`**    | This is an example and test file for simulating occupancy through time and space using our model parameters.                                                                                                                                                                                                                                                                                                | Only Base R used                            |

### **Group 7 - Plotting**

These scripts are used to plot occupancy and its processes, and ultimately create the figures in the manuscrpt. 
All output figures are saved into the `./plots` folder.

| File                                             | Description                                                                                                                                                                      | Required Packages            |
| ------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------- |
| **`./R/random_colonization_plotting.R`**       | Plots random colonization probabilities (γ<sub>i,t</sub>) as a function of the three spatial covariates (urbanization, distance to water, and developed open space), and season. | `bbplot`, `dplyr`, `MCMCvis` |
| **`./R/neighborhood_colonization_plotting.R`** | Plots the neighborhood colonization probabilities (d<sub>i,t</sub>) as function of how many occupied neighboring sites were present in the previous timestep.                    | `bbplot`, `dplyr`, `MCMCvis` |
| **`./R/persistence.R`**                        | Plots the probability of persistence (ϕ<sub>i,t</sub>) as a function of the three spatial covariates (urbanization, distance to water, and developed open space).                | `bbplot`, `dplyr`, `MCMCvis` |
| **`./R/simulation_map_plots.R`**               | Plots the simulated occupancy probability over the Chicagoland study area.                                                                                                       | `bbplot`, `prettymapr`, `sf` |
| **`./R/occupancy.R`**                          | Plots the overall occupancy probability as a function of the three spatial covariates (urbanization, distance to water, and developed open space).                               | `bbplot`, `dplyr`, `MCMCvis` |

**Note:** `bbplot` is a plotting package developed by Dr. Juniper Simonis that Dr. Mason Fidino has collaborated on. It can be found at https://github.com/dapperstats/bbplot.

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)


## The nimble folder (`./nimble`)

This folder contains the 6 `nimble` scripts used to fit each model. 

| **`./nimble/intercept_model.R`** | The intercept only model (no spatial or temporal predictors).

| **`./nimble/spatial_covariates_model.R`** | The model parameterized by the three site specific spatial covariates (urbanization metric, distance to water, proportion developed open space; H1 in Table 1 in manuscript).

| **`./nimble/fall_term_model.R`** | The model including a temporally varying dummy variable to account for striped skunk life history (H2 in Table 1 in manuscript).

| **`./nimble/spatial_covariates_fall_model.R`** | The model parameterized by the three site specific spatial covariates and the temporally varying dummy variable (H3 in Table 1 in manuscript).

| **`./nimble/spatial_covariates_inxs_model.R`** | The model parameterized by the three site specific spatial covariates and an additional term to evaluate whether striped skunk have become more or less urban through time (H4 in Table 1 in manuscript).

| **`./nimble/spatial_covariates_fall_urbless_model.R`** | The model parameterized by the three site specific spatial covariates, the temporally varying dummy variable, and the additional term to evaluate whether striped skunk have become more or less urban through time (H5 in Table 1 in manuscript).

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)


## The fuzzy plots folder (`./fuzzy_plots`)

This folder is intentionally left blank, and stored the temporary files for the traceplots of model parameters that allow us to visually assess MCMC chain convergence after model runs are complete.


[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)


## The plots folder (`./plots`)

This folder contains some of the raw figures (which were later cleaned up in Inkscape), as well as other figures that were publication ready. All of the extra "cleaning" was related to the maps and ultimately the creation of Figure 4.

| File                                     | Description                                                           |
| ---------------------------------------- | --------------------------------------------------------------------- |
| **`./plots/bivariate.svg`**            | Panel C of Figure 4 in manuscript.                                    |
| **`./plots/colonization_changes.svg`** | Panel D of Figure 4 in manuscript.                                    |
| **`./plots/county_reference.svg`**     | Overlay of Chicagoland county borders used for Figure 4 map creation. |
| **`./plots/delta_figure.tiff`**        | Figure 3 in manuscript.                                               |
| **`./plots/gamma_figure.tiff`**        | Figure 2 in manuscript.                                               |
| **`./plots/mean_occ_legend.svg`**      | Legend for mean occupancy map from Panel A of Figure 4 in manuscript. |
| **`./plots/mean_occ_map.svg`**         | Only the map from Panel A of Figure 4 in manuscript.                  |
| **`./plots/mean_occ_map.tiff`**        | Raw Panel A of Figure 4 in manuscript.                                |
| **`./plots/mean_occ.svg`**             | Panel A of Figure 4 in manuscript.                                    |
| **`./plots/occ_maps.pdf`**             | Figure 4 in manuscript as a PDF file.                                 |
| **`./plots/occ_maps.png`**             | Figure 4 in manuscripts as a png file.                                |
| **`./plots/occ_maps.tiff`**            | Figure 4 in manuscripts as a tiff file.                               |
| **`./plots/phi_figure.tiff`**          | Figure 1 in manuscript.                                               |
| **`./plots/sd_occ_map.svg`**           | Only the map from Panel B of Figure 4 in manuscript.                  |
| **`./plots/sd_occ_map.tiff`**          | Raw Panel B of Figure 4 in manuscript.                                |
| **`./plots/sd_occ.svg`**               | Panel B of Figure 4 in manuscript.                                    |


[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)


## The results folder (`./skunk_rds`)

This folder contains the 4 RDS model outputs.

|**`intercept_only.RDS`**| The output of the intercept only model.

|**`fall_term.RDS`**| The output of the model including a temporally varying dummy variable to account for striped skunk life history (H2 in Table 1 in manuscript).

|**`spatial_covariates_fall2.RDS`**| The output of the model parameterized by the three site specific spatial covariates and the temporally varying dummy variable (H3 in Table 1 in manuscript).

|**`spatial_covariates_fall_urbless.RDS`**| The output of the model parameterized by the three site specific spatial covariates, the temporally varying dummy variable, and the additional term to evaluate whether striped skunk have become more or less urban through time (H5 in Table 1 in manuscript).

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

