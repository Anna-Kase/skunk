# A repository for:

Kase, A et al. The impact of local and long-distance colonization in a fragmented landscape

## Links to different parts of the readme file

1. [What is in this repository?](#what-is-in-this-repository)
2. [The working directory](#the-working-directory)
3. [The data folder (`./data`)](#the-data-folder-data)
4. [The R folder (`./R`)](#the-r-folder-R)

## What is in this repository?
This repository stores all of the data and code used to fit a dynamic occupancy model with explicit colonization terms, compare model outputs using Brier scores, and forecast model estimates across space and time. The folder organization separates the data (`./data`), figures from the manuscript (`./figures`), nimble models (`./nimble`), and the R code (`./R`).

This document describes all of the files present in the repository.

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

## The working directory
Aside from the aforementioned folders, the working directory here stores the `.gitignore` file for this repository, this README file (`README.md`) and the `.Rproj` file (for if you are using RStudio, `skunks.Rproj`).

[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

## The data folder (`./data`)
This folder has ____ files and 1 sub-folder. 

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


**| `./data/site_covariates.csv` |** The spatial covariate (urbanization score, distance to water, and proportion of managed open space) values for each site. This csv file is generated from the `creating_spatial_points.R`, `dist_water_cov.R`, `managed_lawn_cov.R`, and `scaled_covariates.R` scripts located in the `./R` folder. This csv file has 107 rows and 7 columns.

| Column     | Data Type | Explanation                                                                                                                                                                     |
| ---------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Site       | `factor`  | The site associated to this data point                                                                                                                                          |
| HU10       | `numeric` | The mean housing density per meter<sup>2</sup> within a 1 kilometer radius of a site                                                                                            |
| tree       | `numeric` | The mean proportion of tree canopy cover within a 1 kilometer radius of a site                                                                                                  |
| imperv     | `numeric` | The mean proportion of impervious cover within a 1 kilometer radius of a site                                                                                                   |
| urb        | `numeric` | Urban intensity metric created from a principal component analysis using mean housing density, mean tree cover, and mean impervious cover within a 1 kilometer radius of a site |
| water_dist | `numeric` | The shortest Euclidian distance between a site and a permanent body of water in meters                                                                                          |
| open_dev   | `numeric` | The mean proportion of developed open space (areas in which impervious surfaces account for less than 20 percent of total cover) within a 1 kilometer radius of a site          |


**| `./data/original_z_sim_values.csv` |**

**| `./data/scaled_simulation_covariates.csv` |**

**| `./data/phi_gamma_delta_by_season.RDS` |**

**| `./data/simulation_neighbors.RDS` |**

**| `./data/raw_sim_covariates/dist2water.csv` |**

**| `./data/raw_sim_covariates/point_locs.csv` |**

**| `./data/raw_sim_covariates/urb_covars.csv` |**

**| `./data/raw_sim_covariates/urban_openspace.csv` |**


[Back to table of contents ⤒](#links-to-different-parts-of-the-readme-file)

## The R folder (`./R`)

This folder contains 35 R scripts. We have grouped these scripts based upon the function they perform in the overall scheme of the analysis and the relative order in which they would be run if someone was interested in recreating this analysis.

**Group 1 - Functions**  
These scripts contain various functions that are sourced in other scripts to complete the analysis

| File                    | Description                                                                  | Packages Required     |
| ----------------------- | ---------------------------------------------------------------------------- | --------------------- |
| **init_functions.R**    | Functions to assign initial values to parameters during model run            | Only Base R used      |
| **mcmc_functions.R**    | Functions to extract the MCMC posterior from raw model RDS outputs           | Only Base R used      |
| **raster_extraction.R** | Functions to extract spatial data from raster files and calculate proportion | `cli`, `raster`, `sf` |



**Group 2 - Preparing Spatial Covariates**  
These scripts source and extract spatial data, calculate covariate values for each site, and compile data into the `./data/site_covariates.csv`

| File                          | Description                                                                                                                                                                                                                     | Packages Required                           |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| **creating_spatial_points.R** | Creates spatial points for each site and projects them into UTM                                                                                                                                                                 | `dplyr`, `sf`                               |
| **dist_water_cov.R**          | Generates the shortest Euclidian distances between each site and a permanent body of water                                                                                                                                      | `dplyr`, `sf`                               |
| **managed_lawn_cov.R**        | Generates the mean proportion of developed open space (areas in which impervious surfaces account for less than 20 percent of total cover) within a 1 kilometer radius of a site                                                | `dplyr`, `sf`                               |
| **scaled_covariates.R**       | Extracts the mean housing density, proportion of tree canopy cover, and the proportion of developed open space, and generates the urban intensity metric from a principal component analysis using the aforementioned variables | `dplyr`, `raster`, `sf`, `uwinspatialtools` |

**Group 3 - Data Prep**
These scripts prepare the camera trap detection/non-detection data, spatial covariates, and temporal covariates for the specified model
