# A repository for:

Kase, A et al. The impact of local and long-distance colonization in a fragmented landscape

## Links to different parts of the readme file

1. [What is in this repository?](##what-is-in-this-repository?)
2. [The working directory](##the-working-directory)
3. [The data folder (`./data`)](##the-data-folder-(`./data`))

## What is in this repository?
This repository stores all of the data and code used to fit a dynamic occupancy model with explicit colonization terms, compare model outputs using Brier scores, and forecast model estimates across space and time. The folder organization separates the data (`./data`), figures from the manuscript (`./figures`), nimble models (`./nimble`), and the R code (`./R`).

This document describes all of the files present in the repository.

[Back to table of contents ⤒](##links-to-different-parts-of-the-readme-file)

## The working directory
Aside from the aforementioned folders, the working directory here stores the `.gitignore` file for this repository, this README file (`README.md`) and the `.Rproj` file (for if you are using RStudio, `skunks.Rproj`).

[Back to table of contents ⤒](##links-to-different-parts-of-the-readme-file)

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

