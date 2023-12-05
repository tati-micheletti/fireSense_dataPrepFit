---
title: "fireSense_dataPrepFit"
author: 
  - Ian Eddy
  - Alex Chubaty
date: "September 2020; updated February 2022"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---



# Overview

Prepare data required by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit`.

# Parameters

Provide a summary of user-visible parameters.


|paramName                  |paramClass   |default      |min |max |paramDesc                                                                                                                                                                                                                                                                                                                                       |
|:--------------------------|:------------|:------------|:---|:---|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|areaMultiplier             |numeric,.... |fireSens.... |NA  |NA  |Either a scalar that will buffer `areaMultiplier fireSize` or a function of `fireSize`. See `?fireSenseUtils::bufferToArea`.                                                                                                                                                                                                                    |
|bufferForFireRaster        |numeric      |1000         |0   |NA  |The distance that determine whether separate patches of burned pixels originated from the same fire. Only relevant when 'useRasterizedFireForSpread' is TRUE. This param is separate from 'minBufferSize', which is used to determine the minimum sample of burned and unburned pixels to include in each fire.                                 |
|cutoffForYoungAge          |numeric      |15           |NA  |NA  |Age at and below which pixels are considered 'young' (`young <- age <= cutoffForYoungAge`)                                                                                                                                                                                                                                                      |
|fireYears                  |integer      |2001, 20.... |NA  |NA  |A numeric vector indicating which years should be extracted from the fire databases to use for fitting                                                                                                                                                                                                                                          |
|forestedLCC                |numeric      |1, 2, 3,.... |NA  |NA  |Forested land cover classes. These classes will be excluded from the PCA.                                                                                                                                                                                                                                                                       |
|igAggFactor                |numeric      |40           |1   |NA  |aggregation factor for rasters during ignition prep.                                                                                                                                                                                                                                                                                            |
|ignitionFuelClassCol       |character    |FuelClass    |NA  |NA  |the column in `sppEquiv` that defines unique fuel classes for ignition                                                                                                                                                                                                                                                                          |
|minBufferSize              |numeric      |5000         |NA  |NA  |Minimum number of cells in buffer and nonbuffer. This is imposed after the multiplier on the `bufferToArea` fn                                                                                                                                                                                                                                  |
|missingLCCgroup            |character    |nonFores.... |NA  |NA  |if a pixel is forested but is absent from `cohortData`, it will be grouped in this class. Must be one of the names in `sim$nonForestedLCCGroups`                                                                                                                                                                                                |
|nonflammableLCC            |numeric      |13, 16, .... |NA  |NA  |non-flammable LCC in `sim$rstLCC`.                                                                                                                                                                                                                                                                                                              |
|nonForestCanBeYoungAge     |logical      |TRUE         |NA  |NA  |if TRUE, burned non-forest will be treated as `youngAge`                                                                                                                                                                                                                                                                                        |
|sppEquivCol                |character    |LandR        |NA  |NA  |column name in `sppEquiv` object that defines unique species in `cohortData`                                                                                                                                                                                                                                                                    |
|spreadFuelClassCol         |character    |FuelClass    |NA  |NA  |if using fuel classes for spread, the column in `sppEquiv` that defines unique fuel classes                                                                                                                                                                                                                                                     |
|useCentroids               |logical      |TRUE         |NA  |NA  |Should fire ignitions start at the `sim$firePolygons` centroids (TRUE) or at the ignition points in `sim$firePoints`?                                                                                                                                                                                                                           |
|useRasterizedFireForSpread |logical      |FALSE        |NA  |NA  |Should rasterized fire be used in place of a vectorized fire dataset? This method attributes burned pixels to specific fires, only examines the latest fire in a pixel, and may be subject to temporal error. is therefore more appropriate in areas with low rates of fire, or where the NFDB dataset may be incomplete (ie northern Ontario). |
|whichModulesToPrepare      |character    |fireSens.... |NA  |NA  |Which fireSense fit modules to prep? defaults to all 3                                                                                                                                                                                                                                                                                          |
|.plotInitialTime           |numeric      |NA           |NA  |NA  |Describes the simulation time at which the first plot event should occur.                                                                                                                                                                                                                                                                       |
|.plotInterval              |numeric      |NA           |NA  |NA  |Describes the simulation time interval between plot events.                                                                                                                                                                                                                                                                                     |
|.saveInitialTime           |numeric      |NA           |NA  |NA  |Describes the simulation time at which the first save event should occur.                                                                                                                                                                                                                                                                       |
|.saveInterval              |numeric      |NA           |NA  |NA  |This describes the simulation time interval between save events.                                                                                                                                                                                                                                                                                |
|.studyAreaName             |character    |             |NA  |NA  |`studyArea` name that will be appended to file-backed rasters                                                                                                                                                                                                                                                                                   |
|.useCache                  |logical      |FALSE        |NA  |NA  |Should this entire module be run with caching activated? This is intended for data-type modules, where stochasticity and time are not relevant                                                                                                                                                                                                  |

# Events

Describe what happens for each event type.

## Plotting

Write what is plotted.

## Saving

Write what is saved.

# Data dependencies

## Input data

How to obtain input data, and a description of the data required by the module.
If `sourceURL` is specified, `downloadData("fireSense_dataPrepFit", "..")` may be sufficient.


|objectName               |objectClass |desc                                                                                                                                             |sourceURL                                                                      |
|:------------------------|:-----------|:------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------|
|cohortData2001           |data.table  |Table that defines the cohorts by pixelGroup in 2001                                                                                             |NA                                                                             |
|cohortData2011           |data.table  |Table that defines the cohorts by pixelGroup in 2011                                                                                             |NA                                                                             |
|spreadFirePoints         |list        |named list of spatial points for each fire year with each point denoting an ignition location.                                                   |NA                                                                             |
|firePolys                |list        |List of sf polygon objects representing annual fire polygons.List must be named with followign convention: 'year<numeric year>'                  |NA                                                                             |
|firePolysForAge          |list        |list of fire polygons used to classify `timeSinceDisturbance` in nonforest LCC                                                                   |NA                                                                             |
|historicalFireRaster     |SpatRaster  |a raster with values representing fire year 1985-2020                                                                                            |https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip |
|flammableRTM             |SpatRaster  |RTM without ice/rocks/urban/water. Flammable map with 0 and 1.                                                                                   |NA                                                                             |
|historicalClimateRasters |list        |length-one list of containing a raster stack of historical climate list named after the variable and raster layers named as 'year<numeric year>' |NA                                                                             |
|ignitionFirePoints       |list        |list of sf polygon objects representing annual ignition locations. This includes all fires regardless of size                                    |NA                                                                             |
|nonForestedLCCGroups     |list        |a named list of non-forested landcover groups e.g. list('wetland' = c(19, 23, 32)) These will become covariates in `fireSense_IgnitionFit`       |NA                                                                             |
|pixelGroupMap2001        |SpatRaster  |SpatRaster that defines the `pixelGroups` for cohortData table in 2001                                                                           |NA                                                                             |
|pixelGroupMap2011        |SpatRaster  |SpatRaster that defines the `pixelGroups` for cohortData table in 2011                                                                           |NA                                                                             |
|rasterToMatch            |SpatRaster  |template raster for study area. Assumes some buffering of core area to limit edge effect of fire.                                                |NA                                                                             |
|rstLCC                   |SpatRaster  |Raster of land cover. Defaults to LCC05.                                                                                                         |NA                                                                             |
|sppEquiv                 |data.table  |table of LandR species equivalencies                                                                                                             |NA                                                                             |
|standAgeMap2001          |SpatRaster  |map of stand age in 2001 used to create `cohortData2001`                                                                                         |NA                                                                             |
|standAgeMap2011          |SpatRaster  |map of stand age in 2011 used to create `cohortData2011`                                                                                         |NA                                                                             |
|studyArea                |sf          |studyArea that determines spatial boundaries of all data                                                                                         |NA                                                                             |

## Output data

Description of the module outputs.


|objectName                             |objectClass |desc                                                                                                                                                                                                                                                                                          |
|:--------------------------------------|:-----------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|fireBufferedListDT                     |list        |list of data.tables with fire id, `pixelID`, and buffer status                                                                                                                                                                                                                                |
|firePolys                              |list        |list of sf polygon objects representing annual fires                                                                                                                                                                                                                                          |
|fireSense_annualSpreadFitCovariates    |list        |list of tables with climate covariates, `youngAge`, burn status, `polyID`, and `pixelID`                                                                                                                                                                                                      |
|fireSense_escapeCovariates             |data.table  |ignition covariates with added column of escapes                                                                                                                                                                                                                                              |
|fireSense_escapeFormula                |character   |formula for escape, using fuel classes and landcover, as character                                                                                                                                                                                                                            |
|fireSense_ignitionCovariates           |data.table  |table of aggregated ignition covariates with annual ignitions                                                                                                                                                                                                                                 |
|fireSense_ignitionFormula              |character   |formula for ignition, using climate and vegetation covariates, as character                                                                                                                                                                                                                   |
|fireSense_nonAnnualSpreadFitCovariates |list        |list of two tables with vegetation covariates, burn status, polyID, and `pixelID`                                                                                                                                                                                                             |
|fireSense_spreadFormula                |character   |formula for spread, using climate and vegetation covariates, as character                                                                                                                                                                                                                     |
|ignitionFitRTM                         |SpatRaster  |A (template) raster with information with regards to the spatial resolution and geographical extent of `fireSense_ignitionCovariates`. Used to pass this information onto `fireSense_ignitionFitted` Needs to have number of non-NA cells as attribute (`attributes(ignitionFitRTM)$nonNAs`). |
|landcoverDT                            |data.table  |data.table with `pixelID` and relevant landcover classes that is used by predict functions.                                                                                                                                                                                                   |
|nonForest_timeSinceDisturbance2001     |SpatRaster  |time since burn for non-forested pixels in 2001                                                                                                                                                                                                                                               |
|nonForest_timeSinceDisturbance2011     |SpatRaster  |time since burn for non-forested pixels in 2011                                                                                                                                                                                                                                               |
|spreadFirePoints                       |list        |Named list of `sf` polygon objects representing annual fire centroids. This only includes fires that escaped (e.g. `size > res(flammableRTM)`.                                                                                                                                                |

# Links to other modules

Outputs used by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit`.
