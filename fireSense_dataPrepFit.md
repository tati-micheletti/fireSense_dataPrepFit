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


|paramName               |paramClass   |default      |min |max |paramDesc                                                                                                                                      |
|:-----------------------|:------------|:------------|:---|:---|:----------------------------------------------------------------------------------------------------------------------------------------------|
|.plotInitialTime        |numeric      |NA           |NA  |NA  |Describes the simulation time at which the first plot event should occur.                                                                      |
|.plotInterval           |numeric      |NA           |NA  |NA  |Describes the simulation time interval between plot events.                                                                                    |
|.saveInitialTime        |numeric      |NA           |NA  |NA  |Describes the simulation time at which the first save event should occur.                                                                      |
|.saveInterval           |numeric      |NA           |NA  |NA  |This describes the simulation time interval between save events.                                                                               |
|.studyAreaName          |character    |             |NA  |NA  |studyArea name that will be appended to file-backed rasters                                                                                    |
|.useCache               |logical      |FALSE        |NA  |NA  |Should this entire module be run with caching activated? This is intended for data-type modules, where stochasticity and time are not relevant |
|areaMultiplier          |numeric,.... |fireSens.... |NA  |NA  |Either a scalar that will buffer areaMultiplier fireSize or a function of fireSize. Default is 2. See fireSenseUtils::bufferToArea for help    |
|climateGCM              |character    |13GCMs_e.... |NA  |NA  |Global Circulation Model to use for climate projections: currently '13GCMs_ensemble' or 'CanESM5'.                                             |
|climateSSP              |numeric      |370          |NA  |NA  |SSP emissions scenario for `climateGCM`: one of 245, 370, or 585.                                                                              |
|cutoffForYoungAge       |numeric      |15           |NA  |NA  |Age at and below which pixels are considered 'young' --> young <- age <= cutoffForYoungAge                                                     |
|fireYears               |integer      |2001, 20.... |NA  |NA  |A numeric vector indicating which years should be extracted from the fire databases to use for fitting                                         |
|forestedLCC             |numeric      |1, 2, 3,.... |NA  |NA  |Forested land cover classes. These classes will be excluded from the PCA.                                                                      |
|igAggFactor             |numeric      |40           |1   |NA  |aggregation factor for rasters during ignition prep                                                                                            |
|ignitionFuelClassCol    |character    |FuelClass    |NA  |NA  |the column in sppEquiv that defines unique fuel classes for ignition                                                                           |
|minBufferSize           |numeric      |5000         |NA  |NA  |Minimum size of buffer and nonbuffer. This is imposed after multiplier on the bufferToArea fn                                                  |
|missingLCCgroup         |character    |nonFores.... |NA  |NA  |if a pixel is forested but is absent from cohortData, it will be grouped in this class. Must be one of the names in `sim$nonForestedLCCGroups` |
|nonflammableLCC         |numeric      |13, 16, .... |NA  |NA  |non-flammable LCC in `sim$rstLCC`.                                                                                                             |
|PCAcomponentsForClimate |numeric      |1            |1   |NA  |number of PCA components to include from climate variables                                                                                     |
|PCAcomponentsForTerrain |numeric      |1            |1   |NA  |currently unused - may be needed if using separate terrain and veg PCAs                                                                        |
|PCAcomponentsForVeg     |numeric      |10           |1   |NA  |number of veg and terrain components to include in GLM                                                                                         |
|PCAcomponentsFromGLM    |numeric      |5            |0   |NA  |the number of components to select from GLM model of `burn ~ PCAcomponents`.                                                                   |
|plotPCA                 |logical      |TRUE         |NA  |NA  |plot the PCA components with a heat map                                                                                                        |
|sppEquivCol             |character    |LandR        |NA  |NA  |column name in sppEquiv object that defines unique species in cohortData                                                                       |
|spreadFuelClassCol      |character    |FuelClass    |NA  |NA  |if using fuel classes for spread, the column in sppEquiv that defines unique fuel classes                                                      |
|useCentroids            |logical      |TRUE         |NA  |NA  |Should fire ignitions start at the `sim$firePolygons` centroids (TRUE) or at the ignition points in `sim$firePoints`?                          |
|usePCA                  |logical      |TRUE         |NA  |NA  |use PCA approach to covariates, as opposed to fuel class approach                                                                              |
|whichModulesToPrepare   |character    |fireSens.... |NA  |NA  |Which fireSense fit modules to prep? defaults to all 3                                                                                         |

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


|objectName               |objectClass              |desc                                                                                                                                     |sourceURL |
|:------------------------|:------------------------|:----------------------------------------------------------------------------------------------------------------------------------------|:---------|
|cohortData2001           |data.table               |Table that defines the cohorts by pixelGroup in 2001                                                                                     |NA        |
|cohortData2011           |data.table               |Table that defines the cohorts by pixelGroup in 2011                                                                                     |NA        |
|spreadFirePoints         |list                     |list of spatialPointsDataFrame for each fire yearwith each point denoting an ignition location                                           |NA        |
|firePolys                |list                     |List of SpatialPolygonsDataFrames representing annual fire polygons.List must be named with followign convention: 'year<numeric year>'   |NA        |
|firePolysForAge          |list                     |firePolys used to classify timeSinceDisturbance in nonforest LCC                                                                         |NA        |
|flammableRTM             |RasterLayer              |RTM without ice/rocks/urban/water. Flammable map with 0 and 1.                                                                           |NA        |
|historicalClimateRasters |list                     |list of historical climate variables in raster stack form, name according to variable                                                    |NA        |
|ignitionFirePoints       |list                     |list of spatialPolygonDataFrame objects representing annual ignition locations. This includes all fires regardless of size               |NA        |
|nonForestedLCCGroups     |list                     |a named list of non-forested landcover groups e.g. list('wetland' = c(19, 23, 32)) These will become covariates in fireSense_IgnitionFit |NA        |
|pixelGroupMap2001        |RasterLayer              |RasterLayer that defines the pixelGroups for cohortData table in 2001                                                                    |NA        |
|pixelGroupMap2011        |RasterLayer              |RasterLayer that defines the pixelGroups for cohortData table in 2011                                                                    |NA        |
|rasterToMatch            |RasterLayer              |template raster for study area. Assumes some buffering of core area to limit edge effect of fire.                                        |NA        |
|rstLCC                   |RasterLayer              |Raster of land cover. Defaults to LCC05.                                                                                                 |NA        |
|sppEquiv                 |data.table               |table of LandR species equivalencies                                                                                                     |NA        |
|standAgeMap2001          |RasterLayer              |map of stand age in 2001 used to create cohortData2001                                                                                   |NA        |
|standAgeMap2011          |RasterLayer              |map of stand age in 2011 used to create cohortData2011                                                                                   |NA        |
|studyArea                |SpatialPolygonsDataFrame |studyArea that determines spatial boundaries of all data                                                                                 |NA        |
|terrainCovariates        |RasterStack              |a raster stack of terrain covariates; defaults are elev, aspect, slope, TRI, TWI                                                         |NA        |

## Output data

Description of the module outputs.


|objectName                             |objectClass |desc                                                                                                                                                                                                                                                                                         |
|:--------------------------------------|:-----------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|climateComponentsToUse                 |character   |names of the climate components or variables needed for FS models                                                                                                                                                                                                                            |
|coefficientPrintOut                    |data.table  |Coefficients from the logit model                                                                                                                                                                                                                                                            |
|componentPrintOut                      |data.table  |A data.table showing the PCA axes and their loadings with the different covariates, e.g., fuel, TPI, HLI, etc.                                                                                                                                                                               |
|fireBufferedListDT                     |list        |list of data.tables with fire id, pixelID, and buffer status                                                                                                                                                                                                                                 |
|firePolys                              |list        |list of spatialPolygonDataFrame objects representing annual fires                                                                                                                                                                                                                            |
|fireSense_annualSpreadFitCovariates    |list        |list of tables with climate PCA components, burn status, polyID, and pixelID                                                                                                                                                                                                                 |
|fireSense_escapeCovariates             |data.table  |ignition covariates with added column of escapes                                                                                                                                                                                                                                             |
|fireSense_escapeFormula                |character   |formula for escape, using fuel classes and landcover, as character                                                                                                                                                                                                                           |
|fireSense_ignitionCovariates           |data.table  |table of aggregated ignition covariates with annual ignitions                                                                                                                                                                                                                                |
|fireSense_ignitionFormula              |character   |formula for ignition, using fuel classes and landcover, as character                                                                                                                                                                                                                         |
|fireSense_nonAnnualSpreadFitCovariates |list        |list of two tables with veg PCA components, burn status, polyID, and pixelID                                                                                                                                                                                                                 |
|fireSense_spreadFormula                |character   |formula for spread, using climate and terrain components, as character                                                                                                                                                                                                                       |
|ignitionFitRTM                         |RasterLayer |A (template) raster with information with regards to the spatial resolution and geographical extent of fireSense_ignitionCovariates. Used to pass this information onto fireSense_ignitionFitted Needs to have number of non-NA cells as attribute (`ignitionFitRTM@data@attributes$nonNAs`) |
|landcoverDT                            |data.table  |data.table with pixelID and relevant landcover classes that is used by predict functions                                                                                                                                                                                                     |
|nonForest_timeSinceDisturbance2001     |RasterLayer |time since burn for non-forested pixels in 2001                                                                                                                                                                                                                                              |
|nonForest_timeSinceDisturbance2011     |RasterLayer |time since burn for non-forested pixels in 2011                                                                                                                                                                                                                                              |
|PCAclimate                             |prcomp      |PCA model for climate covariates, needed for fireSensePredict                                                                                                                                                                                                                                |
|PCAcoeffPlot                           |gglot       |ggplot with PCA loadings for axes used to predict spread                                                                                                                                                                                                                                     |
|PCAveg                                 |prcomp      |PCA model for veg and LCC covariates, needed for FS models                                                                                                                                                                                                                                   |
|spreadFirePoints                       |list        |list of spatialPolygonDataFrame objects representing annual fire centroids. This only includes fires that escaped (e.g. size > flammableRTM resolution                                                                                                                                       |
|terrainDT                              |data.table  |data.table with pixelID and relevant terrain variables used by predict models                                                                                                                                                                                                                |
|vegComponentsToUse                     |character   |names of the veg components to use in ignition, escape, and spread predict models                                                                                                                                                                                                            |

# Links to other modules

Outputs used by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit`.
