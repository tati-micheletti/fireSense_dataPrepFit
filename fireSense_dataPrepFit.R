## Everything in this file gets sourced during `simInit()`,
## and all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used without `sim$` as they are namespaced to the module,
## just like functions in R packages.
## If exact location is required, functions will be: `sim$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "fireSense_dataPrepFit",
  description = "",
  keywords = "",
  authors = structure(list(list(given = c("Ian"), family = "Eddy", role = c("aut", "cre"), email = "ian.eddy@canada.ca", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.4.9003", fireSense_dataPrepFit = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "fireSense_dataPrepFit.Rmd")),
  reqdPkgs = list('raster', 'sf', 'sp', 'data.table', 'PredictiveEcology/fireSenseUtils (>=0.0.2)', 'parallel', 'fastDummies'),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated? This is intended",
                          "for data-type modules, where stochasticity and time are not relevant")),
    defineParameter('areaMultiplier', class = 'numeric', 2, NA, NA,
                    desc = paste('Either a scalar that will buffer areaMultiplier * fireSize or a function',
                    'of fireSize. Default is 2. See fireSenseUtils::bufferToArea for help')),
    defineParameter(name = "fireYears", class = "integer", default = 1991:2017,
                    desc = "A numeric vector indicating which years should be extracted
                    from the fire databases to use for fitting"),
    defineParameter(name = 'forestedLCC', class = 'numeric', default = c(1:15, 20, 34, 35), NA, NA,
                    desc = paste0('forested land cover classes. If using LCC 2005, this should also include burn classes 34 and 35.',
                                  'These classes will be excluded from the PCA')),
    defineParameter(name = 'nonflammableLCC', class = 'numeric', c(0, 33, 36, 37, 38, 39), NA, NA,
                    desc = 'non-flammable LCC in sim$rstLCC'),
    defineParameter(name = 'minBufferSize', 'numeric', 500, NA, NA,
                    desc = "Minimum size of buffer and nonbuffer. This is imposed after multiplier on the bufferToArea fn"),
    defineParameter(name = 'PCAcomponentsForClimate', 'numeric', 1, 1, NA,
                    desc = 'number of PCA components to include from climate variables'),
    defineParameter(name = 'PCAcomponentsForTerrain', 'numeric', 1, 1, NA,
                    desc = 'currently unused - may be needed if using separate terrain and veg PCAs'),
    defineParameter(name = 'PCAcomponentsForVeg', 'numeric', 10, 1, NA,
                    desc = 'number of veg and terrain components to include in GLM'),
    defineParameter(name = 'PCAcomponentsFromGLM', 'numeric', 4, 0, NA,
                    desc = 'the number of components to select from GLM model of burn ~ PCAcomponents' ),
    defineParameter(name = "useCentroids", class = "logical", default = TRUE,
                    desc = paste("Should fire ignitions start at the sim$firePolygons centroids (TRUE)",
                                 "or at the ignition points in sim$firePoints?")),
    defineParameter("whichModulesToPrepare", "character",
                    default = c("fireSense_IgnitionFit", "fireSense_IgnitionPredict", "fireSense_EscapeFit"),
                    NA, NA, desc = "Which fireSense fit modules to prep? defaults to all 3")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "cohortData2001", objectClass = "data.table", sourceURL = NA,
                 desc = paste0("Table that defines the cohorts by pixelGroup in 2001")),
    expectsInput(objectName = "cohortData2011", objectClass = "data.table", sourceURL = NA,
                 desc = paste0("Table that defines the cohorts by pixelGroup in 2011")),
    expectsInput(objectName = 'firePoints', objectClass = 'list', sourceURL = NA,
                 desc = paste0("list of spatialPointsDataFrame for each fire year",
                               "with each point denoting an ignition location")),
    expectsInput(objectName = "firePolys", objectClass = "list", sourceURL = NA,
                 desc = paste0("List of SpatialPolygonsDataFrames representing annual fire polygons.",
                               "List must be named with followign convention: 'year<numeric year>'")),
    expectsInput(objectName = "flammableRTM", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
    expectsInput(objectName = 'nonForestedLCCGroups', objectClass = 'list',
                 desc = paste("a named list of non-forested landcover groups",
                              "e.g. list('wetland' = c(19, 23, 32))")),
    expectsInput(objectName = "pixelGroupMap2001", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RasterLayer that defines the pixelGroups for cohortData table in 2001"),
    expectsInput(objectName = "pixelGroupMap2011", objectClass = "RasterLayer",
                 desc = "RasterLayer that defines the pixelGroups for cohortData table in 2011"),
    expectsInput(objectName = "rstLCC", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "Raster of land cover. Defaults to LCC05."),
    expectsInput(objectName = 'historicalClimateRasters', objectClass = 'list', sourceURL = NA,
                 desc = 'list of historical climate variables in raster stack form, name according to variable'),
    expectsInput(objectName = 'rasterToMatch', objectClass = 'RasterLayer', sourceURL = NA,
                 desc = "template raster for study area"),
    expectsInput(objectName = 'studyArea', objectClass = 'SpatialPolygonsDataFrame', sourceURL = NA,
                 desc = "studyArea that determines spatial boundaries of all data"),
    expectsInput(objectName = 'terrainCovariates', objectClass = 'RasterStack', sourceURL = NA,
                 desc = 'a raster stack of terrain covariates; defaults are elev, aspect, slope, TRI, TWI')
  ),
  outputObjects = bindrows(
    createsOutput(objectName = 'firePolys', objectClass = 'list',
                  desc = 'list of spatialPolygonDataFrame objects representing annual fires'),
    createsOutput(objectName = 'firePoints', objectClass = 'list',
                  desc = 'list of spatialPolygonDataFrame objects representing annual fire centroids'),
    createsOutput(objectName = 'fireSense_annualFitCovariates', objectClass = 'data.table',
                  desc = 'table of climate PCA components, burn status, polyID, and pixelID'),
    createsOutput(objectName = 'fireSense_nonAnnualFitCovariates', objectClass = 'data.table',
                  desc = 'table of veg PCA components, burn status, polyID, and pixelID'),
    createsOutput(objectName = 'fireSense_spreadLogitModel', objectClass = 'glm',
                  desc = 'GLM with burn as dependent variable and PCA components as covariates'),
    createsOutput(objectName = 'landcoverDT', 'data.table',
                  desc = paste('data.table with pixelID and relevant landcover classes',
                               'that is used by predict functions')),
    createsOutput(objectName = 'PCAclimate', objectClass = 'prcomp',
                  desc = 'PCA model for climate covariates, needed for fireSensePredict'),
    createsOutput(objectName = 'PCAveg', objectClass = 'prcomp',
                  desc = 'PCA model for veg and LCC covariates, needed for fireSensePredict'),
    createsOutput(objectName = 'vegComponentsToUse', 'character',
                  desc = 'names of the veg components to use in ignition, escape, and spread predict models')
  )
))

## event types
#   - type `init` is required for initialization

doEvent.fireSense_dataPrepFit = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      if ("fireSense_IgnitionFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepIgnitionFitData")
      if ("fireSense_EscapeFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepEscapeFitData")
      if ("fireSense_SpreadFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepSpreadFitData")

    },
    prepIgnitionFitData = {
      #fill as needed
    },
    prepEscapeFitData = {
      #fill as needed
    },
    prepSpreadFitData = {
      #this will handle the specific needs of spread not done in Init already

    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  ####prep terrain for PCA####
  terrainDT <- lapply(names(sim$terrainCovariates), FUN = function(x){
    y <- data.table(getValues(sim$terrainCovariates[[x]]))
  }) %>%
    dplyr::bind_cols(.)
  setnames(terrainDT, names(sim$terrainCovariates))
  setDT(terrainDT) #needed to pre-allocate space for new columns

  set(terrainDT, j = 'pixelID', value = 1:ncell(sim$pixelGroupMap2001))
  set(terrainDT, j = 'flammable', value = getValues(sim$flammableRTM))
  terrainDT <- terrainDT[flammable == 1] %>%
    set(., NULL, 'flammable', NULL) %>%
    na.omit(.)
  sim$terrainDT <- terrainDT #unclear if terrain will be separate PCA or not.
  #if separate, just keep the PCA object, terrain won't change. If not, keep DT for predictions
  #need to make this more flexible so terrain can be combined or not
  rm(terrainDT)

  ####prep landcover for PCA####
  #this can be made into a fireSenseUtils function, but it is only done once - here.
  lcc <- data.table(pixelID = 1:ncell(sim$rstLCC),
                    lcc = getValues(sim$rstLCC),
                    flammable = getValues(sim$flammableRTM)) %>%
    .[!lcc == 0,] %>%
    .[!is.na(flammable),] %>%
    .[flammable == 1,] %>%
    .[, lcc := as.factor(lcc)]
  set(lcc, NULL, 'flammable', NULL)
  lcc <- Cache(fastDummies::dummy_cols, .data = lcc,
               select_columns = 'lcc',
               remove_selected_columns = TRUE,
               remove_first_dummy = FALSE, #no need if all forest LC is removed anyway
               ignore_na = TRUE,
               userTags = c("dummy_cols"))

  forestedLCC <- paste0('lcc_', P(sim)$forestedLCC)
  forestedLCC <- colnames(lcc)[colnames(lcc) %in% forestedLCC]
  set(lcc, NULL, j = forestedLCC, NULL) #removes columns of forested LCC

  #Group dummy columns into similar landcovers
  lccColsPreGroup <- colnames(lcc)[!colnames(lcc) %in% c('pixelID')]
  setDT(lcc) #pre-allocate space for new columns

  for (i in names(sim$nonForestedLCCGroups)) {
    class <- paste0('lcc_', sim$nonForestedLCCGroups[[i]])
    set(lcc, NULL,  eval(i), rowSums(lcc[, .SD, .SDcols = class]))
  }
  rm(class)
  set(lcc, NULL, lccColsPreGroup, NULL)
  sim$landcoverDT <- lcc
  #save the lcc - it will be used by predict models and there's no reason to re-do these steps

  #do veg PCA
  cohorts2001 <- castCohortData(cohortData = sim$cohortData2001,
                                pixelGroupMap = sim$pixelGroupMap2001,
                                lcc = lcc,
                                terrainDT = sim$terrainDT)
  set(cohorts2001, NULL, 'year', 2001)

  cohorts2011 <- castCohortData(cohortData = sim$cohortData2011,
                                pixelGroupMap = sim$pixelGroupMap2011,
                                lcc = lcc,
                                terrainDT = sim$terrainDT)
  set(cohorts2011, NULL, 'year', 2011)

  vegPCAdat <- rbind(cohorts2001, cohorts2011)
  rm(cohorts2001, cohorts2011, lcc)

  ###*** move this chunk to fireSenseUtils - it should do PCA and convert to int ***###
  vegTerrainPCA <- prcomp(vegPCAdat[, .SD, .SDcols = !c("pixelGroup", 'pixelID', 'year')],
                          center = TRUE, scale. = TRUE, rank = 10)
  sim$PCAveg <- vegTerrainPCA

  #store as Integer
  vegComponents <- as.data.table(vegTerrainPCA$x * 1000)
  vegComponents <- vegComponents[, lapply(.SD, asInteger), .SDcols = colnames(vegComponents)]
  vegComponents[, pixelID := vegPCAdat$pixelID]
  ###*** end of chunk ***###

  #year is preserved in fit, but not predict
  vegComponents[, year := vegPCAdat$year] #need to preserve for logistic

  rm(vegPCAdat)

  components <- paste0('PC', 1:P(sim)$PCAcomponentsForVeg)
  vegComponents <- vegComponents[, .SD, .SDcols = c(components, 'pixelID', 'year')]
  #rename components so climate/veg components distinguishable
  setnames(vegComponents, old = components, new = paste0("veg", components))
  rm(components)

  ####prep Climate components####
  flammableIndex <- data.table(index = 1:ncell(sim$flammableRTM), value = getValues(sim$flammableRTM)) %>%
    .[value == 1,] %>%
    .$index

  #this involves a large object (27 years of raster data converted to 3 columns), year, pixelID, value
  #need it in one dt for PCA, but it is less efficient so we convert back to list of annual dts
  climatePCAdat <- Cache(climateRasterToDataTable,
                         historicalClimateRasters = sim$historicalClimateRasters,
                         Index = flammableIndex,
                         userTags = c('climateRasterToDataTable'))
  rm(flammableIndex)

  if (length(climatePCAdat) > 1) {
    #havent planned this out yet, do call-merge  isnt working - we want a wide layout of n variables
    warning("running fireSense_dataPrepFit with two climate components is still in development")
    climatePCAdat$by <- c("pixelID", "year")
    #this merge fails - we want to go from lists of list of data.tables to a list of data.tables, with one variable for each stack
    #TODO: need to solve this with dummy data
    climatePCAdat <- do.call(what = merge.data.table, args = climatePCAdat)
    climatePCA <- prcomp(climatePCAdat[, .SD, .SDcols = !c("pixelID", "year")],
                         center = TRUE,
                         scale. = TRUE)
    sim$climatePCA <- climatePCA
    climateComponents <- as.data.table(climatePCA$x * 1000)
    climateComponents <- climateComponents[, lapply(.SD, asInteger), .SDcols = colnames(climateComponents)]
    set(climateComponents, NULL,'pixelID', climatePCAdat$pixelID)
    set(climateComponents, NULL,'year', climatePCAdat$year)
    components <- paste0('PC', 1:P(sim)$PCAcomponentsForClimate)
    climateComponents <- climateComponents[, .SD, .SDcols = c(components, 'pixelID', 'year')]
    setnames(climateComponents, old = components, new = paste0("climate", components))
    rm(components, climatePCA)
    } else {
      #if there is only one climate variable, no PCA
      climateComponents <-  climatePCAdat[[1]]
    }

    #put back into list form to reduce object size
    years <- sort(unique(climateComponents$year))
    fireSense_annualFitCovariates <- lapply(years, FUN = function(x){
      thisYear <- climateComponents[year == x,]
      thisYear <- thisYear[, .SD, .SDcols = -c("year")]
      return(thisYear)
    })
    names(fireSense_annualFitCovariates) <- years


  ####prep fire data####
  if (is.null(sim$firePolys[[1]]$FIRE_ID)) {
    stop("firePolys needs a numeric FIRE_ID column")
  }

  if (!is.numeric(sim$firePolys[[1]]$FIRE_ID)) {
    message("need numeric FIRE_ID column in fire polygons. Coercing to numeric...")
    #this is true of the current NFBB
    origNames <- names(sim$firePolys)
    PointsAndPolys <- lapply(names(sim$firePolys), FUN = function(year, polys = sim$firePolys, points = sim$firePoints) {
      polys <- polys[[year]]
      points <- points[[year]]
      #ensure matching IDs
      points <- points[points$FIRE_ID %in% polys$FIRE_ID,]
      polys <- polys[polys$FIRE_ID %in% points$FIRE_ID,]
      points$FIRE_ID <- as.numeric(as.factor(points$FIRE_ID))
      polys$FIRE_ID <- as.numeric(as.factor(polys$FIRE_ID))
      return(list(polys = polys, points = points))
    })
    sim$firePoints <- lapply(PointsAndPolys, FUN = function(x) {return(x[['points']])})
    sim$firePolys <- lapply(PointsAndPolys, FUN = function(x) {return(x[['polys']])})
    rm(PointsAndPolys)
    names(sim$firePolys) <- origNames
    names(sim$firePoints) <- origNames
  }

  fireBufferedListDT <- Cache(bufferToArea,
                              poly = sim$firePolys,
                              polyName = names(sim$firePolys),
                              rasterToMatch = sim$flammableRTM,
                              verb = TRUE, areaMultiplier = P(sim)$areaMultiplier,
                              field = "FIRE_ID",
                              cores = length(sim$firePolys),
                              minSize = P(sim)$minBufferSize,
                              userTags = c("bufferToArea"),
                              omitArgs = "cores")

  # Post buffering, new issues --> must make sure points and buffers match
  sim$firePoints <- Cache(harmonizeBufferAndPoints, cent = sim$firePoints,
                          buff = fireBufferedListDT,
                          ras = sim$flammableRTM,
                          idCol = "FIRE_ID",
                          userTags = c("harmonizeBufferAndPoints"))

  # get pixelIDs pre2005 and post2005
  #these will become lists of stacks
  pre2005 <- paste0('year', min(P(sim)$fireYears):2005)
  pre2005Indices <- fireBufferedListDT[names(fireBufferedListDT) %in% pre2005]
  post2005Indices <- fireBufferedListDT[!names(fireBufferedListDT) %in% pre2005]

  #takes veg PCA, optionally terrain PCA, and optional index to build covariates
  logisticCovariatesPre2005 <- rbindlist(pre2005Indices) %>%
    vegComponents[year < 2005][., on = c("pixelID")]
  logisticCovariatesPost2005 <- rbindlist(post2005Indices) %>%
    vegComponents[year > 2005][., on = c("pixelID")]

  # we want to collapse both time steps for logistic regression
  #but need to preserve structure of named lists for spreadFit
  fireSenseVegData <- rbind(logisticCovariatesPost2005, logisticCovariatesPre2005)
  # there are some NAs due to postProcess  - burned cells on border of terrain covariates
  setnames(fireSenseVegData, 'buffer', 'burned')

  ####run logistic regression####
  #build logistic formula
  logitFormula <- grep('veg*', names(fireSenseVegData), value = TRUE) %>%
    paste0(., collapse = ' +') %>%
    paste0('burned ~ ', .)

  fireSenseLogit <- glm(formula = eval(logitFormula),
                        data = fireSenseVegData,
                        family = 'binomial',
                        na.action = na.exclude)

  #take largest coeffiecients as they are mean-centered and scaled, number determined by param
  bestComponents <- sort(abs(fireSenseLogit$coefficients[2:length(fireSenseLogit$coefficients)]),
                         decreasing = TRUE)[1:P(sim)$PCAcomponentsFromGLM]
  sim$vegComponentsToUse <- names(bestComponents) #the values will be wrong due to abs, just take names
  sim$fireSense_spreadLogitModel <- fireSenseLogit

  colsToExtract <- c('year', 'pixelID', 'burned', 'ids', sim$vegComponentsToUse)

  #save two data.tables, one with climate, one with veg
  sim$fireSense_nonAnnualFitCovariates <- fireSenseVegData[, .SD, .SDcols = colsToExtract]
  sim$fireSense_annualFitCovariates <- fireSense_annualFitCovariates

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  #Plot(sim$object)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("studyArea", sim)) {
    stop("Please supply study area - this object is key")
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- prepInputsLCC(destinationPath = dPath,
                                       studyArea = sim$studyArea,
                                       useCache = TRUE)
  }

  if (!all(suppliedElsewhere('cohortData2011', sim),
           suppliedElsewhere("pixelGroupMap2011", sim),
           suppliedElsewhere("pixelGroupMap2001", sim),
           suppliedElsewhere('cohortData2001', sim))) {
    stop("Stop - need cohortData and pixelGroupMap objects - contact module creators")
  }

  if (!suppliedElsewhere("firePolys", sim)){

    sim$firePolys <- Cache(fireSenseUtils::getFirePolygons, years = P(sim)$fireYears,
                           studyArea = sim$studyArea,
                           destinationPath = dPath,
                           useInnerCache = TRUE,
                           userTags = c('firePolys', paste0("years:", range(P(sim)$fireYears))))
  }
  if (isTRUE(P(sim)$useCentroids)) {
    if (!suppliedElsewhere("firePoints", sim)){
      message("... preparing polyCentroids")

      centerFun <- function(x){
        if (is.null(x)) {
          return(NULL)
        } else {
          ras <- x
          ras$ID <- 1:NROW(ras)
          centCoords <- rgeos::gCentroid(ras, byid = TRUE)
          cent <- SpatialPointsDataFrame(centCoords,
                                         as.data.frame(ras))
          return(cent)
        }
      }

      sim$firePoints <- Cache(FUN = parallel::mclapply,
                              X = sim$firePolys,
                              mc.cores = pemisc::optimalClusterNum(2e3, maxNumClusters = length(sim$firePolys)),
                              centerFun, #don't specify FUN argument or Cache will mistake it.
                              userTags = c(currentModule(sim), 'firePoints'),
                              omitArgs = c("userTags", "mc.cores", "useCloud", "cloudFolderID"))
      names(sim$firePoints) <- names(sim$firePolys)

    }
  } else {

    if (!suppliedElsewhere("firePoints", sim)){
      sim$firePoints <- Cache(getFirePoints_NFDB,
                              url = "http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point.zip",
                              studyArea = sim$studyArea,
                              rasterToMatch = sim$rasterToMatch,
                              NFDB_pointPath = file.path(Paths$inputPath, "NFDB_point"),
                              years = P(sim)$fireYears,
                              userTags = c("what:firePoints", "forWhat:fireSense_SpreadFit"))
      crs(sim$firePoints) <- crs(sim$rasterToMatch)
      names(sim$firePoints) <- names(sim$firePolys)
    }
  }

  if (all(!is.null(sim$firePoints), !is.null(sim$firePolys))) { #may be NULL if passed by objects - add to Init?
  #this is necessary because centroids may be fewer than fires if fire polys were small
    min1Fire <- lapply(sim$firePoints, length) > 0
    sim$firePoints <- sim$firePoints[min1Fire]
    sim$firePolys <- sim$firePolys[min1Fire]
  }
  if (length(sim$firePolys) != length(sim$firePoints)) {
    stop("mismatched years between firePolys and firePoints")
    #need to implement a better approach that matches each year's IDS
  }

  if (!suppliedElsewhere("rstLCC", sim)){
    sim$rstLCC <- prepInputsLCC(destinationPath = dPath,
                                studyArea = sim$studyArea,
                                useCache = TRUE)
  }

  if (!suppliedElsewhere("historicalClimateRasters", sim)) {
    stop("please supply sim$historicalClimateRasters")
  }

  if (!suppliedElsewhere('terrainCovariates', sim)) {
    sim$terrainCovariates <- Cache(fireSenseUtils::prepTerrainCovariates,
                                   studyArea = sim$studyArea,
                                   rasterToMatch = sim$rasterToMatch,
                                   dPath = dPath,
                                   userTags = c('terrainCovariates'))
  }

  if (!suppliedElsewhere("flammableRTM", sim)) {
    sim$flammableRTM <- LandR::defineFlammable(sim$rstLCC,
                                               nonFlammClasses = P(sim)$nonflammableLCC,
                                               mask = sim$rasterToMatch,
                                               filename2 = file.path(dPath, 'flammableRTM.tif'))
  }

  if (!suppliedElsewhere('nonForestedLCCGroups', sim)) {
    sim$nonForestedLCCGroups <- list(
      'shrubland' = c(16, 22),
      'wetland' = c(19, 23, 32),
      'cropland' = c(26, 27, 28, 29),
      'grassland' = c(17, 18, 21, 24, 26, 30)
    )
  }

  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
