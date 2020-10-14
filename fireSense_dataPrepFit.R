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
  reqdPkgs = list('raster', 'sf', 'sp', 'data.table', 'PredictiveEcology/fireSenseUtils (>=0.0.2)', 'parallel'),
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
    defineParameter('areaMultiplier', class = 'numeric', 1, NA, NA,
                    desc = paste('Either a scalar that will buffer areaMultiplier * fireSize or a function',
                    'of fireSize. Default is 1. See fireSenseUtils::bufferToArea for help')),
    defineParameter(name = "fireYears", class = "integer", default = 1991:2017,
                    desc = "A numeric vector indicating which years should be extracted
                    from the fire databases to use for fitting"),
    defineParameter(name = 'minBufferSize', 'numeric', 500, NA, NA,
                    desc = "Minimum size of buffer and nonbuffer. This is imposed after multiplier on the bufferToArea fn"),
    defineParameter(name = "useCentroids", class = "logical", default = TRUE,
                    desc = paste("Should fire ignitions start at the sim$firePolygons centroids (TRUE)",
                                 "or at the ignition points in sim$firePoints?")),
    defineParameter("whichModulesToPrepare", "character",
                    default = c("fireSense_IgnitionFit", "fireSense_IgnitionPredict", "fireSense_EscapeFit"),
                    NA, NA, desc = "Which fireSense fit modules to prep? defaults to all 3")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = 'DEM', objectClass = 'RasterLayer',
                 sourceURL = 'https://drive.google.com/file/d/121x_CfWy2XP_-1av0cYE7sxUfb4pmsup/view?usp=sharing',
                 desc = "DEM for deriving terrain metrics"),
    expectsInput(objectName = "cohortData2001", objectClass = "data.table", sourceURL = NA,
                 desc = paste0("Table that defines the cohorts by pixelGroup in 2001")),
    expectsInput(objectName = "cohortData2011", objectClass = "data.table", sourceURL = NA,
                 desc = paste0("Table that defines the cohorts by pixelGroup in 2011")),
    expectsInput(objectName = 'DEM', objectClass = 'RasterLayer',
                 sourceURL = 'https://drive.google.com/file/d/121x_CfWy2XP_-1av0cYE7sxUfb4pmsup/view?usp=sharing',
                 desc = "DEM for deriving terrain metrics"),
    expectsInput(objectName = 'firePoints', objectClass = 'list', sourceURL = NA,
                 desc = paste0("list of spatialPointsDataFrame for each fire year",
                               "with each point denoting an ignition location")),
    expectsInput(objectName = "firePolys", objectClass = "list", sourceURL = NA,
                 desc = paste0("List of SpatialPolygonsDataFrames representing annual fire polygons.",
                               "This defaults to https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/ and uses ",
                               "the most current versions of the database (Nov or Sept 2019).",
                               "List must be named with followign convention: 'year<numeric year>'")),
    expectsInput(objectName = "flammableRTM", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
    expectsInput(objectName = "pixelGroupMap2001", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RasterLayer that defines the pixelGroups for cohortData table in 2001"),
    expectsInput(objectName = "pixelGroupMap2011", objectClass = "RasterLayer",
                 desc = "RasterLayer that defines the pixelGroups for cohortData table in 2011"),
    expectsInput(objectName = "rstLCC", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "Raster of land cover. Defaults to LCC05."),
    expectsInput(objectName = 'historicalClimateRasters', objectClass = 'RasterStack', sourceURL = NA,
                 desc = 'historical climate variables corresponding to "fireYears" parameter'),
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
    createsOutput(objectName = 'fireSense_terrainCovariates', objectClass = 'rasterStack',
                  desc = 'terrain covariates used by both fitting and predicting'),
    createsOutput(objectName = 'fireSense_escapeFitCovariates', objectClass = 'data.table', desc = 'WIP'),
    createsOutput(objectName = 'fireSense_ignitionFitCovariates', objectClass = 'data.table', desc = 'WIP'),
    createsOutput(objectName = 'fireSense_spreadFitCovariates', objectClass = 'data.table',
                  desc = 'WIP - a data.table with a bunch of covariates for burned and non-burned classes')
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

  # some data prep that is essential to all three modules will happen here - e.g. reclassifying LCC?
  #every fitting will need the cohortData objects split into cohortDataLong, with a wide layout with biomass and landcover
  #Each fitting event will prepare the data tables of relevant events (ignition, spread, or escape) and MDC
  #then, each event should join the tables
  if (is.null(sim$firePolys[[1]]$FIRE_ID)) {
    stop("firePolys needs a numeric FIRE_ID column")
  }

  if (!is.numeric(sim$firePolys[[1]]$FIRE_ID)) {
    message("need numeric FIRE_ID column in fire polygons. Coercing to numeric...")
    #this is an annoying trait in the current NFBB
    origNames <- names(sim$firePolys)
    PointsAndPolys <- lapply(names(sim$firePolys), FUN = function(year, polys = sim$firePolys, points = sim$firePoints) {
      polys <- polys[[year]]
      points <- points[[year]] #get matching year
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
                              rasterToMatch = sim$rasterToMatch,
                              verb = TRUE, areaMultiplier = P(sim)$areaMultiplier,
                              field = "FIRE_ID",
                              cores = length(sim$firePolys),
                              minSize = P(sim)$minBufferSize,
                              userTags = c("bufferToArea"),
                              omitArgs = "cores")

  ###################################################
  # Post buffering, new issues --> must make sure points and buffers match
  ###################################################
  sim$firePoints <- Cache(harmonizeBufferAndPoints, cent = sim$firePoints,
                          buff = fireBufferedListDT,
                          ras = sim$rasterToMatch,
                          idCol = "FIRE_ID",
                          userTags = c("harmonizeBufferAndPoints"))

  # get pixelIDs pre2005 and post2005
  browser()
  pre2005 <- names(fireBufferedListDT) %in% paste0('year', min(P(sim)$fireYears):2005)

  pre2005Indices <- fireBufferedListDT[pre2005]
  pre2005Climate <- unstack(sim$historicalClimateRasters)[pre2005] # names are wrong because prepInputs isn't being used to source them :(
  post2005Indices <- fireBufferedListDT[!pre2005]
  post2005Climate <- unstack(sim$historicalClimateRasters)[!pre2005]


  cohorts2001 <- castCohortData(cohortData = sim$cohortData2001,
                                pixelGroupMap = sim$pixelGroupMap2001,
                                climateRasters = pre2005Climate,
                                terrainRasters = sim$terrainCovariates,
                                index = pre2005Indices)

  cohorts2011 <- fireSenseUtils::castCohortData(cohortData = sim$cohortData2011)




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

### template for your event1
Event1 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event2
Event2 <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test

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

  if (!suppliedElsewhere('landcoverGroups', sim)) {
    sim$landcoverGroups <- list(
      'forest' = c(1:15),
      'wetland' = c(19, 31, 32),
      'nonTree' = c(16:18, 20:30),
      'nonVeg' = c(33, 36:38) #wait what are 34/35
    )
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


  #this is important because centroids may be fewer than fires if fire polys were small
  min1Fire <- lapply(sim$firePoints, length) > 0
  sim$firePoints <- sim$firePoints[min1Fire]
  sim$firePolys <- sim$firePolys[min1Fire]
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
                                   landcover = sim$rstLCC,
                                   userTags = c('terrainCovariates')
    )

  }


  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
