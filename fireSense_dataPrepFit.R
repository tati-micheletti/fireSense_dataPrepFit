defineModule(sim, list(
  name = "fireSense_dataPrepFit",
  description = "",
  keywords = "",
  authors = structure(list(list(given = c("Ian"), family = "Eddy", role = c("aut", "cre"),
                                email = "ian.eddy@canada.ca", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.4.9003", fireSense_dataPrepFit = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "fireSense_dataPrepFit.Rmd")),
  reqdPkgs = list("raster", "sf", "sp", "data.table", "PredictiveEcology/fireSenseUtils (>=0.0.3.0000)",
                  "parallel", "fastDummies", "spatialEco", "snow"),
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
    defineParameter(".studyAreaName", "character", NULL, NA, NA,
                    desc = "studyArea name that will be appended to file-backed rasters"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated? This is intended",
                          "for data-type modules, where stochasticity and time are not relevant")),
    defineParameter("areaMultiplier", class = "numeric", 2, NA, NA,
                    desc = paste("Either a scalar that will buffer areaMultiplier * fireSize or a function",
                    "of fireSize. Default is 2. See fireSenseUtils::bufferToArea for help")),
    defineParameter(name = "fireYears", class = "integer", default = 2001:2019,
                    desc = "A numeric vector indicating which years should be extracted
                    from the fire databases to use for fitting"),
    defineParameter(name = "forestedLCC", class = "numeric", default = c(1:15, 20, 32, 34, 35), NA, NA,
                    desc = paste0("forested land cover classes. If using LCC 2005, this should also include burn classes 34 and 35.",
                                  "These classes will be excluded from the PCA")),
    defineParameter(name = "minBufferSize", "numeric", 500, NA, NA,
                    desc = "Minimum size of buffer and nonbuffer. This is imposed after multiplier on the bufferToArea fn"),
    defineParameter(name = 'missingLCCgroup', class = 'character', 'nonForest_highFlam', NA, NA,
                    desc = paste("if a pixel is forested but is absent from cohortData, it will be grouped in this class.",
                                 "Must be one of the names in sim$nonForestedLCCGroups")),
    defineParameter(name = "nonflammableLCC", class = "numeric", c(0, 25, 30, 33, 36, 37, 38, 39), NA, NA,
                    desc = "non-flammable LCC in sim$rstLCC"),
    defineParameter(name = "PCAcomponentsForClimate", "numeric", 1, 1, NA,
                    desc = "number of PCA components to include from climate variables"),
    defineParameter(name = "PCAcomponentsForTerrain", "numeric", 1, 1, NA,
                    desc = "currently unused - may be needed if using separate terrain and veg PCAs"),
    defineParameter(name = "PCAcomponentsForVeg", "numeric", 10, 1, NA,
                    desc = "number of veg and terrain components to include in GLM"),
    defineParameter(name = "PCAcomponentsFromGLM", "numeric", 5, 0, NA,
                    desc = "the number of components to select from GLM model of burn ~ PCAcomponents" ),
    defineParameter(name = "sppEquivCol", class = "character", default = "LandR", NA, NA,
                    desc = "column name in sppEquiv object that defines unique species in cohortData"),
    defineParameter(name = "useCentroids", class = "logical", default = TRUE,
                    desc = paste("Should fire ignitions start at the sim$firePolygons centroids (TRUE)",
                                 "or at the ignition points in sim$firePoints?")),
    defineParameter(name = "whichModulesToPrepare", class = "character",
                    default = c("fireSense_IgnitionFit", "fireSense_IgnitionPredict", "fireSense_EscapeFit"),
                    NA, NA, desc = "Which fireSense fit modules to prep? defaults to all 3")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "cohortData2001", objectClass = "data.table", sourceURL = NA,
                 desc = paste0("Table that defines the cohorts by pixelGroup in 2001")),
    expectsInput(objectName = "cohortData2011", objectClass = "data.table", sourceURL = NA,
                 desc = paste0("Table that defines the cohorts by pixelGroup in 2011")),
    expectsInput(objectName = "spreadFirePoints", objectClass = "list", sourceURL = NA,
                 desc = paste0("list of spatialPointsDataFrame for each fire year",
                               "with each point denoting an ignition location")),
    expectsInput(objectName = "firePolys", objectClass = "list", sourceURL = NA,
                 desc = paste0("List of SpatialPolygonsDataFrames representing annual fire polygons.",
                               "List must be named with followign convention: 'year<numeric year>'")),
    expectsInput(objectName = 'firePolysForAge', objectClass = 'list', sourceURL = NA,
                 desc=  'firePolys used to classify timeSinceDisturbance in nonforest LCC'),
    expectsInput(objectName = "flammableRTM", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
    expectsInput(objectName = "historicalClimateRasters", objectClass = "list", sourceURL = NA,
                 desc = "list of historical climate variables in raster stack form, name according to variable"),
    expectsInput(objectName = "ignitionFirePoints", objectClass = "list", sourceURL = NA,
                  desc = paste("list of spatialPolygonDataFrame objects representing annual ignition locations.",
                               "This includes all fires regardless of size")),
    expectsInput(objectName = "nonForestedLCCGroups", objectClass = "list",
                 desc = paste("a named list of non-forested landcover groups",
                              "e.g. list('wetland' = c(19, 23, 32))",
                              "These will become covariates in fireSense_IgnitionFit")),
    expectsInput(objectName = "pixelGroupMap2001", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "RasterLayer that defines the pixelGroups for cohortData table in 2001"),
    expectsInput(objectName = "pixelGroupMap2011", objectClass = "RasterLayer",
                 desc = "RasterLayer that defines the pixelGroups for cohortData table in 2011"),
    expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "template raster for study area"),
    expectsInput(objectName = 'rasterToMatch', objectClass = 'RasterLayer', sourceURL = NA,
                 desc = "template raster for study area"),
    expectsInput(objectName = "rstLCC", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "Raster of land cover. Defaults to LCC05."),
    expectsInput(objectName = "sppEquiv", objectClass = "data.table", sourceURL = NA,
                 desc = "table of LandR species equivalencies"),
    expectsInput(objectName = 'standAgeMap2001', objectClass = "RasterLayer", sourceURL = NA,
                 desc = "map of stand age in 2001 used to create cohortData2001"),
    expectsInput(objectName = 'standAgeMap2011', objectClass = "RasterLayer", sourceURL = NA,
                 desc = "map of stand age in 2011 used to create cohortData2011"),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", sourceURL = NA,
                 desc = "studyArea that determines spatial boundaries of all data"),
    expectsInput(objectName = "terrainCovariates", objectClass = "RasterStack", sourceURL = NA,
                 desc = "a raster stack of terrain covariates; defaults are elev, aspect, slope, TRI, TWI")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "climateComponentsToUse", objectClass = "character",
                  desc = "names of the climate components or variables needed for FS models"),
    createsOutput(objectName = "fireBufferedListDT", objectClass = "list",
                  desc = "list of data.tables with fire id, pixelID, and buffer status"),
    createsOutput(objectName = "firePolys", objectClass = "list",
                  desc = "list of spatialPolygonDataFrame objects representing annual fires"),
    createsOutput(objectName = "fireSense_annualSpreadFitCovariates", objectClass = "list",
                  desc = "list of tables with climate PCA components, burn status, polyID, and pixelID"),
    createsOutput(objectName = "fireSense_ignitionCovariates", objectClass = "data.table",
                  desc = "table of aggregated ignition covariates with annual ignitions"),
    createsOutput(objectName = "fireSense_ignitionFormula", objectClass = "character",
                  desc = "formula for ignition, using fuel classes and landcover, as character"),
    createsOutput(objectName = "fireSense_nonAnnualSpreadFitCovariates", objectClass = "list",
                  desc = "list of two tables with veg PCA components, burn status, polyID, and pixelID"),
    createsOutput(objectName = "fireSense_spreadFormula", objectClass = "character",
                  desc = "formula for spread, using climate and terrain components, as character"),
    createsOutput(objectName = "fireSense_spreadLogitModel", objectClass = "glm",
                  desc = "GLM with burn as dependent variable and PCA components as covariates"),
    createsOutput(objectName = "landcoverDT", "data.table",
                  desc = paste("data.table with pixelID and relevant landcover classes",
                               "that is used by predict functions")),
    createsOutput(objectName = "nonForest_timeSinceDisturbance", objectClass = "RasterLayer",
                  desc = "time since burn for non-forested pixels"),
    createsOutput(objectName = "PCAclimate", objectClass = "prcomp",
                  desc = "PCA model for climate covariates, needed for fireSensePredict"),
    createsOutput(objectName = "PCAveg", objectClass = "prcomp",
                  desc = "PCA model for veg and LCC covariates, needed for FS models"),
    createsOutput(objectName = "spreadFirePoints", objectClass = "list",
                  desc = paste("list of spatialPolygonDataFrame objects representing annual fire centroids.",
                               "This only includes fires that escaped (e.g. size > flammableRTM resolution")),
    createsOutput(objectName = "terrainDT", "data.table",
                  desc = "data.table with pixelID and relevant terrain variables used by predict models"),
    createsOutput(objectName = "vegComponentsToUse", "character",
                  desc = "names of the veg components to use in ignition, escape, and spread predict models")
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

      if (!all(P(sim)$whichModulesToPrepare %in% c("fireSense_SpreadFit",
                                                   "fireSense_IgnitionFit",
                                                   "fireSense_EscapeFit"))){
        stop("unrecognized module to prepare - review parameter whichModulesToPrepare")
        #the camelcase is still different with FS from LandR Biomass
      }
      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      if ("fireSense_IgnitionFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepIgnitionFitData")
      if ("fireSense_EscapeFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepEscapeFitData")
      if ("fireSense_SpreadFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepSpreadFitData")

      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "plotAndMessage", eventPriority = 9)
      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "cleanUp", eventPriority = 10) #cleans up Mod objects

    },
    prepIgnitionFitData = {
      sim <- prepare_IgnitionFit(sim)
    },
    prepEscapeFitData = {
      #fill as needed
    },
    prepSpreadFitData = {
      sim <- prepare_SpreadFit(sim)
    },

    plotAndMessage = {
      sim <- plotAndMessage(sim)
    },

    cleanUp = {
      sim <- cleanUpMod(sim)
    },
    warning(paste("Undefined event type: \"", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {

  ####prep fire data ####
  if (is.null(sim$firePolys[[1]]$FIRE_ID)) {
    stop("firePolys needs a numeric FIRE_ID column")
  }

  if (!is.numeric(sim$firePolys[[1]]$FIRE_ID)) {
    message("need numeric FIRE_ID column in fire polygons. Coercing to numeric...")
    #this is true of the current NFBB
    origNames <- names(sim$firePolys)
    PointsAndPolys <- lapply(names(sim$firePolys), FUN = function(year, polys = sim$firePolys, points = sim$spreadFirePoints) {
      polys <- polys[[year]]
      points <- points[[year]]
      #ensure matching IDs
      points <- points[points$FIRE_ID %in% polys$FIRE_ID,]
      polys <- polys[polys$FIRE_ID %in% points$FIRE_ID,]
      points$FIRE_ID <- as.numeric(as.factor(points$FIRE_ID))
      polys$FIRE_ID <- as.numeric(as.factor(polys$FIRE_ID))
      return(list(polys = polys, points = points))
    })
    sim$spreadFirePoints <- lapply(PointsAndPolys, FUN = function(x) {return(x[["points"]])})
    sim$firePolys <- lapply(PointsAndPolys, FUN = function(x) {return(x[["polys"]])})
    rm(PointsAndPolys)
    names(sim$firePolys) <- origNames
    names(sim$spreadFirePoints) <- origNames
  }

  nCores <- ifelse(grep('*Windows', osVersion), 1, length(sim$firePolys))
  fireBufferedListDT <- Cache(bufferToArea,
                              poly = sim$firePolys,
                              polyName = names(sim$firePolys),
                              rasterToMatch = sim$flammableRTM,
                              verb = TRUE, areaMultiplier = P(sim)$areaMultiplier,
                              field = "FIRE_ID",
                              cores = nCores,
                              minSize = P(sim)$minBufferSize,
                              userTags = c("bufferToArea"),
                              omitArgs = "cores")

  # Post buffering, new issues --> must make sure points and buffers match
  sim$spreadFirePoints <- Cache(harmonizeBufferAndPoints, cent = sim$spreadFirePoints,
                                buff = fireBufferedListDT,
                                ras = sim$flammableRTM,
                                idCol = "FIRE_ID", #this is different from default
                                userTags = c("harmonizeBufferAndPoints"))

  ####prep terrain for PCA####
  terrainDT <- lapply(names(sim$terrainCovariates), FUN = function(x){
    y <- data.table(getValues(sim$terrainCovariates[[x]]))
  }) %>%
    dplyr::bind_cols(.)
  setnames(terrainDT, names(sim$terrainCovariates))
  setDT(terrainDT) #needed to pre-allocate space for new columns

  set(terrainDT, j = "pixelID", value = 1:ncell(sim$pixelGroupMap2001))
  set(terrainDT, j = "flammable", value = getValues(sim$flammableRTM))
  terrainDT <- terrainDT[flammable == 1,] %>%
    set(., NULL, "flammable", NULL) %>%
    na.omit(.)
  sim$terrainDT <- terrainDT #unclear if terrain will be separate PCA or not.
  #if separate, just keep the PCA object, terrain won't change.
  #it should be combined with lcc for a single table that can be referenced during predict
  rm(terrainDT)

  ## prep landcover for PCA
  # this can be made into a fireSenseUtils function, but it is only done once - here.
  lcc <- data.table(pixelID = 1:ncell(sim$rstLCC),
                    lcc = getValues(sim$rstLCC),
                    flammable = getValues(sim$flammableRTM)) %>%
    .[lcc != 0,] %>%
    .[!is.na(flammable),] %>%
    .[flammable == 1,] %>%
    .[, lcc := as.factor(lcc)]
  set(lcc, NULL, "flammable", NULL)
  lcc <- Cache(fastDummies::dummy_cols, .data = lcc,
               select_columns = "lcc",
               remove_selected_columns = TRUE,
               remove_first_dummy = FALSE, #no need if all forest LC is removed anyway
               ignore_na = TRUE,
               userTags = c("dummy_cols"))

  forestedLCC <- paste0("lcc_", P(sim)$forestedLCC)
  forestedLCC <- colnames(lcc)[colnames(lcc) %in% forestedLCC]
  set(lcc, NULL, j = forestedLCC, NULL) #removes columns of forested LCC

  #Group dummy columns into similar landcovers
  lccColsPreGroup <- colnames(lcc)[!colnames(lcc) %in% c("pixelID")]
  setDT(lcc) #pre-allocate space for new columns

  for (i in names(sim$nonForestedLCCGroups)) {
    classes <- paste0("lcc_", sim$nonForestedLCCGroups[[i]])
    classes <- classes[classes %in% colnames(lcc)]
    set(lcc, NULL,  eval(i), rowSums(lcc[, .SD, .SDcols = classes]))
  }
  rm(classes)
  set(lcc, NULL, lccColsPreGroup, NULL)
  sim$landcoverDT <- lcc
  #save the lcc - it will be used by predict models

# cannot merge because before subsetting due to column differences over time

  mod$firePolysForAge <- lapply(sim$firePolysForAge[lengths(sim$firePolysForAge) > 0], FUN = sf::st_as_sf) %>%
    lapply(., FUN = function(x){
        x <- x[, "YEAR"]
      }) %>%
    do.call(rbind, .)

  cohorts2001 <-  castCohortData(cohortData = sim$cohortData2001,
                                 pixelGroupMap = sim$pixelGroupMap2001,
                                 year = 2001,
                                 ageMap = NULL,
                                 lcc = sim$landcoverDT,
                                 terrainDT = sim$terrainDT,
                                 missingLCC = P(sim)$missingLCC)

  cohorts2011 <-  castCohortData(cohortData = sim$cohortData2011,
                                 pixelGroupMap = sim$pixelGroupMap2011,
                                 year = 2011,
                                 ageMap = NULL,
                                 lcc = sim$landcoverDT,
                                 terrainDT = sim$terrainDT,
                                 missingLCC = P(sim)$missingLCC)

  vegPCAdat <- rbind(cohorts2001, cohorts2011)

  rm(cohorts2001, cohorts2011, lcc)

  #Clean up missing pixels - this is a temporary fix
  #we will always have NAs because of edge pixels - will be an issue when predicting
  origNames <- names(fireBufferedListDT)
  fireBufferedListDT <- lapply(fireBufferedListDT, FUN = function(year){
    missingPixels <- year[!pixelID %in% vegPCAdat$pixelID]
    if (nrow(missingPixels) > 0) {
      year <- year[!pixelID %in% missingPixels$pixelID]
    }
    return(year)
  })
  names(fireBufferedListDT) <- origNames
  rm(origNames)

  ###*predict will run castCohortData and then makeVegTerrainPCA for predicting
  vegList <- makeVegTerrainPCA(dataForPCA = vegPCAdat)
  vegComponents <- vegList$vegComponents
  sim$PCAveg <- vegList$vegTerrainPCA
  rm(vegList, vegPCAdat)

  components <- paste0("PC", 1:P(sim)$PCAcomponentsForVeg)
  vegComponents <- vegComponents[, .SD, .SDcols = c(components, "pixelID", "year", "youngAge")]
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
                         userTags = c("climateRasterToDataTable"))
  rm(flammableIndex)

  if (length(climatePCAdat) > 1) {
    warning("running fireSense_dataPrepFit with two climate components is still in development")
    ## TODO: this is untested
    climatePCAdat <- Reduce(x = climatePCAdat, function(x, y, ...) merge(x, y , ...))
    climatePCA <- prcomp(climatePCAdat[, .SD, .SDcols = !c("pixelID", "year")],
                         center = TRUE,
                         scale. = TRUE)
    sim$PCAclimate <- climatePCA
    climateComponents <- as.data.table(climatePCA$x * 1000)
    climateComponents <- climateComponents[, lapply(.SD, asInteger), .SDcols = colnames(climateComponents)]
    set(climateComponents, NULL,"pixelID", climatePCAdat$pixelID)
    set(climateComponents, NULL,"year", climatePCAdat$year)
    components <- paste0("PC", 1:P(sim)$PCAcomponentsForClimate)
    climateComponents <- climateComponents[, .SD, .SDcols = c(components, "pixelID", "year")]
    setnames(climateComponents, old = components, new = paste0("climate", components))
    rm(components, climatePCA)
  } else {
    #if there is only one climate variable, no PCA
    climateComponents <- climatePCAdat[[1]]
    setDT(climateComponents)
    climCol <- names(climateComponents)[!colnames(climateComponents) %in% c("pixelID", "year")]
    #scale climCol to have unit variance and mean center
    set(climateComponents, NULL, "climPCA1",
        scale(climateComponents[, .SD, .SDcols = climCol], center = TRUE, scale = TRUE))
    climateComponents[, climPCA1 := asInteger(climPCA1 * 1000)]
    climateComponents <- climateComponents[, .SD, .SDcols = c("pixelID", "climPCA1", "year")]
  }
  #this is to construct the formula,
  #whether there are multiple climate components or a single non-transformed variable
  sim$climateComponentsToUse <- names(climateComponents)[!names(climateComponents) %in% c("pixelID", "year")]

  # get pixelIDs pre2011 and post2011
  #these will become lists of stacks
  pre2011 <- paste0("year", min(P(sim)$fireYears):2010)
  pre2011Indices <- fireBufferedListDT[names(fireBufferedListDT) %in% pre2011]
  post2011Indices <- fireBufferedListDT[!names(fireBufferedListDT) %in% pre2011]

  logisticCovariatesPre2011 <- rbindlist(pre2011Indices) %>%
    vegComponents[year < 2011][., on = c("pixelID")]
  logisticCovariatesPre2011[is.na(year), year := 2001] #these are non-flammable indices

  logisticCovariatesPost2011 <- rbindlist(post2011Indices) %>%
    vegComponents[year >= 2011][., on = c("pixelID")]
  logisticCovariatesPost2011[is.na(year), year := 2011]

  #Some pixels will be NA because the polygon includes non-flammable cells
  #As long as these pixels are also NA in climate data, no issue

  # we want to collapse both time steps for logistic regression
  #but need to preserve structure of named lists for spreadFit
  fireSenseVegData <- rbind(logisticCovariatesPost2011, logisticCovariatesPre2011)
  # there are some NAs due to postProcess  - burned cells on border of terrain covariates
  setnames(fireSenseVegData, "buffer", "burned")

  ####run logistic regression####
  #build logistic formula
  logitFormula <- grep("veg*", names(fireSenseVegData), value = TRUE) %>%
    paste0(., collapse = " +") %>%
    paste0("burned ~ ", .)

  fireSenseLogit <- glm(formula = eval(logitFormula),
                        data = fireSenseVegData,
                        family = "binomial",
                        na.action = na.exclude)

  #take largest coeffiecients as they are mean-centered and scaled, number determined by param
  bestComponents <- sort(abs(fireSenseLogit$coefficients[2:length(fireSenseLogit$coefficients)]),
                         decreasing = TRUE)[1:P(sim)$PCAcomponentsFromGLM]
  sim$vegComponentsToUse <- names(bestComponents) #the values will be wrong due to abs, just take names
  sim$fireSense_spreadLogitModel <- fireSenseLogit

  RHS <- paste0("youngAge +",
                paste(sim$climateComponentsToUse, collapse = " + "),
                " + ",
                paste(sim$vegComponentsToUse, collapse = " + "))

  mod$fireSenseVegData <- fireSenseVegData
  mod$climateComponents <- climateComponents #needed by prep spread
  #sim$fireSense_spreadFormula <- as.formula(paste("~0 + ", paste(RHS)))
  sim$fireSense_spreadFormula <- paste("~0 + ", paste(RHS))
  sim$fireBufferedListDT <- fireBufferedListDT #needed by DEOptim

  return(invisible(sim))
}

prepare_SpreadFit <- function(sim) {

  #Put in format for DEOptim
  #Prepare annual spread fit covariates
  #this is a funny way to get years but avoids years with 0 fires
  years <- paste0("year", P(sim)$fireYears)
  yearsWithFire <- years[paste0("year", P(sim)$fireYears) %in% names(sim$firePolys)]
  pre2011 <- yearsWithFire[yearsWithFire %in% paste0("year", min(P(sim)$fireYears):2010)]
  post2011 <- yearsWithFire[yearsWithFire %in% paste0("year", 2011:max(P(sim)$fireYears))]

  #currently there are NAs in climate due to non-flammable pixels in fire buffer
  fireSense_annualSpreadFitCovariates <- lapply(years, FUN = function(x) {
    thisYear <- mod$climateComponents[year == x,]
    thisYear <- thisYear[, .SD, .SDcols = -c("year")]
    annualFireSenseDT <- sim$fireBufferedListDT[[x]]
    thisYear <- thisYear[pixelID %in% annualFireSenseDT$pixelID,]
    thisYear <- na.omit(thisYear)
    thisYear <- as.data.table(thisYear)
    #NAs are possible from raster projection issues and non-flammable buffer pixels
    return(thisYear)
  })
  names(fireSense_annualSpreadFitCovariates) <- years


  #prepare non-annual spread fit covariates
  pre2011Indices <- sim$fireBufferedListDT[names(sim$fireBufferedListDT) %in% pre2011]
  post2011Indices <- sim$fireBufferedListDT[!names(sim$fireBufferedListDT) %in% pre2011]
  colsToExtract <- c("pixelID", "youngAge", sim$vegComponentsToUse)

  #remove age from the annual covariates - it will come from TSD
  nonAnnualPre2011 <- mod$fireSenseVegData[year < 2011, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[, youngAge := NULL]
  nonAnnualPost2011 <- mod$fireSenseVegData[year >= 2011, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[, youngAge := NULL]


  #Create one universal TSD map for each initial time period combining stand age/ time since burn
  TSD2001 <- makeTSD(year = 2001, firePolys = sim$firePolysForAge,
                     standAgeMap = sim$standAgeMap2001,
                     lcc = sim$landcoverDT)
  TSD2011 <- makeTSD(year = 2011, firePolys = sim$firePolysForAge,
                     standAgeMap = sim$standAgeMap2011,
                     lcc = sim$landcoverDT)
  #the function will do this below, and then use the data.table with location of non-forest to fill in those ages
  annualCovariates <- Map(f = calcYoungAge,
                          years = list(c(2001:2010), c(2011:max(P(sim)$fireYears))),
                          annualCovariates = list(fireSense_annualSpreadFitCovariates[pre2011],
                                                  fireSense_annualSpreadFitCovariates[post2011]),
                          standAgeMap = list(TSD2001, TSD2011),
                          MoreArgs = list(fireBufferedListDT = sim$fireBufferedListDT)
  )
  sim$fireSense_annualSpreadFitCovariates <- do.call(c, annualCovariates)

  sim$fireSense_nonAnnualSpreadFitCovariates <- list(nonAnnualPre2011, nonAnnualPost2011)
  names(sim$fireSense_nonAnnualSpreadFitCovariates) <- c(paste(names(pre2011Indices), collapse = "_"),
                                                         paste(names(post2011Indices), collapse = "_"))

  return(invisible(sim))
}

prepare_IgnitionFit <- function(sim) {

  #first put landcover into raster stack - it will be aggregated
  putBackIntoRaster <- function(lcc, landcoverDT, templateRas) {
    lccRas <- raster(templateRas)
    lccRas[landcoverDT$pixelID] <- landcoverDT[, get(lcc)]
    return(lccRas)
  }

  LCCs <- Cache(lapply,
                names(sim$nonForestedLCCGroups),
                putBackIntoRaster,
                landcoverDT = sim$landcoverDT,
                templateRas = sim$rasterToMatch,
                userTags = c("putBackIntoRaster"))  %>%
    stack(.) %>%
    raster::aggregate(., fact = 25, fun = mean)
  names(LCCs) <- names(sim$nonForestedLCCGroups)

  fuelClasses <- Map(f = classifyCohortsFireSenseSpread,
                     cohortData = list(sim$cohortData2001, sim$cohortData2011),
                     yearCohort = list(2001, 2011),
                     pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                     MoreArgs = list(flammableMap = sim$flammableRTM,
                                     sppEquiv = sim$sppEquiv,
                                     sppEquivCol = P(sim)$sppEquivCol)) %>%
    lapply(., FUN = raster::brick) %>%
    lapply(., FUN = raster::aggregate, fact = 25, fun = mean)
  names(fuelClasses) <- c("year2001", "year2011")

  #Next aggregate fire data? for each year? Then get climate values. mix and match. done.
  climate <- Cache(lapply,
                   sim$historicalClimateRasters,
                   raster::aggregate,
                   fact = 25,
                   fun = mean,
                   userTags = c("climate", "aggregate"))
  #now we want a table with the climate values of each firePoint at each year
  if (length(climate) > 1) {
    stop("need to fix ignition for multiple climate variables. contact module developers")
    #for now - will fix later
  } else {
    climate <- raster::stack(climate[[1]])
  }
   #ignition won't have same years as spread so we do not use names of init objects
   #The reason is some years may have no significant fires, e.g. 2001 in RIA
   pre2011 <- paste0("year", min(P(sim)$fireYears):2010)
   post2011 <- paste0("year", 2011:max(P(sim)$fireYears))

   sim$fireSense_ignitionCovariates <- Map(f = stackAndExtract,
                     years = list(pre2011, post2011),
                     fuel = list(fuelClasses$year2001, fuelClasses$year2011),
                     MoreArgs = list(LCC = LCCs,
                                      fires = sim$ignitionFirePoints,
                                      climate = climate,
                                      climVar = names(sim$historicalClimateRasters))) %>%
   rbindlist(.) %>%
     as.data.frame(.) #TODO: confirm this
   #Formula naming won't work with >1 climate variable, regardless a stop is upstream

   # recall that fuelClass 1 is actually youngAge.
   RHS <- names(sim$fireSense_ignitionCovariates) %>%
     .[!. %in% c("cells", "nFires", names(sim$historicalClimateRasters))] %>%
     paste0(., ":", names(sim$historicalClimateRasters)) %>%
     paste(., collapse = " + ")
   #review formula
   sim$fireSense_ignitionFormula <- paste("nFires ~0 + ", paste(RHS))

  return(sim)
}

cleanUpMod <- function(sim){
  #because prep for ignition/escape/spread isn't necessarily run, clean up here
  mod$climateComponents <- NULL # remove for memory sake
  mod$firePolysForAge <- NULL
  mod$fireSenseVegData <- NULL

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
plotAndMessage <- function(sim) {

  components <- as.data.table(sim$PCAveg$rotation)
  setnames(components, old = colnames(components), new = paste0("veg", colnames(components)))
  components[, covariate := row.names(sim$PCAveg$rotation)]
  setcolorder(components, neworder = 'covariate')
  componentPrintOut <- components[, .SD, .SDcol = c('covariate', sim$vegComponentsToUse)]

  message("the loading of the components most correlated with fire are:")
  print(componentPrintOut)

  coefficientPrintOut <- sim$fireSense_spreadLogitModel$coefficients[sim$vegComponentsToUse]
  message("the coefficients of these components are:")
  print(coefficientPrintOut)

  return(invisible(sim))
}

.inputObjects <- function(sim) {
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

  if (!all(suppliedElsewhere("cohortData2011", sim),
           suppliedElsewhere("pixelGroupMap2011", sim),
           suppliedElsewhere("pixelGroupMap2001", sim),
           suppliedElsewhere("cohortData2001", sim))) {
    stop("Stop - need cohortData and pixelGroupMap objects - contact module creators")
  }

  if (!suppliedElsewhere("firePolys", sim)) {
    if (suppliedElsewhere("firePolysForAgeMap", sim)) { #don't want to needlessly postProcess the same firePolys objects
      #Maybe this should all be moved to init - Then we source years from 15 years prior to P(sim)$fireYears
      sim$firePolys <- Cache(fireSenseUtils::getFirePolygons, years = P(sim)$fireYears,
                             studyArea = sim$studyArea,
                             destinationPath = dPath,
                             useInnerCache = TRUE,
                             userTags = c("firePolys", paste0("years:", range(P(sim)$fireYears))))
    } else {

      allFirePolys <- Cache(fireSenseUtils::getFirePolygons,
                             years = c(min(P(sim)$fireYears - 15):max(P(sim)$fireYears)), #get enough data for the before years
                             studyArea = sim$studyArea,
                             destinationPath = dPath,
                             useInnerCache = TRUE,
                             userTags = c("firePolys", paste0("years:", range(P(sim)$fireYears))))

      sim$firePolys <- allFirePolys[names(allFirePolys) %in% paste0("year", P(sim)$fireYears)]
      sim$firePolysForAge <- allFirePolys
    }
  }

  if (!suppliedElsewhere("firePolys", sim)) {
    if (is.null(sim$firePolysForAge)) { #this only happens if firePolys supplied but not firePolysForAge
      sim$allFirePolys <- Cache(fireSenseUtils::getFirePolygons,
                            years = c(min(P(sim)$fireYears) - 15):max(P(sim)$fireYears), #get enough data for the before years
                            studyArea = sim$studyArea,
                            destinationPath = dPath,
                            useInnerCache = TRUE,
                            userTags = c("firePolys", paste0("years:", range(P(sim)$fireYears))))
    }
  }

  if (!suppliedElsewhere("standAgeMap2001", sim)) {
    sim$standAgeMap2001 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch,
                                 studyArea = sim$studyArea,
                                 destinationPath = dPath,
                                 filename2 = 'standAgeMap2001.tif',
                                 startTime = 2001,
                                 userTags = c(cacheTags, 'prepInputsStandAgeMap2001'))
  }

  if (!suppliedElsewhere("standAgeMap2011", sim)) {
    sim$standAgeMap2011 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch,
                                 studyArea = sim$studyArea,
                                 destinationPath = dPath,
                                 filename2 = 'standAgeMap2011.tif',
                                 startTime = 2011,
                                 userTags = c(cacheTags, 'prepInputsStandAgeMap2011'))
  }

  if (!suppliedElsewhere("spreadFirePoints", sim)) {
    message("... preparing polyCentroids; starting up parallel R threads")
    centerFun <- function(x) {
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
    # suppress any startup messages
    mc <- pemisc::optimalClusterNum(2e3, maxNumClusters = length(sim$firePolys))
    clObj <- parallel::makeCluster(type = "SOCK", mc)
    a <- parallel::clusterEvalQ(cl = clObj, {library(raster); library(rgeos)})
    clusterExport(cl = clObj, list("firePolys"), envir = sim)
    sim$spreadFirePoints <- Cache(FUN = parallel::clusterApply,
                                  x = sim$firePolys,
                                  cl = clObj,
                                  fun = centerFun, #don't specify FUN argument or Cache will mistake it.
                                  userTags = c(currentModule(sim), "spreadFirePoints"),
                                  omitArgs = c("userTags", "mc.cores", "useCloud", "cloudFolderID"))
    stopCluster(clObj)
    names(sim$spreadFirePoints) <- names(sim$firePolys)
  }

  if (all(!is.null(sim$spreadFirePoints), !is.null(sim$firePolys))) {
  #may be NULL if passed by objects - add to Init?
  #this is necessary because centroids may be fewer than fires if fire polys were small
    min1Fire <- lapply(sim$spreadFirePoints, length) > 0
    sim$spreadFirePoints <- sim$spreadFirePoints[min1Fire]
    sim$firePolys <- sim$firePolys[min1Fire]
  }
  if (length(sim$firePolys) != length(sim$spreadFirePoints)) {
    stop("mismatched years between firePolys and firePoints")
    #need to implement a better approach that matches each year's IDS
    #these are mostly edge cases if a user passes only one of spreadFirePoints/firePolys
  }

  if (!suppliedElsewhere("ignitionFirePoints", sim)) {
    sim$ignitionFirePoints <- Cache(
      fireSenseUtils::getFirePoints_NFDB_V2,
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      years = P(sim)$fireYears,
      NFDB_pointPath = dPath,
      userTags = c("ignitionFirePoints", P(sim)$.studyAreaName))
    #TODO: what should we set arg redownloadIn to?
  }

  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- prepInputsLCC(destinationPath = dPath,
                                studyArea = sim$studyArea,
                                filename2 = file.path(dPath, paste0("rstLCC_",
                                                                    P(sim)$.studyAreaName,
                                                                    ".tif")),
                                useCache = TRUE)
  }

  if (!suppliedElsewhere("historicalClimateRasters", sim)) {
    stop("please supply sim$historicalClimateRasters")
  }

  if (!suppliedElsewhere("terrainCovariates", sim)) {
     sim$terrainCovariates <- Cache(fireSenseUtils::prepTerrainCovariates,
                                   studyArea = sim$studyArea,
                                   rasterToMatch = sim$rasterToMatch,
                                   destinationPath = dPath,
                                   userTags = c("terrainCovariates"))
  }

  if (!suppliedElsewhere("flammableRTM", sim)) {
    sim$flammableRTM <- LandR::defineFlammable(sim$rstLCC,
                                               nonFlammClasses = P(sim)$nonflammableLCC,
                                               mask = sim$rasterToMatch,
                                               filename2 = file.path(dPath, paste0("flammableRTM_",
                                                                                   P(sim)$.studyAreaName,
                                                                                   ".tif"))
                                               )
  }

  if (!suppliedElsewhere("nonForestedLCCGroups", sim)) {
    sim$nonForestedLCCGroups <- list(
      "nonForest_highFlam" = c(16, 17, 18, 19, 22),
      "nonForest_lowFlam" = c(21, 23, 24, 26, 27, 28, 29, 31)
      #0, 25, 30, 33, 36, 37, 38, 39 non flammable
      # 1:15, 20, 32, 34, 35 forest

    )
  }

  return(invisible(sim))
}
