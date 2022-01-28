defineModule(sim, list(
  name = "fireSense_dataPrepFit",
  description = "",
  keywords = "",
  authors = c(
    person("Ian", "Eddy", role = c("aut", "cre"), email = "ian.eddy@canada.ca"),
    person(c("Alex", "M"), "Chubaty", role = c("ctb"), email = "achubaty@for-cast.ca")
  ),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.4.9003", fireSense_dataPrepFit = "0.0.0.9002"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "fireSense_dataPrepFit.Rmd")),
  reqdPkgs = list("data.table", "fastDummies", "ggplot2", "purrr", "SpaDES.tools",
                  "PredictiveEcology/SpaDES.core@development (>= 1.0.6.9016)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.0.5.9012)",
                  "parallel", "raster", "sf", "sp", "spatialEco", "snow"),
  parameters = bindrows(
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
                    "studyArea name that will be appended to file-backed rasters"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated? This is intended",
                          "for data-type modules, where stochasticity and time are not relevant")),
    defineParameter("areaMultiplier", c("numeric", "function"), fireSenseUtils::multiplier, NA, NA,
                    paste("Either a scalar that will buffer areaMultiplier * fireSize or a function",
                          "of fireSize. Default is 2. See fireSenseUtils::bufferToArea for help")),
    defineParameter("climateGCM", "numeric", "13GCMs_ensemble", NA, NA,
                    paste("Global Circulation Model to use for climate projections:",
                          "currently '13GCMs_ensemble' or 'CanESM5'.")),
    defineParameter("climateSSP", "numeric", 370, NA, NA,
                    "SSP emissions scenario for `climateGCM`: one of 245, 370, or 585."),
    defineParameter("cutoffForYoungAge", class = "numeric", 15, NA, NA,
                    "Age at and below which pixels are considered 'young' --> young <- age <= cutoffForYoungAge"),
    defineParameter("fireYears", "integer", 2001:2019, NA, NA,
                    paste("A numeric vector indicating which years should be extracted",
                          "from the fire databases to use for fitting")),
    defineParameter("forestedLCC", "numeric", c(1:6), NA, NA,
                    "Forested land cover classes. These classes will be excluded from the PCA."),
    defineParameter(name = "igAggFactor", "numeric", 40, 1, NA,
                    desc = "aggregation factor for rasters during ignition prep"),
    defineParameter(name = "ignitionFuelClassCol", class = "character", default = "FuelClass",
                    desc = "the column in sppEquiv that defines unique fuel classes for ignition"),
    defineParameter(name = "minBufferSize", class = "numeric", 5000, NA, NA,
                    desc = "Minimum size of buffer and nonbuffer. This is imposed after multiplier on the bufferToArea fn"),
    defineParameter(name = 'missingLCCgroup', class = 'character', 'nonForest_highFlam', NA, NA,
                    desc = paste("if a pixel is forested but is absent from cohortData, it will be grouped in this class.",
                                 "Must be one of the names in `sim$nonForestedLCCGroups`")),
    defineParameter(name = "nonflammableLCC", class = "numeric", c(13, 16, 17, 18, 19), NA, NA,
                    desc = "non-flammable LCC in sim$rstLCC"),
    defineParameter(name = "PCAcomponentsForClimate", "numeric", 1, 1, NA,
                    desc = "number of PCA components to include from climate variables"),
    defineParameter(name = "PCAcomponentsForTerrain", "numeric", 1, 1, NA,
                    desc = "currently unused - may be needed if using separate terrain and veg PCAs"),
    defineParameter(name = "PCAcomponentsForVeg", "numeric", 10, 1, NA,
                    desc = "number of veg and terrain components to include in GLM"),
    defineParameter(name = "PCAcomponentsFromGLM", "numeric", 5, 0, NA,
                    desc = "the number of components to select from GLM model of burn ~ PCAcomponents" ),
    defineParameter(name = "plotPCA", class = "logical", default = TRUE, NA, NA,
                    desc = "plot the PCA components with a heat map"),
    defineParameter(name = "sppEquivCol", class = "character", default = "LandR", NA, NA,
                    desc = "column name in sppEquiv object that defines unique species in cohortData"),
    defineParameter(name = "spreadFuelClassCol", class = "character", default = "FuelClass",
                    desc = "if using fuel classes for spread, the column in sppEquiv that defines unique fuel classes"),
    defineParameter(name = "useCentroids", class = "logical", default = TRUE,
                    desc = paste("Should fire ignitions start at the sim$firePolygons centroids (TRUE)",
                                 "or at the ignition points in sim$firePoints?")),
    defineParameter(name = "usePCA", class = "logical", default = TRUE, NA, NA,
                    desc = "use PCA approach to covariates, as opposed to fuel class approach"),
    defineParameter(name = "whichModulesToPrepare", class = "character",
                    default = c("fireSense_IgnitionFit", "fireSense_SpreadFit", "fireSense_EscapeFit"),
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
                 desc = "firePolys used to classify timeSinceDisturbance in nonforest LCC"),
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
                 desc = "template raster for study area. Assumes some buffering of core area to limit edge effect of fire."),
    expectsInput(objectName = "rstLCC", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "Raster of land cover. Defaults to LCC05."),
    expectsInput(objectName = "sppEquiv", objectClass = "data.table", sourceURL = NA,
                 desc = "table of LandR species equivalencies"),
    expectsInput(objectName = "standAgeMap2001", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "map of stand age in 2001 used to create cohortData2001"),
    expectsInput(objectName = "standAgeMap2011", objectClass = "RasterLayer", sourceURL = NA,
                 desc = "map of stand age in 2011 used to create cohortData2011"),
    expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", sourceURL = NA,
                 desc = "studyArea that determines spatial boundaries of all data"),
    expectsInput(objectName = "terrainCovariates", objectClass = "RasterStack", sourceURL = NA,
                 desc = "a raster stack of terrain covariates; defaults are elev, aspect, slope, TRI, TWI")
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "climateComponentsToUse", objectClass = "character",
                  desc = "names of the climate components or variables needed for FS models"),
    createsOutput(objectName = "coefficientPrintOut", "data.table",
                  desc = paste("Coefficients from the logit model")),
    createsOutput(objectName = "componentPrintOut", "data.table",
                  desc = paste("A data.table showing the PCA axes and their loadings with the different ",
                               "covariates, e.g., fuel, TPI, HLI, etc.")),
    createsOutput(objectName = "fireBufferedListDT", objectClass = "list",
                  desc = "list of data.tables with fire id, pixelID, and buffer status"),
    createsOutput(objectName = "firePolys", objectClass = "list",
                  desc = "list of spatialPolygonDataFrame objects representing annual fires"),
    createsOutput(objectName = "fireSense_annualSpreadFitCovariates", objectClass = "list",
                  desc = "list of tables with climate PCA components, burn status, polyID, and pixelID"),
    createsOutput(objectName = "fireSense_escapeCovariates", objectClass = "data.table",
                  desc = "ignition covariates with added column of escapes"),
    createsOutput(objectName = "fireSense_escapeFormula", objectClass = "character",
                  desc = "formula for escape, using fuel classes and landcover, as character"),
    createsOutput(objectName = "fireSense_ignitionCovariates", objectClass = "data.table",
                  desc = "table of aggregated ignition covariates with annual ignitions"),
    createsOutput(objectName = "fireSense_ignitionFormula", objectClass = "character",
                  desc = "formula for ignition, using fuel classes and landcover, as character"),
    createsOutput(objectName = "fireSense_nonAnnualSpreadFitCovariates", objectClass = "list",
                  desc = "list of two tables with veg PCA components, burn status, polyID, and pixelID"),
    createsOutput(objectName = "fireSense_spreadFormula", objectClass = "character",
                  desc = "formula for spread, using climate and terrain components, as character"),
    createsOutput(objectName = "ignitionFitRTM", objectClass = "RasterLayer",
                  desc = paste("A (template) raster with information with regards to the spatial resolution and geographical extent of",
                               "fireSense_ignitionCovariates. Used to pass this information onto fireSense_ignitionFitted",
                               "Needs to have number of non-NA cells as attribute (ignitionFitRTM@data@attributes$nonNAs)")),
    # createsOutput(objectName = "fireSense_spreadLogitModel", objectClass = "glm",
    #               desc = "GLM with burn as dependent variable and PCA components as covariates"),
    createsOutput(objectName = "landcoverDT", "data.table",
                  desc = paste("data.table with pixelID and relevant landcover classes",
                               "that is used by predict functions")),
    createsOutput(objectName = "nonForest_timeSinceDisturbance2001", objectClass = "RasterLayer",
                  desc = "time since burn for non-forested pixels in 2001"),
    createsOutput(objectName = "nonForest_timeSinceDisturbance2011", objectClass = "RasterLayer",
                  desc = "time since burn for non-forested pixels in 2011"),
    createsOutput(objectName = "PCAclimate", objectClass = "prcomp",
                  desc = "PCA model for climate covariates, needed for fireSensePredict"),
    createsOutput(objectName = "PCAcoeffPlot", objectClass = "gglot",
                  desc = "ggplot with PCA loadings for axes used to predict spread"),
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

doEvent.fireSense_dataPrepFit = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      if (!all(P(sim)$whichModulesToPrepare %in%
               c("fireSense_SpreadFit", "fireSense_IgnitionFit", "fireSense_EscapeFit"))) {
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

      if (P(sim)$plotPCA) {
        sim <- scheduleEvent(sim, end(sim), "fireSense_dataPrepFit", "plotAndMessage", eventPriority = 9)
      }
      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "cleanUp", eventPriority = 10) #cleans up Mod objects

    },
    prepIgnitionFitData = {
      sim <- prepare_IgnitionFit(sim)
    },
    prepEscapeFitData = {
      sim <- prepare_EscapeFit(sim)
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
  doAssertion <- getOption("fireSenseUtils.assertions", TRUE)

  mod$vegFile <- file.path(outputPath(sim),
                           paste0("fireSense_SpreadFit_veg_coeffs_", P(sim)$.studyAreaName, ".txt"))

  ####prep fire data ####
  if (is.null(sim$firePolys[[1]]$FIRE_ID)) {
    stop("firePolys needs a numeric FIRE_ID column")
  }

  if (!is.numeric(sim$firePolys[[1]]$FIRE_ID)) {
    message("need numeric FIRE_ID column in fire polygons. Coercing to numeric...")
    #this is true of the current NFBB
    origNames <- names(sim$firePolys)
    PointsAndPolys <- lapply(names(sim$firePolys),
                             function(year, polys = sim$firePolys, points = sim$spreadFirePoints) {
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

  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1, length(sim$firePolys))
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
  pointsIDColumn <- "FIRE_ID"
  sim$spreadFirePoints <- Cache(harmonizeBufferAndPoints,
                                cent = sim$spreadFirePoints,
                                buff = fireBufferedListDT,
                                ras = sim$flammableRTM,
                                idCol = pointsIDColumn, #this is different from default
                                userTags = c("harmonizeBufferAndPoints"))

  ## Also 2 other problems:
  ## 1. Big fire, but ignition is in non-flammable pixels e.g., lake -- bad;
  ##    solution -- pick nearest pixel in burned polygon
  ## 2. Small fire, ignition in non-flammable pixel, but NO pixel in burned polygon
  ##    is actually flammable -- remove this from data
  out22 <- Map(fp = sim$spreadFirePoints, fpoly = fireBufferedListDT, function(fp, fpoly) {
    isFlammable <- raster::extract(sim$flammableRTM, fp)
    isFlammable[which(is.na(isFlammable))] <- FALSE
    if (any(!isFlammable)) {
      badStarts <- fp[!isFlammable,][[pointsIDColumn]]
      badStartsPixels <- cbind(ids = badStarts,
                               raster::extract(sim$flammableRTM, fp[!isFlammable,], cellnumbers = TRUE))
      polysWBadStarts <- fpoly[ids %in% badStarts,]
      cells <- polysWBadStarts[polysWBadStarts$buffer == 1, ]
      flammableInPolys <- sim$flammableRTM[][cells$pixelID] == 1
      cells <- cells[flammableInPolys,]
      rmFireIDs <- setdiff(badStarts, unique(cells$ids))
      newSp <- numeric()
      if (any(flammableInPolys)) {
        xyPolys <- cbind(id = cells$ids,
                         pixelID = cells$pixelID,
                         raster::xyFromCell(sim$flammableRTM, cells$pixelID))
        xyPoints <- cbind(id = badStartsPixels[, "ids"],
                          #pixelID = badStartsPixels[, "cells"],
                          raster::xyFromCell(sim$flammableRTM, badStartsPixels[, "cells"]))
        dd <- as.data.table(distanceFromEachPoint(to = xyPolys, from = xyPoints))
        nearestPixels <- dd[, .SD[which.min(dists)], by = "id"]
        idsToChange <- unique(nearestPixels$id)

        # Rm bad points that are just not in fires
        # 1. create new SpatialPointsDataFrame with shifted coordinates
        df <- as.data.frame(fp[fp[[pointsIDColumn]] %in% idsToChange,])
        df <- df[, !colnames(df) %in% c("x", "y")]
        newSp <- SpatialPointsDataFrame(nearestPixels[, .(x,y)], data = df, proj4string = crs(fp))
      }
      # 2. rm bad points
      fp <- fp[!fp[[pointsIDColumn]] %in% badStarts,]
      # 3. rbind these two together
      if (length(newSp))
        fp <- rbind(fp, newSp)

      fpoly <- fpoly[!fpoly$ids %in% rmFireIDs,]
    }
    list(SpatialPoints = fp, FireBuffered = fpoly)

  })
  out22 <- purrr::transpose(out22)
  sim$spreadFirePoints <- out22$SpatialPoints
  fireBufferedListDT <- out22$FireBuffered

  #### prep terrain for PCA ####
  sim$terrainDT <- buildTerrainDT(terrainCovariates = sim$terrainCovariates,
                                  flammableRTM = sim$flammableRTM)

  #### prep landcover for PCA ####
  #this is needed whether PCA or not
  sim$landcoverDT <- makeLandcoverDT(rstLCC = sim$rstLCC,
                                     flammableRTM = sim$flammableRTM,
                                     forestedLCC = P(sim)$forestedLCC,
                                     nonForestedLCCGroups = sim$nonForestedLCCGroups)

  # cannot merge because before subsetting due to column differences over time

  mod$firePolysForAge <- lapply(sim$firePolysForAge[lengths(sim$firePolysForAge) > 0],
                                FUN = sf::st_as_sf) %>%
    lapply(., FUN = function(x) {
        x <- x[, "YEAR"]
      }) %>%
    do.call(rbind, .)

  #fuelClasses only uses this object for comparing the logistic model
  #otherwise it could be inside the if
  origDTThreads <- data.table::setDTthreads(2)
  cohorts2001 <-  Cache(castCohortData,
                        cohortData = sim$cohortData2001,
                        pixelGroupMap = sim$pixelGroupMap2001,
                        year = 2001,
                        ageMap = NULL,
                        cutoffForYoungAge = P(sim)$cutoffForYoungAge,
                        lcc = sim$landcoverDT,
                        terrainDT = sim$terrainDT,
                        missingLCC = P(sim)$missingLCC)
  cohorts2011 <-  Cache(castCohortData, cohortData = sim$cohortData2011,
                        pixelGroupMap = sim$pixelGroupMap2011,
                        year = 2011,
                        cutoffForYoungAge = P(sim)$cutoffForYoungAge,
                        ageMap = NULL,
                        lcc = sim$landcoverDT,
                        terrainDT = sim$terrainDT,
                        missingLCC = P(sim)$missingLCC)

  vegPCAdat <- rbindlist(list(cohorts2001, cohorts2011), use.names = TRUE)
  mod$vegPCAdat <- vegPCAdat

  if (P(sim)$usePCA) {
    nonTreeNames <- c(names(sim$landcoverDT), names(sim$terrainDT))
    nonTreeNames <- nonTreeNames[!nonTreeNames %in% "pixelID"]
    # Next line can't handle memoising well
    opts <- options("reproducible.useMemoise" = FALSE)
    vegList <- Cache(makeVegTerrainPCA,
                     dontAlter = nonTreeNames, #in case default lcc change
                     dataForPCA = vegPCAdat)
    options(opts)
    vegComponents <- vegList$vegComponents

    # To see how much variance explained https://stats.stackexchange.com/questions/254592/calculating-pca-variance-explained/254598
    sim$PCAveg <- vegList$vegTerrainPCA
    if (isTRUE(doAssertion))  {
      ss <- summary(sim$PCAveg)$importance
      ss <- data.frame("component" = rownames(ss), ss)
      reproducible::messageDF(ss, round = 3, colour = "green")
    }
    #if there are fewer components than parameter, take them all
    components <- paste0("PC", 1:min(as.integer(P(sim)$PCAcomponentsForVeg), length(sim$PCAveg$scale)))

    removeCols <- setdiff(colnames(vegComponents), c(components, "pixelID", "year", "youngAge"))
    if (length(removeCols))
      set(vegComponents, NULL, removeCols, NULL)
    # vegComponents <- vegComponents[, .SD, .SDcols = c(components, "pixelID", "year", "youngAge")]
    #rename components so climate/veg components distinguishable
    setnames(vegComponents, old = components, new = paste0("veg", components))

  } else {
    #fuel class approach
    if (any(is.na(sim$sppEquiv[[P(sim)$spreadFuelClassCol]]))) {
      stop("All species must have a spread fuel class defined in sppEquiv.")
    }

    vegComponents <- Map(f = cohortsToFuelClasses,
                         cohortData = list(sim$cohortData2001, sim$cohortData2011),
                         yearCohort = list(2001, 2011),
                         pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                         MoreArgs = list(sppEquiv = sim$sppEquiv,
                                         sppEquivCol = P(sim)$sppEquivCol,
                                         flammableRTM = sim$flammableRTM,
                                         fuelClassCol = P(sim)$spreadFuelClassCol,
                                         cutoffForYoungAge = -1)) #youngAge will be calculated every year

    vegComponents <- lapply(vegComponents, FUN = function(x){
      dt <- as.data.table(getValues(x))
      dt[, pixelID := 1:ncell(x)]
      # dt <- dt[!is.na(youngAge)] unlike in prepare_ignition, we do not remove youngAge
      #it will be added annually by spreadFit
      return(dt)
    })

    vegComponents[[1]][, year := 2001]
    vegComponents[[2]][, year := 2011]

    #joining with landcoverDT eliminates non-flammable pixels - every remaining pixel MUST have a value
    #join landcover separately, as the forest extent changes between 2001 and 2011
    vegComponents[[1]] <- vegComponents[[1]][sim$landcoverDT, on = c("pixelID")]
    vegComponents[[2]] <- vegComponents[[2]][sim$landcoverDT, on = c("pixelID")]
    vegComponents <- rbindlist(vegComponents)
    setcolorder(vegComponents, c("pixelID", "year"))

    #set 'orphaned' pixels as P(sim)$missingLCC - the forest that is not LandR forest
    lccNames <- setdiff(names(vegComponents), c("pixelID", "year"))
    vegComponents[, missingLCC := rowSums(vegComponents[, .SD, .SDcols = lccNames])]
    vegComponents[missingLCC == 0, eval(P(sim)$missingLCC) := 1]
    #TODO: when we add assertions, assert that there are no rows where missingLCC = 2
    vegComponents[, missingLCC := NULL]
  }

  # The next 1 line replaces the 8 or so lines after
  fireBufferedListDT <- Cache(rmMissingPixels, fireBufferedListDT, vegComponents$pixelID)
  ####prep Climate components####
  flammableIndex <- data.table(index = 1:ncell(sim$flammableRTM), value = getValues(sim$flammableRTM)) %>%
    .[value == 1,] %>%
    .$index

  #this involves a large object (27 years of raster data converted to 3 columns), year, pixelID, value
  #need it in one dt for PCA, but it is less efficient so we convert back to list of annual dts
  climatePCAdat <- Cache(climateRasterToDataTable,
                         historicalClimateRasters = sim$historicalClimateRasters,
                         Index = flammableIndex, userTags = c("climateRasterToDataTable"))
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
    #don't rename, don't rescale
    climateComponents <- climatePCAdat[[1]]
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
  vars <- names(fireSenseVegData)[!names(fireSenseVegData) %in% c("year", "pixelID", "youngAge", "burned", "ids")]


  grepCol <- ifelse(P(sim)$usePCA, "veg*", paste0(vars, collapse = "|"))
  logitFormula <- grep(grepCol, names(fireSenseVegData), value = TRUE) %>%
    paste0(., collapse = " + ") %>%
    paste0("burned ~ ", .)

  if (isTRUE(doAssertion)) {
    ttt <- table(fireSenseVegData$burned)
    ratioZeroToOne <- ttt[1]/ttt[2]
    if (ratioZeroToOne < 5) {
      stop("The number of pixels in the fire buffers should be at least 5x the number of burned pixels\n",
           "Please create larger buffers around fires in fireBufferedListDT, e.g., via ",
           "fireSenseUtils::bufferToArea(..., areaMultiplier = multiplier)")
    }

    fireSenseLogit <- glm(formula = eval(logitFormula),
                          data = fireSenseVegData,
                          family = "binomial",
                          na.action = na.exclude)

    ## create vegFile outputs
    #this is a "null model" to compare with the PCA approach -
    fb1 <- logisticCovariatesPre2011[, .(pixelID, year, buffer)][vegPCAdat, on = c("pixelID", "year")]
    fb2 <- logisticCovariatesPost2011[, .(pixelID, year, buffer)][vegPCAdat, on = c("pixelID", "year")]
    vv <- rbind(fb1, fb2)
    #left join - so must remove the unbuffered pixels
    vv <- vv[!is.na(buffer)]

    set(vv, NULL, c("year", "pixelID", "pixelGroup", "youngAge"), NULL)
    gg <- glm(buffer ~ ., data = vv, family = "binomial", na.action = na.exclude)
    ggs <- summary(gg)
    coefs1 <- ggs$coefficients
    coefs1 <- data.frame("term" = rownames(coefs1), coefs1)
    pseudoR2_vegDirect <- 1 - gg$deviance / gg$null.deviance
    pseudoR2_vegPCA <- 1 - fireSenseLogit$deviance / fireSenseLogit$null.deviance

    ## print to screen
    message("Vegetation model direct: ")
    messageDF(coefs1)
    message("R2 with vegetation direct: ", round(pseudoR2_vegDirect, 3))
    message("R2 with fireSense version: ", round(pseudoR2_vegPCA, 3))

    ## print to file
    cat(paste("Vegetation model direct:\n"), file = mod$vegFile, sep = "\n", append = FALSE)
    cat(capture.output(coefs1), file = mod$vegFile, sep = "\n", append = TRUE)
    cat(paste("R2 with vegetation direct: ", round(pseudoR2_vegDirect, 3)),
        file = mod$vegFile, sep = "\n", append = TRUE)
    cat(paste("R2 with fireSense version: ", round(pseudoR2_vegPCA, 3), "\n"),
        file = mod$vegFile, sep = "\n", append = TRUE)
  }

  #take largest coeffiecients as they are mean-centered and scaled, number determined by param
  bestComponents <- sort(abs(fireSenseLogit$coefficients[2:length(fireSenseLogit$coefficients)]),
                         decreasing = TRUE)[1:P(sim)$PCAcomponentsFromGLM]
  if (P(sim)$usePCA) {
    bestComponents <- sort(abs(fireSenseLogit$coefficients[2:length(fireSenseLogit$coefficients)]),
                           decreasing = TRUE)[1:P(sim)$PCAcomponentsFromGLM]
    sim$vegComponentsToUse <- names(bestComponents) #the values will be wrong due to abs, just take names
  } else {
   sim$vegComponentsToUse <-  setdiff(names(fireSenseLogit$coefficients), "(Intercept)")
  }
  # mod$fireSense_spreadLogitModel <- fireSenseLogit
  sim$coefficientPrintOut <- fireSenseLogit$coefficients[sim$vegComponentsToUse]

  RHS <- paste0("youngAge + ",
                paste(sim$climateComponentsToUse, collapse = " + "),
                " + ",
                paste(sim$vegComponentsToUse, collapse = " + "))

  mod$fireSenseVegData <- fireSenseVegData
  mod$climateComponents <- climateComponents #needed by prep spread
  #sim$fireSense_spreadFormula <- as.formula(paste("~0 + ", paste(RHS)))
  sim$fireSense_spreadFormula <- paste("~0 +", paste(RHS))
  sim$fireBufferedListDT <- fireBufferedListDT #needed by DEOptim

  return(invisible(sim))
}

prepare_SpreadFit <- function(sim) {

  ## Put in format for DEOptim
  ## Prepare annual spread fit covariates
  ## this is a funny way to get years but avoids years with 0 fires
  years <- paste0("year", P(sim)$fireYears)
  yearsWithFire <- years[paste0("year", P(sim)$fireYears) %in% names(sim$firePolys)]
  pre2011 <- yearsWithFire[yearsWithFire %in% paste0("year", min(P(sim)$fireYears):2010)]
  post2011 <- yearsWithFire[yearsWithFire %in% paste0("year", 2011:max(P(sim)$fireYears))]

  #currently there are NAs in climate due to non-flammable pixels in fire buffer
  # names(years) <- years

  fbl <- rbindlist(sim$fireBufferedListDT, idcol = "year")
  rmCols <- setdiff(colnames(fbl), c("pixelID", "year"))
  set(fbl, NULL, rmCols, NULL)
  fbl <- mod$climateComponents[fbl, on = c("year", "pixelID"), nomatch = NULL]
  fireSense_annualSpreadFitCovariates <- split(fbl, by = "year", keep.by = FALSE)

  #prepare non-annual spread fit covariates
  pre2011Indices <- sim$fireBufferedListDT[names(sim$fireBufferedListDT) %in% pre2011]
  post2011Indices <- sim$fireBufferedListDT[!names(sim$fireBufferedListDT) %in% pre2011]

  colsToExtract <- c("pixelID", sim$vegComponentsToUse)

  #remove age from the annual covariates - it will come from TSD
  nonAnnualPre2011 <- mod$fireSenseVegData[year < 2011, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID),]

  nonAnnualPost2011 <- mod$fireSenseVegData[year >= 2011, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID)] # remove duplicates from same pixel diff year

  # Create one universal TSD map for each initial time period combining stand age/ time since burn
  TSD2001 <- makeTSD(year = 2001, firePolys = sim$firePolysForAge,
                     standAgeMap = sim$standAgeMap2001, lcc = sim$landcoverDT,
                     cutoffForYoungAge = P(sim)$cutoffForYoungAge)
  TSD2011 <- makeTSD(year = 2011, firePolys = sim$firePolysForAge,
                     standAgeMap = sim$standAgeMap2011,
                     lcc = sim$landcoverDT,
                     cutoffForYoungAge = P(sim)$cutoffForYoungAge)
  #the function will do this below, and then use the data.table with location of non-forest to fill in those ages
  # pmap allows for internal debugging when there are large lists that are passed in; Map does not
  annualCovariates <- Cache(purrr::pmap, .l = list(
    years = list(c(2001:2010), c(2011:max(P(sim)$fireYears))),
    annualCovariates = list(fireSense_annualSpreadFitCovariates[pre2011],
                            fireSense_annualSpreadFitCovariates[post2011]),
    standAgeMap = list(TSD2001, TSD2011)
  ), .f = calcYoungAge, fireBufferedListDT = sim$fireBufferedListDT,
  cutoffForYoungAge = P(sim)$cutoffForYoungAge)

  sim$fireSense_annualSpreadFitCovariates <- do.call(c, annualCovariates)

  if (start(sim) != 2011) {
    stop("This module is assuming that the start year is 2011. ",
         "It must start in 2011, or some of the assumptions would be incorrect")
  }
  sim$nonForest_timeSinceDisturbance2011 <- TSD2011
  sim$nonForest_timeSinceDisturbance2001 <- TSD2001

  sim$fireSense_nonAnnualSpreadFitCovariates <- list(nonAnnualPre2011, nonAnnualPost2011)
  names(sim$fireSense_nonAnnualSpreadFitCovariates) <- c(paste(names(pre2011Indices), collapse = "_"),
                                                         paste(names(post2011Indices), collapse = "_"))

  return(invisible(sim))
}

prepare_IgnitionFit <- function(sim) {

  if (any(is.na(sim$sppEquiv[[P(sim)$ignitionFuelClassCol]]))) {
    stop("All species must have an ignition fuelClass defined.")
  }

  #correct ignitions that fall on non-flammable pixels
  #if aggregating, still seems like it is an important step

  #account for forested pixels that aren't in cohortData
  sim$landcoverDT[, rowSums := rowSums(.SD), .SD = setdiff(names(sim$landcoverDT), "pixelID")]
  forestPix <- sim$landcoverDT[rowSums == 0,]$pixelID
  problemPix2001 <- forestPix[is.na(sim$pixelGroupMap2001[forestPix])]
  problemPix2011 <- forestPix[is.na(sim$pixelGroupMap2011[forestPix])]
  set(sim$landcoverDT, NULL, 'rowSums', NULL)

  #The non-forests aren't the same between years, due to cohortData being different
  landcoverDT2001 <- copy(sim$landcoverDT)
  landcoverDT2001[pixelID %in% problemPix2001, eval(P(sim)$missingLCC) := 1]
  landcoverDT2011 <- copy(sim$landcoverDT)
  landcoverDT2011[pixelID %in% problemPix2011, eval(P(sim)$missingLCC) := 1]
  #TODO Work this into assertions
  #all 'problem pixels' should be forest cover classes in LCC raster
  #all(unique(sim$rstLCC[problemPix2001] %in% P(sim)$forestedLCC))
  #all(unique(sim$rstLCC[problemPix2011] %in% P(sim)$forestedLCC))
  #first put landcover into raster stack


  #non-flammable pixels require zero values for non-forest landcover, not NA
  putBackIntoRaster <- function(lcc, landcoverDT, flammableMap) {
    lccRasters <- list()
    for (i in 1:length(lcc)) {
      lccRasters[[i]] <- raster(flammableMap)
      lccRasters[[i]][!is.na(flammableMap)] <- 0 #all non-NA pixels must be 0
      lccRasters[[i]][landcoverDT$pixelID] <- landcoverDT[, get(lcc[i])]
    }
    names(lccRasters) <- lcc
    return(lccRasters)
  }

  LCCras <- Cache(Map,
                  f = putBackIntoRaster,
                  landcoverDT = list(landcoverDT2001,landcoverDT2011),
                  MoreArgs = list(lcc = names(sim$nonForestedLCCGroups),
                                  flammableMap = sim$flammableRTM),
                  userTags = c("putBackIntoRaster"))  %>%
    lapply(., FUN = brick) %>%
    lapply(., aggregate, fact = P(sim)$igAggFactor, fun = mean)
  names(LCCras) <- c("year2001", "year2011")

  fuelClasses <- Map(f = cohortsToFuelClasses,
                     cohortData = list(sim$cohortData2001, sim$cohortData2011),
                     yearCohort = list(2001, 2011),
                     pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                     MoreArgs = list(sppEquiv = sim$sppEquiv,
                                     sppEquivCol = P(sim)$sppEquivCol,
                                     flammableRTM = sim$flammableRTM,
                                     fuelClassCol = P(sim)$ignitionFuelClassCol,
                                     cutoffForYoungAge = P(sim)$cutoffForYoungAge)) %>%
    lapply(., FUN = raster::brick) %>%
    lapply(., aggregate, fact = P(sim)$igAggFactor, fun = mean)
  names(fuelClasses) <- c("year2001", "year2011")

  climate <- sim$historicalClimateRasters
  if (length(climate) > 1) {
    stop("need to fix ignition for multiple climate variables. contact module developers")
    #for now - fix when priority
  } else {
    climVar <- names(climate)
    climate <- raster::stack(climate[[1]]) %>%
      aggregate(., fact = P(sim)$igAggFactor, fun = mean)
  }
  #ignition won't have same years as spread so we do not use names of init objects
  #The reason is some years may have no significant fires, e.g. 2001 in RIA
  pre2011 <- paste0("year", min(P(sim)$fireYears):2010)
  post2011 <- paste0("year", 2011:max(P(sim)$fireYears))

  compareRaster(fuelClasses$year2001, fuelClasses$year2011, LCCras$year2001, LCCras$year2011, climate)
  #this is joining fuel class, LCC, and climate, subsetting to flamIndex, calculating n of ignitions
  fireSense_ignitionCovariates <- Map(f = fireSenseUtils::stackAndExtract,
                                      years = list(pre2011, post2011),
                                      fuel = list(fuelClasses$year2001, fuelClasses$year2011),
                                      LCC = list(LCCras$year2001, LCCras$year2011),
                                      MoreArgs = list(climate = climate,
                                                      fires = sim$ignitionFirePoints,
                                                      climVar = climVar
                                      ))

  fireSense_ignitionCovariates <- rbindlist(fireSense_ignitionCovariates)

  #remove any pixels that are 0 for all classes
  fireSense_ignitionCovariates[, coverSums := rowSums(.SD), .SD = setdiff(names(fireSense_ignitionCovariates),
                                                                          c("MDC", "cells", "ignitions", "year"))]

  fireSense_ignitionCovariates <- fireSense_ignitionCovariates[coverSums > 0]
  set(fireSense_ignitionCovariates, NULL, "coverSums", NULL)

  #rename cells to pixelID - though aggregated raster is not saved
  setnames(fireSense_ignitionCovariates, old = "cells", new = "pixelID")
  fireSense_ignitionCovariates[, year := as.numeric(year)]
  setcolorder(fireSense_ignitionCovariates, neworder = c("pixelID", "ignitions", climVar, 'youngAge'))
  sim$fireSense_ignitionCovariates <- as.data.frame(fireSense_ignitionCovariates) #avoid potential conflict in ignition


  #make new ignition object, ignitionFitRTM
  sim$ignitionFitRTM <- raster(fuelClasses$year2001)
  sim$ignitionFitRTM@data@attributes$nonNAs <- nrow(sim$fireSense_ignitionCovariates)

  #build formula
  igCovariates <- names(sim$fireSense_ignitionCovariates)
  igCovariates <- igCovariates[!igCovariates %in% c(climVar, "year", "ignitions", "pixelID")]
  pwNames <- abbreviate(igCovariates, minlength = 3, use.classes = TRUE, strict = FALSE)
  interactions <- paste0(igCovariates, ":", climVar)
  pw <- paste0(igCovariates, ":", "pw(", climVar, ", k_", pwNames, ")")
  #sanity check for base::abbreviate
  if (!all(length(unique(pw)), length(unique(interactions)) == length(igCovariates))) {
    warning("automated ignition formula construction needs review")
  }
  sim$fireSense_ignitionFormula <- paste0("ignitions ~ ", paste0(interactions, collapse = " + "), " + ",
                                          paste0(pw, collapse  = " + "), "- 1")

  return(invisible(sim))
}

prepare_EscapeFit <- function(sim) {
  if (is.null(sim$fireSense_ignitionCovariates)) {
    #the datasets are essentially the same, with one column difference
    stop("Please include ignitionFit in parameter 'whichModulesToPrepare' if running EscapeFit")
  }

  escapeThreshHa <- prod(res(sim$flammableRTM))/10000
  escapes <- sim$ignitionFirePoints[sim$ignitionFirePoints$SIZE_HA > escapeThreshHa,]

  #make a template aggregated raster - values are irrelevant, only need pixelID
  aggregatedRas <- aggregate(sim$historicalClimateRasters[[1]][[1]],
                             fact = P(sim)$igAggFactor, fun = mean)

  escapeDT <- raster::extract(aggregatedRas, escapes, cellnumber = TRUE) %>%
    as.data.table(.) %>%
    .[, year := escapes$YEAR] %>%
    .[, .(year, cells)] %>%
    .[, .(.N), .(year, cells)] %>%
    setnames(., c("N", "cells"), c("escapes", "pixelID"))

  escapeDT[, year := as.numeric(year)]
  escapeDT <- escapeDT[sim$fireSense_ignitionCovariates, on = c("pixelID", "year")]
  escapeDT[is.na(escapes), escapes := 0]

  sim$fireSense_escapeCovariates <- escapeDT


  escapeVars <- names(escapeDT)[!names(escapeDT) %in% c("year", "pixelID", "escapes", "ignitions")]
  LHS <- paste0("cbind(escapes, ignitions - escapes) ~ ")
  RHS <- paste0(escapeVars, collapse = " + ")
  sim$fireSense_escapeFormula <- paste0(LHS, RHS, " - 1")

  return(invisible(sim))
}

cleanUpMod <- function(sim) {
  ## because prep for ignition/escape/spread isn't necessarily run, clean up here
  mod$climateComponents <- NULL # remove for memory sake
  mod$firePolysForAge <- NULL
  mod$fireSenseVegData <- NULL
  mod$vegPCAdat <- NULL

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
  if (P(sim)$usePCA){
  checkPath(file.path(outputPath(sim), "figures"), create = TRUE)
  components <- as.data.table(sim$PCAveg$rotation)
  setnames(components, old = colnames(components), new = paste0("veg", colnames(components)))
  components[, covariate := row.names(sim$PCAveg$rotation)]
  setcolorder(components, neworder = "covariate")
  sim$componentPrintOut <- components[, .SD, .SDcol = c("covariate", sim$vegComponentsToUse)]

  ## to screen
  messageDF(sim$componentPrintOut)

  ## to file
  cat(capture.output(sim$componentPrintOut), file = mod$vegFile, sep = "\n", append = TRUE)

  components <- melt.data.table(data = sim$componentPrintOut,
                                id.vars = c("covariate"),
                                variable.name = "component",
                                value.name = "loading")
  components[, loading := round(loading, digits = 2)]

  # coefficientPrintOut <- fireSenseLogit$coefficients[sim$vegComponentsToUse]
  # coefficientPrintOut <- mod$fireSense_spreadLogitModel$coefficients[sim$vegComponentsToUse]
  coefficientSigns <- data.table(val = sim$coefficientPrintOut, component = names(sim$coefficientPrintOut))
  components <- coefficientSigns[components, on = "component"]
  components[, sign := ifelse(val < 0, "-", "+")]
  components[, component := paste0(component, sign)]

  sim$PCAcoeffPlot <- ggplot(data = components, aes(x = component, y = covariate, fill = loading)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    geom_text(data = components, label = components$loading) +
    ggtitle(paste("loading of components most correlated with fire for", P(sim)$.studyAreaName))

  if (!is.na(P(sim)$.plotInitialTime)) {
    plot(sim$PCAcoeffPlot)
  }
  ggsave(file.path(outputPath(sim), "figures", paste0("PCAcoeffLoadings_", P(sim)$.studyAreaName, ".png")),
         sim$PCAcoeffPlot)
  }
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("studyArea", sim)) {
    stop("Please supply study area - this object is key")
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- LandR::prepInputsStandAgeMap(studyArea = sim$studyArea,
                                                      destinationPath = dPath)
  }

  if (!all(suppliedElsewhere("cohortData2011", sim),
           suppliedElsewhere("pixelGroupMap2011", sim),
           suppliedElsewhere("pixelGroupMap2001", sim),
           suppliedElsewhere("cohortData2001", sim))) {
    stop("Stop - need cohortData and pixelGroupMap objects - contact module creators")
  }

  if (!suppliedElsewhere("firePolys", sim) | !suppliedElsewhere("firePolysForAge", sim)) {
    browser()
    # don't want to needlessly postProcess the same firePolys objects
    allFirePolys <- Cache(fireSenseUtils::getFirePolygons,
                          years = c(min(P(sim)$fireYears - P(sim)$cutoffForYoungAge):max(P(sim)$fireYears)),
                          studyArea = sim$studyArea,
                          destinationPath = dPath,
                          useInnerCache = TRUE,
                          userTags = c("firePolys", paste0("years:", range(P(sim)$fireYears))))
  }

  if (!suppliedElsewhere("firePolys", sim)) {
    sim$firePolys <- allFirePolys[names(allFirePolys) %in% paste0("year", P(sim)$fireYears)]
  }

  if (!suppliedElsewhere("firePolysForAge", sim)) {
    sim$firePolysForAge <- allFirePolys
  }

  if (!suppliedElsewhere("standAgeMap2001", sim)) {
    sim$standAgeMap2001 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch,
                                 studyArea = sim$studyArea,
                                 destinationPath = dPath,
                                 filename2 = "standAgeMap2001.tif",
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
        cent <- SpatialPointsDataFrame(centCoords, as.data.frame(ras))
        return(cent)
      }
    }
    # suppress any startup messages
    mc <- pemisc::optimalClusterNum(2e3, maxNumClusters = length(sim$firePolys))
    clObj <- parallel::makeCluster(type = "SOCK", mc)

    # library of LandR is a workaround to deal with weird
    #  Function found when exporting methods from the namespace 'raster' which is not S4 generic: 'all.equal'
    a <- parallel::clusterEvalQ(cl = clObj, {library(LandR); library(raster); library(rgeos)})
    clusterExport(cl = clObj, list("firePolys"), envir = sim)
    # debug(reproducible:::dealWithRasters)
    # debug(reproducible:::dealWithRastersOnRecovery)
    # on.exit({
    #   undebug(reproducible:::dealWithRasters)
    #   undebug(reproducible:::dealWithRastersOnRecovery)
    # }, add = TRUE)
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
    ignitionFirePoints <- Cache(
      fireSenseUtils::getFirePoints_NFDB_V2,
      studyArea = sim$studyArea,
      rasterToMatch = sim$rasterToMatch,
      years = P(sim)$fireYears,
      NFDB_pointPath = dPath,
      userTags = c("ignitionFirePoints", P(sim)$.studyAreaName),
      plot = !is.na(P(sim)$.plotInitialTime)
    ) ## TODO: what should we set arg redownloadIn to?

    sim$ignitionFirePoints <- ignitionFirePoints[ignitionFirePoints$CAUSE == "L",]
  }

  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- prepInputsLCC(
      year = 2010,
      destinationPath = dPath,
      studyArea = sim$studyArea,
      filename2 = file.path(dPath, paste0("rstLCC_", P(sim)$.studyAreaName, ".tif")),
      useCache = TRUE)
  }

  if (!suppliedElsewhere("historicalClimateRasters", sim)) {
    stop("please supply sim$historicalClimateRasters")
  }

  if (!suppliedElsewhere("terrainCovariates", sim)) {
     sim$terrainCovariates <- Cache(
       fireSenseUtils::prepTerrainCovariates,
       studyArea = sim$studyArea,
       rasterToMatch = sim$rasterToMatch,
       destinationPath = dPath,
       userTags = c("terrainCovariates")
      )
  }

  if (!suppliedElsewhere("flammableRTM", sim)) {
    sim$rstLCC[] <- as.integer(sim$rstLCC[])
    sim$flammableRTM <- LandR::defineFlammable(
      sim$rstLCC,
      nonFlammClasses = P(sim)$nonflammableLCC,
      mask = sim$rasterToMatch,
      filename2 = file.path(dPath, paste0("flammableRTM_", P(sim)$.studyAreaName, ".tif"))
    )
  }

  if (!suppliedElsewhere("nonForestedLCCGroups", sim)) {
    sim$nonForestedLCCGroups <- list(
      "nonForest_highFlam" = c(8, 10, 14),#shrubland, grassland, wetland
      "nonForest_lowFlam" = c(11, 12, 15) #shrub-lichen-moss + cropland. 2 barren classes are nonflam
    )
  }

  return(invisible(sim))
}

rmMissingPixels <- function(fbldt, pixelIDsAllowed)  {
  fbldt <- rbindlist(fbldt, idcol = "year")
  fbldt <- fbldt[pixelID %in% unique(pixelIDsAllowed)]
  fireBufferedListDT <- split(fbldt, by = "year", keep.by = FALSE)
}
