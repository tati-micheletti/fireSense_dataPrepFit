defineModule(sim, list(
  name = "fireSense_dataPrepFit",
  description = "Prepare data required by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit`.",
  keywords = "fireSense",
  authors = c(
    person("Ian", "Eddy", role = c("aut", "cre"), email = "ian.eddy@nrcan-rncan.gc.ca"),
    person(c("Alex", "M"), "Chubaty", role = c("ctb"), email = "achubaty@for-cast.ca")
  ),
  childModules = character(0),
  version = list(SpaDES.core = "1.0.4.9003", fireSense_dataPrepFit = "0.0.0.9001"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.txt", "fireSense_dataPrepFit.Rmd")),
  reqdPkgs = list("data.table", "fastDummies", "ggplot2", "purrr", "SpaDES.tools",
                  "PredictiveEcology/SpaDES.core@development (>= 1.0.6.9016)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.0.5.9043)",
                  "parallel", "raster", "sf", "sp", "spatialEco", "snow"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("areaMultiplier", c("numeric", "function"), fireSenseUtils::multiplier, NA, NA,
                    paste("Either a scalar that will buffer areaMultiplier * fireSize or a function",
                          "of fireSize. Default is 2. See `fireSenseUtils::bufferToArea` for help")),
    defineParameter("bufferForFireRaster", "numeric", 1000, 0, NA,
                    paste("The distance that determine whether separate patches of burned pixels originated",
                          "from the same fire. Only relevant when 'useRasterizedFireForSpread' is TRUE.",
                          "This param is separate from 'minBufferSize', which is used to determine the",
                          "minimum sample of burned and unburned pixels to include in each fire.")),
    defineParameter("cutoffForYoungAge", "numeric", 15, NA, NA,
                    "Age at and below which pixels are considered 'young' --> young <- age <= cutoffForYoungAge"),
    defineParameter("fireYears", "integer", 2001:2020, NA, NA,
                    paste("A numeric vector indicating which years should be extracted",
                          "from the fire databases to use for fitting")),
    defineParameter("forestedLCC", "numeric", c(1:6), NA, NA,
                    "Forested land cover classes. These classes will be excluded from the PCA."),
    defineParameter("igAggFactor", "numeric", 40, 1, NA,
                    "aggregation factor for rasters during ignition prep."),
    defineParameter("ignitionFuelClassCol", "character", "FuelClass", NA, NA,
                    "the column in `sppEquiv` that defines unique fuel classes for ignition"),
    defineParameter("minBufferSize", "numeric", 5000, NA, NA,
                    paste("Minimum number of cells in buffer and nonbuffer. This is imposed after the",
                          "multiplier on the `bufferToArea` fn")),
    defineParameter("missingLCCgroup", "character", "nonForest_highFlam", NA, NA,
                    paste("if a pixel is forested but is absent from `cohortData`, it will be grouped in this class.",
                          "Must be one of the names in `sim$nonForestedLCCGroups`")),
    defineParameter("nonflammableLCC", "numeric", c(13, 16, 17, 18, 19), NA, NA,
                    "non-flammable LCC in `sim$rstLCC`."),
    defineParameter(name = "nonForestCanBeYoungAge", class = "logical", TRUE, NA, NA,
                    "if TRUE, burned non-forest will be treated as youngAge"),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "column name in sppEquiv object that defines unique species in cohortData"),
    defineParameter("spreadFuelClassCol", "character", "FuelClass", NA, NA,
                    "if using fuel classes for spread, the column in `sppEquiv` that defines unique fuel classes"),
    defineParameter("useCentroids", "logical", TRUE, NA, NA,
                    paste("Should fire ignitions start at the `sim$firePolygons` centroids (TRUE)",
                          "or at the ignition points in `sim$firePoints`?")),
    defineParameter("useRasterizedFireForSpread", "logical", FALSE, NA, NA,
                    paste("Should rasterized fire be used in place of a vectorized fire dataset?",
                          "This method attributes burned pixels to specific fires,",
                          "only examines the latest fire in a pixel, and may be subject to temporal error.",
                          "is therefore more appropriate in areas with low rates of fire,",
                          "or where the NFDB dataset may be incomplete (ie northern Ontario).")),
    defineParameter("whichModulesToPrepare", "character",
                    c("fireSense_IgnitionFit", "fireSense_SpreadFit", "fireSense_EscapeFit"),
                    NA, NA, "Which fireSense fit modules to prep? defaults to all 3"),
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
                          "for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput("cohortData2001", "data.table", sourceURL = NA,
                 paste0("Table that defines the cohorts by pixelGroup in 2001")),
    expectsInput("cohortData2011", "data.table", sourceURL = NA,
                 paste0("Table that defines the cohorts by pixelGroup in 2011")),
    expectsInput("spreadFirePoints", "list", sourceURL = NA,
                 paste("named list of `SpatialPointsDataFrame`s for each fire year",
                       "with each point denoting an ignition location.")),
    expectsInput("firePolys", "list", sourceURL = NA,
                 paste0("List of `SpatialPolygonsDataFrame`s representing annual fire polygons.",
                        "List must be named with followign convention: 'year<numeric year>'")),
    expectsInput("firePolysForAge", "list", sourceURL = NA,
                 "firePolys used to classify timeSinceDisturbance in nonforest LCC"),
    expectsInput("historicalFireRaster", "RasterLayer",
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip",
                 "a raster with values representing fire year 1985-2020"),
    expectsInput("flammableRTM", "RasterLayer", sourceURL = NA,
                 "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
    expectsInput("historicalClimateRasters", "list", sourceURL = NA,
                 paste("length-one list of containing a raster stack of historical climate",
                       "list named after the variable and raster layers named as 'year<numeric year>'")),
    expectsInput("ignitionFirePoints", "list", sourceURL = NA,
                 paste("list of spatialPolygonDataFrame objects representing annual ignition locations.",
                       "This includes all fires regardless of size")),
    expectsInput("nonForestedLCCGroups", "list",
                 paste("a named list of non-forested landcover groups",
                       "e.g. list('wetland' = c(19, 23, 32))",
                       "These will become covariates in fireSense_IgnitionFit")),
    expectsInput("pixelGroupMap2001", "RasterLayer", sourceURL = NA,
                 "RasterLayer that defines the pixelGroups for cohortData table in 2001"),
    expectsInput("pixelGroupMap2011", "RasterLayer",
                 "RasterLayer that defines the pixelGroups for cohortData table in 2011"),
    expectsInput("rasterToMatch", "RasterLayer", sourceURL = NA,
                 "template raster for study area. Assumes some buffering of core area to limit edge effect of fire."),
    expectsInput("rstLCC", "RasterLayer", sourceURL = NA,
                 "Raster of land cover. Defaults to LCC05."),
    expectsInput("sppEquiv", "data.table", sourceURL = NA,
                 "table of LandR species equivalencies"),
    expectsInput("standAgeMap2001", "RasterLayer", sourceURL = NA,
                 "map of stand age in 2001 used to create cohortData2001"),
    expectsInput("standAgeMap2011", "RasterLayer", sourceURL = NA,
                 "map of stand age in 2011 used to create cohortData2011"),
    expectsInput("studyArea", "SpatialPolygonsDataFrame", sourceURL = NA,
                 "studyArea that determines spatial boundaries of all data")
  ),
  outputObjects = bindrows(
    createsOutput("fireBufferedListDT", "list",
                  "list of data.tables with fire id, pixelID, and buffer status"),
    createsOutput("firePolys", "list",
                  "list of spatialPolygonDataFrame objects representing annual fires"),
    createsOutput("fireSense_annualSpreadFitCovariates", "list",
                  "list of tables with climate covariates, youngAge, burn status, polyID, and pixelID"),
    createsOutput("fireSense_escapeCovariates", "data.table",
                  "ignition covariates with added column of escapes"),
    createsOutput("fireSense_escapeFormula", "character",
                  "formula for escape, using fuel classes and landcover, as character"),
    createsOutput("fireSense_ignitionCovariates", "data.table",
                  "table of aggregated ignition covariates with annual ignitions"),
    createsOutput("fireSense_ignitionFormula", "character",
                  "formula for ignition, using climate and vegetation covariates, as character"),
    createsOutput("fireSense_nonAnnualSpreadFitCovariates", "list",
                  "list of two tables with veg covariates, burn status, polyID, and pixelID"),
    createsOutput("fireSense_spreadFormula", "character",
                  "formula for spread, using climate and vegetation covariates, as character"),
    createsOutput("ignitionFitRTM", "RasterLayer",
                  paste("A (template) raster with information with regards to the spatial",
                        "resolution and geographical extent of `fireSense_ignitionCovariates`.",
                        "Used to pass this information onto `fireSense_ignitionFitted`",
                        "Needs to have number of non-NA cells as attribute (`ignitionFitRTM@data@attributes$nonNAs`).")),
    createsOutput("landcoverDT", "data.table",
                  paste("data.table with `pixelID` and relevant landcover classes",
                        "that is used by predict functions.")),
    createsOutput("nonForest_timeSinceDisturbance2001", "RasterLayer",
                  "time since burn for non-forested pixels in 2001"),
    createsOutput("nonForest_timeSinceDisturbance2011", "RasterLayer",
                  "time since burn for non-forested pixels in 2011"),
    createsOutput("spreadFirePoints", "list",
                  paste("Named list of `SpatialPolygonDataFrame` objects representing annual fire centroids.",
                        "This only includes fires that escaped (e.g. `size > res(flammableRTM)`.")),
    createsOutput("terrainDT", "data.table",
                  "`data.table` with `pixelID` and relevant terrain variables used by predict models.")
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
      if ("fireSense_SpreadFit" %in% P(sim)$whichModulesToPrepare) {
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepSpreadFitData")
      }
      sim <- scheduleEvent(sim, end(sim), "fireSense_dataPrepFit", "plotAndMessage", eventPriority = 9)
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
  ## TODO: correct this if ignitionFuelClass and spreadFuelClass are used
  igFuels <- sim$sppEquiv[[P(sim)$ignitionFuelClassCol]]
  spreadFuels <- sim$sppEquiv[[P(sim)$spreadFuelClassCol]]

  if (any(c(is.null(spreadFuels), is.na(spreadFuels),
            is.null(igFuels), is.na(igFuels)))) {
    stop("All species must have spread and ignition fuelClasses defined")
  }

  sim$landcoverDT <- makeLandcoverDT(rstLCC = sim$rstLCC,
                                     flammableRTM = sim$flammableRTM,
                                     forestedLCC = P(sim)$forestedLCC,
                                     nonForestedLCCGroups = sim$nonForestedLCCGroups)

  # cannot merge because before subsetting due to column differences over time
  #TODO: firePolysForAge should be lists of sf... can't reliably read NFDB with sp
  firePolysForAge <- lapply(sim$firePolysForAge[lengths(sim$firePolysForAge) > 0],
                            FUN = st_as_sf) %>%
    lapply(., FUN = function(x){
      x <- x[, "YEAR"]
    }) %>%
    do.call(rbind, .)

  #TODO: spreadFit should not use this object if P(sim)$nonForestCanBeYoungAge is TRUE
  if (P(sim)$nonForestCanBeYoungAge){
    #if fireRaster is null, it will use firePolys
    sim$nonForest_timeSinceDisturbance2001 <- makeTSD(year = 2001,
                                                      fireRaster = sim$historicalFireRaster,
                                                      firePolys = sim$firePolysForAge,
                                                      standAgeMap = sim$standAgeMap2001, lcc = sim$landcoverDT,
                                                      cutoffForYoungAge = P(sim)$cutoffForYoungAge)
    sim$nonForest_timeSinceDisturbance2011 <- makeTSD(year = 2011, fireRaster = sim$historicalFireRaster,
                                                      firePolys = sim$firePolysForAge,
                                                      standAgeMap = sim$standAgeMap2011,
                                                      lcc = sim$landcoverDT,
                                                      cutoffForYoungAge = P(sim)$cutoffForYoungAge)
  }
  origDTThreads <- data.table::setDTthreads(2)

  #TODO: untill we standardize the youngAge treatment there is no point in this
  #standardizing meaning calculating youngAge over the whole landscape every year.
  #currently we do this in spread but only for fire buffers, not whole landscape.

  # cohorts2001 <-  Cache(castCohortData,
  #                       cohortData = sim$cohortData2001,
  #                       pixelGroupMap = sim$pixelGroupMap2001,
  #                       year = 2001,
  #                       ageMap = NULL,
  #                       cutoffForYoungAge = P(sim)$cutoffForYoungAge,
  #                       lcc = sim$landcoverDT,
  #                       missingLCC = P(sim)$missingLCC)
  # cohorts2011 <-  Cache(castCohortData, cohortData = sim$cohortData2011,
  #                       pixelGroupMap = sim$pixelGroupMap2011,
  #                       year = 2011,
  #                       cutoffForYoungAge = P(sim)$cutoffForYoungAge,
  #                       ageMap = NULL,
  #                       lcc = sim$landcoverDT,
  #                       missingLCC = P(sim)$missingLCC)
  #
  # vegData <- rbindlist(list(cohorts2001, cohorts2011), use.names = TRUE)

  flammableIndex <- data.table(index = 1:ncell(sim$flammableRTM), value = getValues(sim$flammableRTM)) %>%
    .[value == 1,] %>%
    .$index
  mod$climateDT <- Cache(climateRasterToDataTable,
                         historicalClimateRasters = sim$historicalClimateRasters,
                         Index = flammableIndex, userTags = c("climateRasterToDataTable"))
  rm(flammableIndex)

  #needed by prep spread
  return(invisible(sim))
}

prepare_SpreadFit <- function(sim) {
  ## Put in format for DEOptim that distinguishes annual and nonannual covariates
  ## Prepare annual spread fit covariates

  ####prep veg data####
  doAssertion <- getOption("fireSenseUtils.assertions", TRUE)

  ## sanity check the inputs
  compareRaster(sim$rasterToMatch, sim$flammableRTM, sim$rstLCC,
                sim$standAgeMap2001, sim$standAgeMap2011)
  lapply(sim$historicalClimateRasters, compareRaster, x = sim$rasterToMatch)


  #prep the fire data
  if (P(sim)$useRasterizedFire){
    sim <- prepare_SpreadFitFire_Raster(sim)
  } else {
    sim <- prepare_SpreadFitFire_Vector(sim)
  }

  ## output filenames ------------------------------------------------------------------------------
  mod$vegFile <- file.path(outputPath(sim),
                           paste0("fireSense_SpreadFit_veg_coeffs_", P(sim)$.studyAreaName, ".txt"))
  vegData <- Map(f = cohortsToFuelClasses,
                       cohortData = list(sim$cohortData2001, sim$cohortData2011),
                       yearCohort = list(2001, 2011),
                       pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                       MoreArgs = list(sppEquiv = sim$sppEquiv,
                                       sppEquivCol = P(sim)$sppEquivCol,
                                       flammableRTM = sim$flammableRTM,
                                       fuelClassCol = P(sim)$spreadFuelClassCol,
                                       cutoffForYoungAge = -1)) #youngAge will be calculated every year

  vegData <- lapply(vegData, FUN = function(x){
    dt <- as.data.table(getValues(x))
    dt[, pixelID := 1:ncell(x)]
    return(dt)
  })

  vegData[[1]][, year := 2001]
  vegData[[2]][, year := 2011]

  #joining with landcoverDT eliminates non-flammable pixels - every remaining pixel MUST have a value
  #join landcover separately, as the forest extent changes between 2001 and 2011
  vegData[[1]] <- vegData[[1]][sim$landcoverDT, on = c("pixelID")]
  vegData[[2]] <- vegData[[2]][sim$landcoverDT, on = c("pixelID")]
  vegData <- rbindlist(vegData)
  setcolorder(vegData, c("pixelID", "year"))

  #set 'orphaned' pixels as P(sim)$missingLCC - the forest that is not LandR forest
  lccNames <- setdiff(names(vegData), c("pixelID", "year"))
  vegData[, missingLCC := rowSums(vegData[, .SD, .SDcols = lccNames])]
  vegData[missingLCC == 0, eval(P(sim)$missingLCC) := 1]
  ## TODO: when we add assertions, assert that there are no rows where missingLCC = 2
  vegData[, missingLCC := NULL]


  ####join fire and veg data ####
  pre2011 <- paste0("year", min(P(sim)$fireYears):2010)
  pre2011Indices <- sim$fireBufferedListDT[names(sim$fireBufferedListDT) %in% pre2011] %>%
    rbindlist(.) %>%
    vegData[year < 2011][., on = c("pixelID")]
  ## TODO, review why is.na is included here? what would be NA? pixels in fireBufferedList that aren't in vegData?
  pre2011Indices[is.na(year), year := 2001]

  post2011Indices <- sim$fireBufferedListDT[!names(sim$fireBufferedListDT) %in% pre2011] %>%
    rbindlist(.) %>%
    vegData[year >= 2011][., on = c("pixelID")]
  post2011Indices[is.na(year), year := 2011]

  ## Some pixels will be NA because the polygon includes non-flammable cells
  ## As long as these pixels are also NA in climate data, no issue

  fireSenseVegData <- rbind(pre2011Indices, post2011Indices)
  setnames(fireSenseVegData, "buffer", "burned")

  vegCols <- setdiff(names(fireSenseVegData), c("pixelID", "burned", "ids", "year"))
  dropCols <- names(which(apply(fireSenseVegData[, ..vegCols], 2, sum) == 0))

  ## spreadFit will fail if there are empty (all zero) columns
  if (length(dropCols) > 0) {
    message("Dropping column(s) from spreadFit covariate table: ",
            paste(dropCols, collapse = ", "))
    vegCols <- vegCols[!vegCols %in% dropCols]
    set(fireSenseVegData, NULL, dropCols, NULL)
  }

  if (isTRUE(doAssertion)) {
    ttt <- table(fireSenseVegData$burned)
    ratioZeroToOne <- ttt[1]/ttt[2]
    if (ratioZeroToOne < 5)
      stop("The number of pixels in the fire buffers should be at least 5x the number of burned pixels\n",
           "Please create larger buffers around fires in fireBufferedListDT, e.g., via ",
           "fireSenseUtils::bufferToArea(..., areaMultiplier = multiplier)")
  }

  RHS <- paste(paste0(names(sim$historicalClimateRasters)), "youngAge",
                paste0(vegCols, collapse = " + "), sep =  " + ")

  ## this is a funny way to get years but avoids years with 0 fires
  years <- paste0("year", P(sim)$fireYears)

  yearsWithFire <- years[years %in% names(sim$fireBufferedListDT)]

  pre2011int <- as.integer(min(P(sim)$fireYears):2010)
  post2011int <- as.integer(2011:max(P(sim)$fireYears))

  pre2011 <- yearsWithFire[yearsWithFire %in% paste0("year", pre2011int)]
  post2011 <- yearsWithFire[yearsWithFire %in% paste0("year", post2011int)]

  fbl <- rbindlist(sim$fireBufferedListDT, idcol = "year")
  rmCols <- setdiff(colnames(fbl), c("pixelID", "year"))
  set(fbl, NULL, rmCols, NULL)
  fbl <- mod$climateDT[fbl, on = c("year", "pixelID"), nomatch = NULL]
  fireSense_annualSpreadFitCovariates <- split(fbl, by = "year", keep.by = FALSE)

  ## prepare non-annual spread fit covariates by getting the youngAge
  pre2011Indices <- sim$fireBufferedListDT[names(sim$fireBufferedListDT) %in% pre2011]
  post2011Indices <- sim$fireBufferedListDT[!names(sim$fireBufferedListDT) %in% pre2011]
  colsToExtract <- c("pixelID", vegCols)

  nonAnnualPre2011 <- fireSenseVegData[year < 2011, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID),]

  nonAnnualPost2011 <- fireSenseVegData[year >= 2011, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID)] ## remove duplicates from same pixel diff year

  ## the function will do this below, and then use the data.table with location of non-forest to fill in those ages
  ## pmap allows for internal debugging when there are large lists that are passed in; Map does not
  annualCovariates <- Cache(
    purrr::pmap,
    .l = list(
      #years = list(c(2001:2010), c(2011:max(P(sim)$fireYears))),
      years = list(pre2011int, post2011int),
      annualCovariates = list(fireSense_annualSpreadFitCovariates[pre2011],
                              fireSense_annualSpreadFitCovariates[post2011]),
      standAgeMap = list(sim$nonForest_timeSinceDisturbance2001,
                         sim$nonForest_timeSinceDisturbance2011)
    ),
    .f = calcYoungAge,
    fireBufferedListDT = sim$fireBufferedListDT,
    cutoffForYoungAge = P(sim)$cutoffForYoungAge
  )

  sim$fireSense_annualSpreadFitCovariates <- do.call(c, annualCovariates)

  sim$fireSense_nonAnnualSpreadFitCovariates <- list(nonAnnualPre2011, nonAnnualPost2011)
  names(sim$fireSense_nonAnnualSpreadFitCovariates) <- c(paste(names(pre2011Indices), collapse = "_"),
                                                         paste(names(post2011Indices), collapse = "_"))

  sim$fireSense_spreadFormula <- paste0("~ 0 + ", RHS)

  return(invisible(sim))
}

prepare_SpreadFitFire_Raster <- function(sim) {


  historicalFireRaster <- sim$historicalFireRaster


  #build initial burn IDs by buffering  - then using clump(raster) or patches(terra)
  if (inherits(historicalFireRaster, "RasterLayer")) {
    historicalFireRaster <- terra::rast(historicalFireRaster)
  } else {
    historicalFireRaster <- historicalFireRaster
  }

  compareGeom(historicalFireRaster, sim$flammableRTM)
  #historical fire Raster is currently not in outputs - if assigned to sim here, it should be added
  # as we modify it by removing non-flammable fires.
  #TODO: discuss the above
  historicalFireRaster[sim$flammableRTM == 0] <- NA

  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1, length(sim$fireYears))

  #this is analogous to buffer to area but for raster datasets as opposed to polygon
  #the inner looping function is very similar - one difference is that non-flammable
  #pixels do not count toward the buffer size, unlike the polygonal version.
  sim$fireBufferedListDT <- Cache(rasterFireBufferDT, years =  P(sim)$fireYears,
                                  fireRaster = historicalFireRaster, flammableRTM = sim$flammableRTM,
                                  bufferForFireRaster = P(sim)$bufferForFireRaster, verb = 1,
                                  areaMultiplier = P(sim)$areaMultiplier, minSize = P(sim)$minBufferSize,
                                  cores = nCores, userTags = c(currentModule(sim), "rasterFireBufferDT"))
  #TODO: test that this is the correct method for missing years
  missingYears <- unlist(lapply(sim$fireBufferedListDT, is.null))

  if (any(missingYears)){
    actualFireYears <- P(sim)$fireYears[!missingYears]
    sim$fireBufferedListDT <- sim$fireBufferedListDT[!missingYears]
  }

  #next up: generate spread fire points. no harmonization is needed with this approach :)
  sim$spreadFirePoints <- rasterFireSpreadPoints(fireBufferDT = sim$fireBufferedListDT[1], flammableRTM = sim$flammableRTM)

  #this is temporary while we migrate out of spatial/raster constructs
  sim$spreadFirePoints <- lappy(sim$spreadFirePoints, as_Spatial)
  names(sim$spreadFirePoints) <- names(sim$fireBufferedListDT)

  return(invisible(sim))
}

prepare_SpreadFitFire_Vector <- function(sim){
  #sanity check
  stopifnot(
    "all annual firePolys are not within studyArea" = all(unlist(lapply(sim$firePolys, function(x) {
      nrow(st_as_sf(x)) == nrow(sf::st_intersection(st_as_sf(x), st_as_sf(sim$studyArea)))
    })))
  )

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
                               ## ensure matching IDs
                               points <- points[points$FIRE_ID %in% polys$FIRE_ID,]
                               polys <- polys[polys$FIRE_ID %in% points$FIRE_ID,]
                               points$FIRE_ID <- as.numeric(as.factor(points$FIRE_ID))
                               polys$FIRE_ID <- as.numeric(as.factor(polys$FIRE_ID))
                               return(list(polys = polys, points = points))
                             })
    sim$spreadFirePoints <- lapply(PointsAndPolys, FUN = function(x) x[["points"]])
    sim$firePolys <- lapply(PointsAndPolys, FUN = function(x) x[["polys"]])
    rm(PointsAndPolys)
    names(sim$firePolys) <- origNames
    names(sim$spreadFirePoints) <- origNames
  }

  ## drop fires less than 1 px in size
  pixSizeHa <- prod(res(sim$flammableRTM)) / 1e4
  sim$spreadFirePoints <- lapply(sim$spreadFirePoints, function(x) {
    x <- subset(x, SIZE_HA > pixSizeHa)
    if (nrow(x) > 0) x else NULL
  })
  sim$spreadFirePoints[sapply(sim$spreadFirePoints, is.null)] <- NULL

  sim$firePolys <- lapply(sim$firePolys, function(x) {
    x <- subset(x, SIZE_HA > pixSizeHa)
    if (nrow(x) > 0) x else NULL
  })
  sim$firePolys[sapply(sim$firePolys, is.null)] <- NULL

  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1, length(sim$firePolys))
  fireBufferedListDT <- Cache(bufferToArea,
                              poly = sim$firePolys,
                              polyName = names(sim$firePolys),
                              rasterToMatch = sim$flammableRTM, ## TODO: use sim$rasterToMatch here?
                              verb = TRUE,
                              areaMultiplier = P(sim)$areaMultiplier,
                              field = "FIRE_ID",
                              cores = nCores,
                              minSize = P(sim)$minBufferSize,
                              userTags = c("bufferToArea", P(sim)$.studyAreaName),
                              omitArgs = "cores")

  ## drop fire years from these lists that don't have any buffer points pre-harmonization
  omitYears <- sapply(fireBufferedListDT, function(x) nrow(x) == 0)
  fireBufferedListDT[omitYears] <- NULL
  sim$firePolys[omitYears] <- NULL
  sim$spreadFirePoints[omitYears] <- NULL

  ## Clean up missing pixels - this is a temporary fix
  ## we will always have NAs because of edge pixels - will be an issue when predicting
  ## The next 1 line replaces the 8 or so lines after
  fireBufferedListDT <- Cache(rmMissingPixels, fireBufferedListDT, vegData$pixelID)

  ## Post buffering, new issues --> must make sure points and buffers match
  pointsIDColumn <- "FIRE_ID"
  sim$spreadFirePoints <- Cache(harmonizeBufferAndPoints,
                                cent = sim$spreadFirePoints,
                                buff = fireBufferedListDT,
                                ras = sim$flammableRTM,
                                idCol = pointsIDColumn,
                                userTags = c("harmonizeBufferAndPoints", P(sim)$.studyAreaName))

  ## drop fire years from these lists that don't have any buffer points post-harmonization
  omitYears <- sapply(sim$spreadFirePoints, is.null)
  fireBufferedListDT[omitYears] <- NULL
  sim$firePolys[omitYears] <- NULL
  sim$spreadFirePoints[omitYears] <- NULL

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
      if (any(flammableInPolys, na.rm = TRUE)) {
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
  sim$fireBufferedListDT <- out22$FireBuffered

  return(invisible(sim))
}

prepare_IgnitionFit <- function(sim) {


  stopifnot(
    "all ignitionFirePoints are not within studyArea" = identical(
      nrow(st_as_sf(sim$ignitionFirePoints)),
      nrow(st_intersection(st_as_sf(sim$ignitionFirePoints), st_as_sf(sim$studyArea)))
    ))

  ## account for forested pixels that aren't in cohortData
  sim$landcoverDT[, rowSums := rowSums(.SD), .SD = setdiff(names(sim$landcoverDT), "pixelID")]
  forestPix <- sim$landcoverDT[rowSums == 0,]$pixelID
  problemPix2001 <- forestPix[is.na(sim$pixelGroupMap2001[forestPix])]
  problemPix2011 <- forestPix[is.na(sim$pixelGroupMap2011[forestPix])]
  set(sim$landcoverDT, NULL, 'rowSums', NULL)

  ## The non-forests aren't the same between years, due to cohortData being different
  landcoverDT2001 <- copy(sim$landcoverDT)
  landcoverDT2001[pixelID %in% problemPix2001, eval(P(sim)$missingLCC) := 1]
  landcoverDT2011 <- copy(sim$landcoverDT)
  landcoverDT2011[pixelID %in% problemPix2011, eval(P(sim)$missingLCC) := 1]
  ## first put landcover into raster stack
  ## non-flammable pixels require zero values for non-forest landcover, not NA
  LCCras <- Cache(Map,
                  f = putBackIntoRaster,
                  landcoverDT = list(landcoverDT2001,landcoverDT2011),
                  MoreArgs = list(lcc = names(sim$nonForestedLCCGroups),
                                  flammableMap = sim$flammableRTM),
                  userTags = c("putBackIntoRaster", P(sim)$.studyAreaName))  %>%
    lapply(., FUN = brick)

  fuelClasses <- Map(f = cohortsToFuelClasses,
                     cohortData = list(sim$cohortData2001, sim$cohortData2011),
                     yearCohort = list(2001, 2011),
                     pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                     MoreArgs = list(sppEquiv = sim$sppEquiv,
                                     sppEquivCol = P(sim)$sppEquivCol,
                                     landcoverDT = sim$landcoverDT,
                                     flammableRTM = sim$flammableRTM,
                                     cutoffForYoungAge = P(sim)$cutoffForYoungAge))

  if (class(fuelClasses[[1]]) == "RasterStack") {
    ## if it is already a brick, running this line will delete the values
    fuelClasses <- lapply(fuelClasses, FUN = raster::brick)
  }

  if (P(sim)$nonForestCanBeYoungAge) {
    ## this modifies the NF landcover by converting some NF to a new YA layer
    ## it must be done before aggregating
    LCCras <- Map(f = calcNonForestYoungAge,
                     landcoverDT = list(landcoverDT2001, landcoverDT2011),
                     NFTSD = list(sim$nonForest_timeSinceDisturbance2001,
                                  sim$nonForest_timeSinceDisturbance2011),
                     LCCras = list(LCCras[[1]], LCCras[[2]]),
                     MoreArgs = list(cutoffForYoungAge = P(sim)$cutoffForYoungAge))

    for (i in c(1:2)) {
      if ("youngAge" %in% names(fuelClasses[[i]])) {
        #brick1[x] + brick2[x] stalled..
        YA1 <- fuelClasses[[i]]$youngAge
        YA2 <- LCCras[[i]]$youngAge
        bothYA <- YA1 + YA2
        fuelClasses[[i]]$youngAge <- bothYA
      }  else {
        fuelClasses[[i]]$youngAge <- LCCras[[i]]$youngAge
      }
      toKeep <- setdiff(names(LCCras[[i]]), "youngAge")
      LCCras[[i]] <- raster::subset(LCCras[[i]], toKeep) ## to avoid double-counting
    }
  }

  LCCras <- lapply(LCCras, FUN = aggregate, fact = P(sim)$igAggFactor, fun = mean)
  names(LCCras) <- c("year2001", "year2011")
  fuelClasses <- lapply(fuelClasses, FUN = aggregate, fact = P(sim)$igAggFactor, fun = mean)
  names(fuelClasses) <- c("year2001", "year2011")

  climate <- sim$historicalClimateRasters
  climVar <- names(climate)
  if (length(climate) > 1) {
    stop("need to fix ignition for multiple climate variables. contact module developers") ## TODO!
  } else {
    climate <- raster::stack(climate[[1]]) %>%
      aggregate(., fact = P(sim)$igAggFactor, fun = mean)
  }
  ## ignition won't have same years as spread so we do not use names of init objects
  ## The reason is some years may have no significant fires, e.g. 2001 in RIA
  pre2011 <- paste0("year", min(P(sim)$fireYears):2010)
  post2011 <- paste0("year", 2011:max(P(sim)$fireYears))

  #this is joining fuel class, LCC, and climate, subsetting to flamIndex, calculating n of ignitions
  fireSense_ignitionCovariates <- Map(f = fireSenseUtils::stackAndExtract,
                                      years = list(pre2011, post2011),
                                      fuel = list(fuelClasses$year2001, fuelClasses$year2011),
                                      LCC = list(LCCras$year2001, LCCras$year2011),
                                      MoreArgs = list(climate = climate,
                                                      fires = sim$ignitionFirePoints,
                                                      climVar = climVar #TODO: this is clunky, rethink
                                      ))

  fireSense_ignitionCovariates <- rbindlist(fireSense_ignitionCovariates)

  #remove any pixels that are 0 for all classes
  fireSense_ignitionCovariates[, coverSums := rowSums(.SD), .SD = setdiff(names(fireSense_ignitionCovariates),
                                                                          c(climVar, "cells", "ignitions", "year"))]
  fireSense_ignitionCovariates <- fireSense_ignitionCovariates[coverSums > 0]
  if (any(fireSense_ignitionCovariates$coverSums > 1)) {
    stop("error with ignition raster aggregation")
  }
  set(fireSense_ignitionCovariates, NULL, "coverSums", NULL)

  #rename cells to pixelID - though aggregated raster is not saved
  setnames(fireSense_ignitionCovariates, old = "cells", new = "pixelID")
  fireSense_ignitionCovariates[, year := as.numeric(year)]
  firstCols <- c("pixelID", "ignitions", climVar, "youngAge")
  firstCols <- firstCols[firstCols %in% names(fireSense_ignitionCovariates)]
  setcolorder(fireSense_ignitionCovariates, neworder = firstCols)
  sim$fireSense_ignitionCovariates <- fireSense_ignitionCovariates


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
  mod$firePolysForAge <- NULL
  mod$fireSenseVegData <- NULL
  mod$climateDT <- NULL

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  sim <- saveFiles(sim)
}

### template for plot events
plotAndMessage <- function(sim) {
  #TODO: this could plot the ignition/spread covariates
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), P(sim)$.studyAreaName)
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("studyArea", sim)) {
    stop("Please supply study area - this object is key")
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- LandR::prepInputsLCC(year = 2005, ## TODO: use 2010
                                              destinationPath = dPath,
                                              studyArea = sim$studyArea,
                                              useCache = TRUE)
  }

  if (!all(suppliedElsewhere("cohortData2011", sim),
           suppliedElsewhere("pixelGroupMap2011", sim),
           suppliedElsewhere("pixelGroupMap2001", sim),
           suppliedElsewhere("cohortData2001", sim))) {
    stop("Stop - need cohortData and pixelGroupMap objects - contact module creators")
  }

  if (!P(sim)$useRasterizedFireForSpread) {
    if (!suppliedElsewhere("firePolys", sim) | !suppliedElsewhere("firePolysForAge", sim)) {
      # don't want to needlessly postProcess the same firePolys objects

      allFirePolys <- Cache(fireSenseUtils::getFirePolygons,
                            years = c(min(P(sim)$fireYears - P(sim)$cutoffForYoungAge):max(P(sim)$fireYears)),
                            studyArea = sim$studyArea,
                            destinationPath = dPath,
                            useInnerCache = TRUE,
                            userTags = c(cacheTags, "firePolys", paste0("years:", range(P(sim)$fireYears))))
    }

    if (!suppliedElsewhere("firePolys", sim)) {
      sim$firePolys <- allFirePolys[names(allFirePolys) %in% paste0("year", P(sim)$fireYears)]
    }

    if (!suppliedElsewhere("firePolysForAge", sim)) {
      sim$firePolysForAge <- allFirePolys
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

      mc <- pemisc::optimalClusterNum(2e3, maxNumClusters = length(sim$firePolys))
      clObj <- parallel::makeCluster(type = "SOCK", mc)
      a <- parallel::clusterEvalQ(cl = clObj, {library(raster); library(rgeos)})
      clusterExport(cl = clObj, list("firePolys"), envir = sim)
      sim$spreadFirePoints <- Cache(FUN = parallel::clusterApply,
                                    x = sim$firePolys,
                                    cl = clObj,
                                    fun = centerFun, #don't specify FUN argument or Cache will mistake it.
                                    userTags = c(cacheTags, "spreadFirePoints"),
                                    omitArgs = c("userTags", "mc.cores", "useCloud", "cloudFolderID"))
      stopCluster(clObj)
      names(sim$spreadFirePoints) <- names(sim$firePolys)
    }

    if (all(!is.null(sim$spreadFirePoints), !is.null(sim$firePolys))) {
      ## may be NULL if passed by objects - add to Init?
      ## this is necessary because centroids may be fewer than fires if fire polys were small
      min1Fire <- lapply(sim$spreadFirePoints, length) > 0
      sim$spreadFirePoints <- sim$spreadFirePoints[min1Fire]
      sim$firePolys <- sim$firePolys[min1Fire]
    }

    if (length(sim$firePolys) != length(sim$spreadFirePoints)) {
      stop("mismatched years between firePolys and firePoints")
      ## TODO: need to implement a better approach that matches each year's IDS
      ## these are mostly edge cases if a user passes only one of spreadFirePoints/firePolys
    }
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

  if (!suppliedElsewhere("flammableRTM", sim)) {
    sim$flammableRTM <- defineFlammable(sim$rstLCC,
                                        nonFlammClasses = P(sim)$nonflammableLCC,
                                        mask = sim$rasterToMatch,
                                        filename2 = file.path(dPath,paste0("flammableRTM_",
                                                                           P(sim)$.studyAreaName,
                                                                           ".tif")))
  }

  if (!suppliedElsewhere("historicalFireRaster", sim)) {
    sim$historicalFireRaster <- Cache(prepInputs,
                                      url = extractURL("historicalFireRaster", sim),
                                      rasterToMatch = sim$flammableRTM,
                                      destinationPath = dPath,
                                      studyArea = sim$studyArea,
                                      method = "ngb",
                                      filename2 = paste0("wildfire_", P(sim)$.studyAreaName, ".tif"))
    #make sure this is near or ngb
  }

  if (!suppliedElsewhere("nonForestedLCCGroups", sim)) {
    sim$nonForestedLCCGroups <- list(
      "nonForest_highFlam" = c(8, 10, 14), #shrubland, grassland, wetland
      "nonForest_lowFlam" = c(11, 12, 15)) #shrub-lichen-moss + cropland. 2 barren classes are non-flammable
  }

  return(invisible(sim))
}

rmMissingPixels <- function(fbldt, pixelIDsAllowed)  {
  fbldt <- rbindlist(fbldt, idcol = "year")
  fbldt <- fbldt[pixelID %in% unique(pixelIDsAllowed)]
  fireBufferedListDT <- split(fbldt, by = "year", keep.by = FALSE)
}
