stackAndExtract <- function(years, climate, fuel, LCC, fires, climVar) {
  ignitionYears <- lapply(years, FUN = function(year){
    fires <- fires[paste0("year", fires$YEAR) %in% year,]
    climate <- climate[[year]]
    names(climate) <- climVar
    dtNames <- c(names(climate), names(fuel), names(LCC))
    yearCovariates <- raster::brick(climate, fuel, LCC)
    fireDT <- raster::extract(x = yearCovariates, y = fires, cellnumbers = TRUE) %>%
      as.data.table(.)
    fireDT <- fireDT[, nFires := .N, .(cells)]
    fireDT <- unique(fireDT)
    noFireDT <- raster::getValues(yearCovariates) %>%
      as.data.table(.) %>%
      .[, cells := 1:ncell(yearCovariates)] %>%
      na.omit(.)
    ignitionYear <- fireDT[noFireDT, on = c('cells', dtNames)]
    ignitionYear[is.na(nFires), nFires := 0]
    return(ignitionYear)
  })
 ignitionYears <- rbindlist(ignitionYears)
 return(ignitionYears)
}
