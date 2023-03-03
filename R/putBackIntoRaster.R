putBackIntoRaster <- function(lcc, landcoverDT, flammableMap) {
  lccRasters <- list()
  for (i in 1:length(lcc)) {
    flamMapVals <- values(flammableMap, mat = FALSE)
    flamMapVals[!is.na(flamMapVals)] <- 0
    lccRasters[[i]] <- rast(flammableMap)
    lccVals <- landcoverDT[, get(lcc[i])]
    flamMapVals[landcoverDT$pixelID] <- lccVals
    lccRasters[[i]] <- setValues(lccRasters[[i]], flamMapVals)
  }
  lccRasters <- rast(lccRasters)
  names(lccRasters) <- lcc

  return(lccRasters)
}
