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
