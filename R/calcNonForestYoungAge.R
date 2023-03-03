calcNonForestYoungAge <- function(landcoverDT, NFTSD, LCCras, cutoffForYoungAge){

  nfLCC <- setdiff(colnames(landcoverDT), "pixelID")
  landcoverDT[, sumRows := rowSums(.SD), .SDcol = nfLCC]
  landcoverDT[, age := NFTSD[pixelID]]
  #this need to be chagned in LCCras and also converted to a youngAge raster
  #as the rasters will be aggregated
  pixToChange <- landcoverDT[sumRows > 0 & age < cutoffForYoungAge]$pixelID

  landcoverDT[, c("sumRows", "age") := NULL]

  #check how terra works - this should change all pixels if length 2+ spat raster
  LCCras[pixToChange] <- 0
  youngAge <- rast(LCCras, nlyr = 1)

  temp <- values(LCCras[[1]])
  temp[!is.na(temp)] <- 0
  temp[pixToChange] <- 1

  youngAge <- setValues(youngAge, temp)

  LCCras$youngAge <- youngAge

  return(LCCras)
}
