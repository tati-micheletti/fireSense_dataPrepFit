calcNonForestYoungAge <- function(landcoverDT, NFTSD, LCCras, cutoffForYoungAge){
  nfLCC <- setdiff(colnames(landcoverDT), "pixelID")
  landcoverDT[, sumRows := rowSums(.SD), .SDcol = nfLCC]
  landcoverDT[, age := NFTSD[pixelID]]
  #this need to be chagned in LCCras and also converted to a youngAge raster
  #as the rasters will be aggregated
  pixToChange <- landcoverDT[sumRows > 0 & age < cutoffForYoungAge]$pixelID
  landcoverDT[, c("sumRows", "age") := NULL]
  LCCras[pixToChange] <- 0
  youngAge <- raster(LCCras)
  youngAge[!is.na(LCCras[[1]][])] <- 0 #NA in the landscape must be 0
  youngAge[pixToChange] <- 1
  LCCras$youngAge <- youngAge
  return(LCCras)
}
