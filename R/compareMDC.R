library(Require)
Require(c("raster", "data.table", "ggplot2"))

#' Download and prepare fire data from National Fire Database
#'
#' @param historicalMDC raster stack of historical MDC
#' @param projectedMDC raster stack of projected MDC
#' @param flammableRTM an optional raster of flammable pixels to subset data
#' @return a ggplot object
#'
#' @export
#' @importFrom data.table data.table rbindlist
#' @importFrom ggplot2 aes ggplot geom_smooth geom_line
#' @importFrom stats median
#' @importFrom raster ncell nlayers getValues
compareMDC <- function(historicalMDC, projectedMDC, flammableRTM = NULL) {

  valfun <- function(x, flamMap = NULL){

    years <- 1:nlayers(x)
    if (!is.null(flamMap)){
      isFlam <- getValues(flamMap)
      flamMap <- c(1:ncell(flamMap))[!is.na(isFlam) & isFlam > 0]
    }
    medVal <- lapply(years, FUN = function(year, index = flamMap) {
      dt <- data.table("year" = year,
                       "MDC" = getValues(x[[year]]))
      if (!is.null(index)){
        dt <- dt[index]
      }
      dt <- dt[, .(MDC = median(MDC, na.rm = TRUE)), .(year)]
      return(dt)
    })
    medVal <- rbindlist(medVal)
    return(medVal)
  }

  proj <- valfun(projectedMDC, flamMap = flammableRTM)
  proj$year <- proj$year + 2010
  proj$year
  proj$stat <- "projected"
  hist <- valfun(historicalMDC, flamMap = flammableRTM)
  hist$year <- hist$year + 1989
  hist$stat <- "historical"
  MDC <- rbind(hist, proj)
  ggplot(data = MDC, aes(y = MDC, x = year, color = stat)) + geom_line() + geom_smooth()
}

compareMDC(historicalMDC = simOutPreamble$historicalClimateRasters$MDC,
           projectedMDC = simOutPreamble$projectedClimateRasters$MDC,
           flammableRTM = fSsimDataPrep$flammableRTM)
