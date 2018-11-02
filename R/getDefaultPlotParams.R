#'getDefaultPlotParams
#'
#' @note: TODO documentation
#'
#' @export getDefaultPlotParams

getDefaultPlotParams <- function(plot.type) {
  #TODO: check parameters are valid

  if(plot.type=="only.exons") {
    pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
    #TODO: Modify any default values to adjust to our needs

    #Add genoploteR specific values
    pp$outer.margin.bases <- 10
    pp$inner.margin.bases <- 5
    pp$margin.between.regions <- 0.02

    return(pp)
  }
  stop("getDefaultPlotParams: plot type not found: ", plot.type)
}

