#'getDefaultPlotParams
#'
#' @note: TODO documentation
#'
#' @export getDefaultPlotParams

getDefaultPlotParams <- function(show.introns=FALSE, compress.introns=TRUE, proportional.introns=TRUE) {

  plot.type <- selectPlotType(show.introns=show.introns, compress.introns=compress.introns, proportional.introns=proportional.introns)

  if(plot.type=="only.exons") {
    pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
    #TODO: Modify any default values to adjust to our needs

    #Add genoploteR specific values
    pp$outer.margin.bases <- 10
    pp$inner.margin.bases <- 5
    pp$margin.between.regions <- 0.02

    return(pp)
  } else if(plot.type=="compressed.introns.proportional") {
    pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
    #TODO: Modify any default values to adjust to our needs

    #Add genoploteR specific values
    pp$outer.margin.bases <- 10
    pp$inner.margin.bases <- 5
    pp$margin.between.regions <- 0

    return(pp)
  } else if(plot.type=="compressed.introns.all.equal") {

    pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
    #TODO: Modify any default values to adjust to our needs

    #Add genoploteR specific values
    pp$outer.margin.bases <- 10
    pp$inner.margin.bases <- 5
    pp$margin.between.regions <- 0
    return(pp)
  } else if(plot.type=="all") {

    pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
    #TODO: Modify any default values to adjust to our needs

    #Add genoploteR specific values
    pp$outer.margin.bases <- 10
    pp$inner.margin.bases <- 5
    pp$margin.between.regions <- 0
    return(pp)
  }
  stop("getDefaultPlotParams: plot type not found: ", plot.type)
}

