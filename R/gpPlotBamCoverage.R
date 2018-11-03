#' gpAddGeneStructure
#'
#' @description
#'
#' Plots the exons and introns structure used to create the GenoPlot
#'
#' @details
#'
#' Given a GenoPlot object
#'
#' @usage
#'
#'
#' @return
#' invisibly returns the given GenoPlot object
#'
#' @seealso \code{\link{plotGene}}
#'
#' @examples
#'
# kp <- plotKaryotype(labels.plotter = NULL)
# gpAddMainTitle(kp, col="red", srt=30)
#'
#' @export gpAddGeneStructure
#'




gpPlotBAMCoverage <- function(genoplot, data=NULL, ymin=NULL, ymax=NULL,
                              data.panel=1, r0=NULL, r1=NULL,
                              col=NULL, border=NA, clipping=TRUE, ...) {
  #TODO: Check parameters

  #TODO: If ymax is NULL, set it to the max coverage of the WHOLE region!!!

  for(i in seq_len(length(genoplot$regions.kp))) {
    kpPlotBAMCoverage(genoplot$regions.kp[[i]], data=data, ymin=ymin, ymax=ymax, data.panel=data.panel, r0=r0, r1=r1, col=col, border=border, clipping=clipping, ...)
  }
}
