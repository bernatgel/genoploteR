#' gpAddGeneLabel
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
#' @export gpAddGeneLabel
#'




gpAddGeneLabel <- function(genoplot, label, r0=NULL, r1=NULL, label.margin=0.01,  ...) {
  #TODO: Check parameters

  kpAddLabels(gp$global.kp, labels=label, r0=r0, r1=r1, label.margin = label.margin, data.panel="ideogram", pos=2, offset=0, ...)

}
