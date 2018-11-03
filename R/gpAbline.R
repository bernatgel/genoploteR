#' gpAbline
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
#' @export gpAbline
#'




gpAbline <- function(genoplot, chr=NULL, h=NULL, v=NULL, ymin=NULL, ymax=NULL, data.panel=1,
                     r0=NULL, r1=NULL, clipping=TRUE, ...) {
  #TODO: Check parameters

  for(i in seq_len(length(genoplot$regions.kp))) {
    kpAbline(genoplot$regions.kp[[i]], chr=chr, h=h, v=v, ymin=ymin, ymax=ymax, r0=r0, r1=r1, clipping=clipping, ...)
  }
}
