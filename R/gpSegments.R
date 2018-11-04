#' gpSegments
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
#' @export gpSegments
#'




gpSegments <- function(genoplot, data=NULL, chr=NULL, x0=NULL, x1=NULL, y0=NULL, y1=NULL,
                       ymin=NULL, ymax=NULL, data.panel=1, r0=NULL, r1=NULL,  clipping=TRUE,  ...) {
  #TODO: Check parameters

  for(i in seq_len(length(genoplot$regions.kp))) {
    kpSegments(genoplot$regions.kp[[i]], data=data, chr=chr, x0=x0, x1=x1, y0=y0, y1=y1,
               ymin=ymin, ymax=ymax, data.panel=1, r0=r0, r1=r1,  clipping=clipping,  ...)
  }
}
