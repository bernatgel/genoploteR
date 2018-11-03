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




gpAddGeneStructure <- function(genoplot, coding.height=1, non.coding.height=0.5, introns.lwd=1, coding.exons.col="black", coding.exons.border="black", non.coding.exons.col="gray25", non.coding.exons.border="gray25", instrons.col="black", ...) {
  #TODO: Check parameters



#Assume there's a "type" mcol in exons with coding non coding
  coding.exons <- gp$exons[gp$exons$type=="coding"]
  non.coding.exons <- gp$exons[gp$exons$type=="non.coding"]
  introns <- gp$introns
  #TODO: enforce type in plotGene or deal with it here with some default
  for(i in seq_len(length(gp$regions.kp))) {
    karyoploteR::kpRect(gp$regions.kp[[i]], data=coding.exons, y0=(1-coding.height)/2, y1=1-(1-coding.height)/2, ymin=0, ymax=1, data.panel="ideogram", col=coding.exons.col, border=coding.exons.border, ...)
    karyoploteR::kpRect(gp$regions.kp[[i]], data=non.coding.exons, y0=(1-non.coding.height)/2, y1=1-(1-non.coding.height)/2, ymin=0, ymax=1, data.panel="ideogram", col=non.coding.exons.col, border=non.coding.exons.border, ...)
    karyoploteR::kpSegments(gp$regions.kp[[i]], data=introns, y0=0.5, y1=0.5, ymin=0, ymax=1, data.panel="ideogram", lwd=introns.lwd, col=instrons.col, ...)
  }
}
