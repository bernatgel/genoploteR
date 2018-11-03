#' gpAddBaseNumbers
#'
#' @description
#'
#' Plots the main title of the plot
#'
#' @details
#'
#' Given a GenoPlot object and a character string, plot the character strings
#' as the main title of the plot. This function is usually automatically
#' called by \code{\link{plotGene}}.
#'
#' @usage gpAddMainTitle(genoplot, main=NULL, ...)
#'
#' @param genoplot    a \code{GenoPlot} object returned by a call to \code{plotGene}
#' @param main (character) the main title of the plot
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#'
#' @return
#' invisibly returns the given GenoPlot object
#'
#' @seealso \code{\link{plotGene}}
#'
#' @examples
#'
#' gp <- plotGene(exons=c("chr1:100-200", "chr2:300-350"))
#' gpAddMainTitle(gp, main="EXAMPLE", col="red", srt=30)
#'
#' @export gpAddBaseNumbers
#'


#TODO: Reimplement here from scratch to:
# - Allow moving the ruler around
# - Getting at least one tick per exon (first base ideally?)
# - Possibility of setting the numbers to only to ticks per exon (first and last base)
# - Selecting between genomic and cDNA bases (genomic could actually use the genomic2cDNA function to compute?)
# - Decide whether to label number introns


gpAddBaseNumbers <- function(genoplot, tick.dist=500, tick.len=5, add.units=FALSE,
                             digits=2, minor.ticks=TRUE,
                             minor.tick.dist=100, minor.tick.len=2,  cex=0.5,
                             tick.col=NULL, minor.tick.col=NULL, clipping=TRUE,  ...) {
  for(i in seq_len(length(genoplot$regions.kp))) {
    kpAddBaseNumbers(genoplot$regions.kp[[i]], tick.dist=tick.dist, tick.len=tick.len, add.units=add.units,
                     digits=digits, minor.ticks=minor.ticks,
                     minor.tick.dist=minor.tick.dist, minor.tick.len=minor.tick.len,  cex=cex,
                     tick.col=tick.col, minor.tick.col=minor.tick.col, clipping=clipping,  ...)
  }
}
