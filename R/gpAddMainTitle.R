#' gpAddMainTitle
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
#' @export gpAddMainTitle
#'




gpAddMainTitle <- function(genoplot, main, ...) {
  #TODO: Check parameters
  karyoploteR::kpAddMainTitle(gp$global.kp, main = main, ...)
}
