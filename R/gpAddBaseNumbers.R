#' gpAddBaseNumbers
#'
#' @description
#' Adds a ruler with base numbers to the plot
#'
#'
#' @details
#' This function adds a ruler to the plot marking the position of the plot elements in bases.
#' It can either label the ends of the exons or create a label every a certain number of bases
#' or both. The base numbers can be plotted any where in the genoplot as any other plot element
#' but it's possible to plot them just below or just above the gene structure setting
#' the parameter \code{position} to  \code{"below.gene"} (default) or \code{"above.gene"}.
#' The representation is quite customizable and it's possible to decide if ticks will be
#' plotted on intronic regions, if we want minor ticks between the labelled ticks,
#' the colors, sizes and distances for major and minor ticks, whether labels should have
#' units...
#'
#' With the parameter \code{col} it's possible to determine the color of all elements, but
#' there are parameters to specify the color of individual elements (labels, ticks...) and
#' these take precendence over \code{col}.
#'
#' @usage gpAddBaseNumbers(genoplot, label.exon.ends=TRUE, label.every.x=FALSE,
#'                            add.horizontal.line=FALSE, add.minor.ticks=TRUE,
#'                            custom.labels=NULL, ticks.in.introns=TRUE,
#'                            position="below.gene", inverted=FALSE,
#'                            add.units=TRUE, digits=2, shorten.labels=TRUE, cDNA=FALSE,
#'                            cex=0.6, label.col=NULL,
#'                            tick.dist=500, tick.len=0.2, tick.col=NULL,
#'                            minor.tick.dist=100, minor.tick.len=0.5, minor.tick.col=NULL,
#'                            r0=0, r1=1, data.panel=1,
#'                            col="black", pos=NULL, ...)
#'
#' @param genoplot (GEnoPlot object) The GenoPlot object returned by the call to \code{\link{plotGene}}
#' @param label.exon.ends (defaults to TRUE)
#' @param label.every.x (defaults to FALSE)
#' @param add.horizontal.line (defaults to FALSE)
#' @param add.minor.ticks (defaults to TRUE)
#' @param custom.labels (defaults to NULL)
#' @param ticks.in.introns (defaults to TRUE)
#' @param position (defaults to "below.gene")
#' @param inverted (defaults to FALSE)
#' @param add.units (defaults to TRUE)
#' @param digits (defaults to 2)
#' @param shorten.labels (defaults to TRUE)
#' @param cDNA (defaults to FALSE)
#' @param cex (defaults to 0.6)
#' @param label.col (defaults to NULL)
#' @param tick.dist (defaults to 500)
#' @param tick.len (defaults to 0.2)
#' @param tick.col (defaults to NULL)
#' @param minor.tick.dist (defaults to 100)
#' @param minor.tick.len (defaults to 0.5)
#' @param minor.tick.col (defaults to NULL)
#' @param r0 (defaults to 0)
#' @param r1 (defaults to 1)
#' @param data.panel (defaults to 1)
#' @param col (defaults to "black")
#' @param pos (defaults to NULL)
#' @param ...
#'
#'
#'
#' @return
#' invisibly returns the given GenoPlot object
#'
#' @seealso \code{\link{plotGene}}
#'
#' @examples
#'
#' gp <- plotGene(exons=toGRanges(c("chr1:1-200", "chr1:300-350")))
#' gpAddBaseNumbers(gp)
#'
#' @export gpAddBaseNumbers
#'
#' @importFrom GenomicRanges start end sort
#' @importFrom IRanges subsetByOverlaps
#' @importFrom assertthat assert_that is.number is.flag





gpAddBaseNumbers <- function(genoplot, label.exon.ends=TRUE, label.every.x=FALSE,
                             add.horizontal.line=FALSE, add.minor.ticks=TRUE,
                             custom.labels=NULL, ticks.in.introns=TRUE,
                             position="below.gene", inverted=FALSE,
                             add.units=TRUE, digits=2, shorten.labels=TRUE, cDNA=FALSE,
                             cex=0.6, label.col=NULL,
                             tick.dist=500, tick.len=0.2, tick.col=NULL,
                             minor.tick.dist=100, minor.tick.len=0.5, minor.tick.col=NULL,
                             r0=0, r1=1, data.panel=1,
                             col="black", pos=NULL, ...) {
    #Check the arguments
    if(!methods::is(genoplot, "GenoPlot")) stop("genoplot is not a GenoPlot object")
    assertthat::assert_that(assertthat::is.flag(label.exon.ends))
    assertthat::assert_that(assertthat::is.flag(label.every.x))
    assertthat::assert_that(assertthat::is.flag(add.horizontal.line))
    assertthat::assert_that(assertthat::is.flag(add.minor.ticks))
    assertthat::assert_that(is.null(custom.labels) || is.numeric(custom.labels))
    assertthat::assert_that(assertthat::is.flag(ticks.in.introns))
    position <- match.arg(position, c("custom", "above.gene", "below.gene"))
    assertthat::assert_that(assertthat::is.flag(inverted))
    assertthat::assert_that(assertthat::is.flag(add.units))
    assertthat::assert_that(assertthat::is.number(digits))
    assertthat::assert_that(assertthat::is.flag(shorten.labels))
    assertthat::assert_that(assertthat::is.flag(cDNA))
    assertthat::assert_that(assertthat::is.number(cex))
    assertthat::assert_that(is.null(label.col) || karyoploteR::is.color(label.col))
    assertthat::assert_that(assertthat::is.number(tick.dist))
    assertthat::assert_that(assertthat::is.number(tick.len))
    assertthat::assert_that(is.null(tick.col) || karyoploteR::is.color(tick.col))
    assertthat::assert_that(assertthat::is.number(minor.tick.dist))
    assertthat::assert_that(assertthat::is.number(minor.tick.len))
    assertthat::assert_that(is.null(minor.tick.col) || karyoploteR::is.color(minor.tick.col))
    assertthat::assert_that(assertthat::is.number(r0) && r0 >=0 && r0 <= 1)
    assertthat::assert_that(assertthat::is.number(r1) && r1 >=0 && r1 <= 1)
    #do not check data panel, since the valid values will usually vary depending on the current plot.type
    assertthat::assert_that(is.null(col) || karyoploteR::is.color(col))
    assertthat::assert_that(is.null(pos) || (assertthat::is.number(pos) && pos >= 1 && pos <= 4))

    #Set scipen to a very high value so large numbers are not written using scientific notation
    old.scipen <- options("scipen")
    options(scipen=999)
    on.exit(options(scipen=old.scipen), add=TRUE)

    #Define colors
    if(is.null(label.col)) label.col <- col
    if(is.null(tick.col)) tick.col <- col
    if(is.null(minor.tick.col)) minor.tick.col <- col


    #Define the x position where labels will be plotted
    lab.pos <- numeric()
    if(!is.null(custom.labels)) {
      lab.pos <- c(lab.pos, custom.labels)
    }
    if(label.exon.ends) {
      lab.pos <- c(lab.pos, GenomicRanges::start(genoplot$exons), GenomicRanges::end(genoplot$exons))
    }
    if(label.every.x) {
      first.tick <- floor(start(genoplot$regions[1])/tick.dist)*tick.dist
      last.tick <- ceiling(end(genoplot$regions[length(genoplot$regions)])/tick.dist)*tick.dist
      lab.pos <- c(lab.pos, seq(from=first.tick, to=last.tick, by=tick.dist))
    }

    lab.pos <- sort(unique(lab.pos)) #This unique should not remove custom label names because they are the first in the vector

    #If no predefined names for the labels, create them
    if(is.null(names(lab.pos)) || any(names(lab.pos)=="")) {
      standard.names <- mapply(lab.pos, FUN = function(x) toLabel(genoplot=genoplot, n=x, digits = digits, add.units = add.units, shorten=shorten.labels, cDNA=cDNA))
      if(is.null(names(lab.pos))) {
        names(lab.pos) <- standard.names
      } else {
        names(lab.pos)[names(lab.pos)==""] <- standard.names[names(lab.pos)==""]
      }
    }


    #build a GRanges with the labels position for later filtering by overlaps
    lab.pos.gr <- regioneR::toGRanges(data.frame(genoplot$chromosome, lab.pos, lab.pos, labels=names(lab.pos), stringsAsFactors=FALSE))

    #Minor ticks
    if(add.minor.ticks) {
      first.tick <- floor(start(genoplot$regions[1])/minor.tick.dist)*minor.tick.dist
      last.tick <- ceiling(end(genoplot$regions[length(genoplot$regions)])/minor.tick.dist)*minor.tick.dist
      minor.ticks <- seq(from=first.tick, to=last.tick, by=minor.tick.dist)
      minor.ticks <- regioneR::toGRanges(genoplot$chromosome, minor.ticks, minor.ticks)
    }


    #Define the vertical positions
    if(position=="below.gene") {
      r1 <- 0
      r0 <- -1 * tick.len
      data.panel <- "ideogram"
      inverted <- FALSE
    } else {
      if(position=="above.gene") {
        r0 <- 1
        r1 <- 1 + tick.len
        data.panel <- "ideogram"
        inverted <- TRUE
      }
      #else {
      #  #Don't change anyting
      #}
    }

    for(i in seq_len(length(genoplot$regions.kp))) {
      kp <- genoplot$regions.kp[[i]]

      local.lab.pos <- IRanges::subsetByOverlaps(lab.pos.gr, kp$plot.region) #We need to subset because we cannot activate the clipping (we are plotting out of data.panels)


      if(add.minor.ticks) {
        local.minor.ticks <- IRanges::subsetByOverlaps(minor.ticks, kp$plot.region)
        if(inverted) {
          r0.min <- r0
          r1.min <- r0+(r1-r0)*minor.tick.len
        } else {
          r0.min <- r1-(r1-r0)*minor.tick.len
          r1.min <- r1
        }
        karyoploteR::kpSegments(kp, chr=genoplot$chromosome, x0=GenomicRanges::start(local.minor.ticks), x1=GenomicRanges::end(local.minor.ticks), y0=0, y1=1, ymin=0, ymax=1, col=minor.tick.col, clipping=FALSE, r0=r0.min, r1=r1.min, data.panel=data.panel, ...)
      }

      karyoploteR::kpSegments(kp, chr=genoplot$chromosome, x0=GenomicRanges::start(local.lab.pos), x1=GenomicRanges::end(local.lab.pos), y0=0, y1=1, ymin=0, ymax=1, col=tick.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, ...)
      if(add.horizontal.line) {
        karyoploteR::kpAbline(kp, h=ifelse(inverted, 0, 1), r0=r0, r1=r1, ymin=0, ymax=1, col=tick.col, data.panel = data.panel, ...)
      }
      karyoploteR::kpSegments(kp, chr=genoplot$chromosome, x0=GenomicRanges::start(local.lab.pos), x1=GenomicRanges::end(local.lab.pos), y0=0, y1=1, ymin=0, ymax=1, col=tick.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, ...)
      if(inverted) {
        karyoploteR::kpText(kp, chr=genoplot$chromosome, x=GenomicRanges::start(local.lab.pos), labels = local.lab.pos$labels, y=1, ymin=0, ymax=1, col=label.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, pos=ifelse(is.null(pos), 3, pos), cex=cex, ...)
      } else {
        karyoploteR::kpText(kp, chr=genoplot$chromosome, x=GenomicRanges::start(local.lab.pos), labels = local.lab.pos$labels, y=0, ymin=0, ymax=1, col=label.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, pos=ifelse(is.null(pos), 1, pos), cex=cex, ...)
      }
    }


}


#internal
#' @keywords internal
toLabel <- function(genoplot, n, add.units=TRUE, digits=2, shorten=TRUE, cDNA=FALSE) {
  if(add.units==TRUE) {
    units <- c("b", "Kb", "Mb")
  } else {
    units <- c("", "", "")
  }
  if(cDNA) n <- genomicTocDNA(genoplot, n)
  if(abs(n) < 1000 | shorten==FALSE) return(paste0(as.character(n), units[1]))
  if(abs(n) < 1000000) return(paste0(as.character(round(n/1000, digits=digits)), units[2])) #Kb
  return(paste0(as.character(round(n/1000000, digits=digits)), units[3])) #Mb
}


