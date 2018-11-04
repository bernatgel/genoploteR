#' plotGene
#'
#' @description
#'
#' Plots the main title of the plot
#'
#' @details
#'
#' @usage plotGene()
#'
#' @param genoplot    a \code{GenoPlot} object returned by a call to \code{plotGene}
#' @param main (character) the main title of the plot
#' @param ...  any additional parameter to be passed to the text plotting. All R base graphics params are passed along.
#'
#' @return
#' invisibly returns the given GenoPlot object
#'
#' @seealso \code{\link{gpAddMainTitle}}
#'
#' @examples
#'
#'
#'
#'
#' @export plotGene


addMargins <- function(regs, inner.margin, outer.margin) {
  regs <- extendRegions(regs, inner.margin, inner.margin)
  regs[1] <- extendRegions(regs[1], extend.start = outer.margin - inner.margin)
  regs[length(regs)] <- extendRegions(regs[length(regs)], extend.end = outer.margin - inner.margin)
  return(regs)
}


selectPlotType <- function(show.introns=FALSE, compress.introns=TRUE, proportional.introns=TRUE) {
  if(show.introns==FALSE) return("only.exons")
  if(compress.introns==FALSE) return("all")
  if(proportional.introns==TRUE) return("compressed.introns.proportional")
  else return("compressed.introns.all.equal")
}


#If the function name is plotGene, regions should be called exons.
#In any case it should accept a list of GRanges (or GRangesList) with
#coding and non-coding exons or a TxDb and transcript name?
#Or maybe the trasformation from transcript should be done outside
#in an exernal function
plotGene <- function(exons, show.introns=FALSE, compress.introns=TRUE, proportional.introns=TRUE, intron.to.exon.ratio=0.1, intron.length=200, main=NULL, plot.params=NULL, gene.plotter=gpAddGeneStructure, ...) {
  #TODO: Check parameters

  plot.type <- selectPlotType(show.introns=show.introns, compress.introns=compress.introns, proportional.introns=proportional.introns)

  #Get plot params given the plot.type
  if(is.null(plot.params)) {
    plot.params <- getDefaultPlotParams(show.introns=show.introns, compress.introns=compress.introns, proportional.introns=proportional.introns)
  }

  exons <- sort(toGRanges(exons)) #So it's possible to give them in any valid format

  total.gene.region <- toGRanges(seqnames(exons[1]), start(exons[1]), end(exons[length(exons)]))
  introns <- subtractRegions(total.gene.region, exons)
  introns$type <- "intron"

  #TODO: If the intron between two exons is shorter than twice the inner margin, join the exons (and the intron) into a single region


  #TODO: To implement a zoom, we should intersect the exons (and introns!) with the zoom region and continue as usual but keeping their mcols

  #Compute the region size values for the selected plot.type
  #Case 1: only exons visible (remove introns)
  if(plot.type=="only.exons") {
    ex.regs <-  regioneR::joinRegions(exons, min.dist = 1)
    ex.regs <- addMargins(ex.regs, inner.margin = plot.params$inner.margin.bases, outer.margin = plot.params$outer.margin.bases)

    internal.margin <- plot.params$margin.between.regions*(length(exons)-1)
    available.space <- 1 - plot.params$leftmargin - plot.params$rightmargin - internal.margin

    exons.space.per.base <- available.space/sum(width(ex.regs))

    mcols(ex.regs) <- DataFrame(space.per.base=exons.space.per.base)
    regions <- ex.regs #Plot only the exons

  } else if(plot.type=="compressed.introns.proportional") {
    internal.margin <- plot.params$margin.between.regions*(length(exons)-1)
    available.space <- 1 - plot.params$leftmargin - plot.params$rightmargin - internal.margin

    ex.regs <- regioneR::joinRegions(exons, min.dist = 1) #Join contiguous coding and non-coding exons into single regions
    ex.regs <- addMargins(ex.regs, inner.margin = plot.params$inner.margin.bases, outer.margin = plot.params$outer.margin.bases)

    in.regs <- extendRegions(introns, -1*plot.params$inner.margin.bases, -1*plot.params$inner.margin.bases)

    exons.bases <- sum(width(ex.regs))
    introns.bases <- sum(width(in.regs))
    adjusted.total.bases <- exons.bases + introns.bases*intron.to.exon.ratio

    exons.space.per.base <- available.space/adjusted.total.bases
    introns.space.per.base <- exons.space.per.base*intron.to.exon.ratio

    mcols(ex.regs) <- DataFrame(space.per.base=exons.space.per.base)
    mcols(in.regs) <- DataFrame(space.per.base=introns.space.per.base)
    regions <- sort(c(ex.regs, in.regs))

  } else if(plot.type=="compressed.introns.all.equal") {
    internal.margin <- plot.params$margin.between.regions*(length(exons)-1)
    available.space <- 1 - plot.params$leftmargin - plot.params$rightmargin - internal.margin

    ex.regs <- regioneR::joinRegions(exons, min.dist = 1) #Join contiguous coding and non-coding exons into single regions
    ex.regs <- addMargins(ex.regs, inner.margin = plot.params$inner.margin.bases, outer.margin = plot.params$outer.margin.bases)

    in.regs <- extendRegions(introns, -1*plot.params$inner.margin.bases, -1*plot.params$inner.margin.bases)

    exons.bases <- sum(width(ex.regs))
    introns.bases <- length(in.regs)*intron.length
    adjusted.total.bases <- exons.bases + introns.bases

    exons.space.per.base <- available.space/adjusted.total.bases
    introns.space.per.base <- exons.space.per.base*intron.length/width(in.regs)

    mcols(ex.regs) <- DataFrame(space.per.base=exons.space.per.base)
    mcols(in.regs) <- DataFrame(space.per.base=introns.space.per.base)
    regions <- sort(c(ex.regs, in.regs))

  } else if(plot.type=="all") {
    internal.margin <- plot.params$margin.between.regions*(length(exons)-1)
    available.space <- 1 - plot.params$leftmargin - plot.params$rightmargin - internal.margin

    ex.regs <- regioneR::joinRegions(exons, min.dist = 1) #Join contiguous coding and non-coding exons into single regions
    ex.regs <- addMargins(ex.regs, inner.margin = plot.params$inner.margin.bases, outer.margin = plot.params$outer.margin.bases)

    in.regs <- extendRegions(introns, -1*plot.params$inner.margin.bases, -1*plot.params$inner.margin.bases)

    exons.bases <- sum(width(ex.regs))
    introns.bases <- sum(width(in.regs))
    adjusted.total.bases <- exons.bases + introns.bases

    exons.space.per.base <- available.space/adjusted.total.bases
    introns.space.per.base <- exons.space.per.base

    mcols(ex.regs) <- DataFrame(space.per.base=exons.space.per.base)
    mcols(in.regs) <- DataFrame(space.per.base=introns.space.per.base)
    regions <- sort(c(ex.regs, in.regs))

  } else {
    stop("Invalid plot.type ", plot.type)
  }


  #At this point we have a "regions" GRanges object with the plotted regions and for each one a space.per.base value

  #Create the GenoPlot object to return
  gp <- list()
  class(gp) <- "GenoPlot"

  #and set some of its elements
  gp$plot.params <- plot.params
  gp$plot.type <- plot.type
  gp$exons <- exons
  gp$introns <- introns
  gp$total.plot.region <- extendRegions(total.gene.region, extend.start = plot.params$outer.margin.bases, extend.end = plot.params$outer.margin.bases)

  gp$regions <- regions


  #TODO: Decide how to manage the genome specification. Should we use a genome? Could be useful if we want to automatically add the sequence...
  #Create a custom genome for karyoploteR with the plot region only (will speed up some functions)
  cust.genome <- gp$total.plot.region

  #Create the plot
  #Start creating an empty karyoplot (zoomed in on the correct chromosome) and get the KaryoPlot object
  gp$global.kp <- karyoploteR::plotKaryotype(zoom=exons[1], plot.type=2, labels.plotter = NULL, ideogram.plotter = NULL, plot.params=plot.params, genome=cust.genome) #The zoom used here is not used to plot anything, but sets kp in the correct state

  #And prepare the list of KaryoPlot objecrts needed to plot
  #Now create a modified copy of this KaryoPlot adapted to each plotted region
  gp$regions.kp <- list()
  #Iteratively make the left margin bigger to plot each region at its position
  last.region.end <- gp$global.kp$plot.params$leftmargin - plot.params$margin.between.regions #The first region will cancel this between regions margin
  for(num.reg in seq_len(length(regions))) {
    r <- regions[num.reg]
    region.size <- width(r)*r$space.per.base
    left.margin <- last.region.end + plot.params$margin.between.regions

    right.margin <- 1 - left.margin - region.size
    last.region.end <- 1 - right.margin

    #DEBUG
    gp$global.kp$beginKpPlot()
    abline(v=left.margin)
    abline(v=1-right.margin)
    text(y=gp$global.kp$chromosome.height/2, x=mean(c(left.margin, (1-right.margin))), labels = num.reg)
    gp$global.kp$endKpPlot()
    #END DEBUG

    #Build the region specific karyoplot
    reg.kp <- gp$global.kp
      #Set the plot.params
        reg.kp$plot.params <- gp$global.kp$plot.params
        reg.kp$plot.params$leftmargin <- left.margin
        reg.kp$plot.params$rightmargin <- right.margin

      #Set the plot region
        names(r) <- as.character(seqnames(r)) #needed for various karyoploteR functions
        reg.kp$plot.region <- r

      #Update the coordinate change functions
        #Warning: using triple semicolon to access the PRIVATE function from karyoploteR.
        #Doing it so it's not necessary for karyoploteR to export that function, which would
        #be only useful in extremely hacky use cases like this one
        reg.kp$coord.change.function <- karyoploteR:::getCoordChangeFunctions(reg.kp)$coordChangeFunction

    #Store the new kp in the list
    gp$regions.kp[[num.reg]] <- reg.kp
  }


  #Set some more elements in the GenoPlot object
  gp$beginGpPlot <- gp$global.kp$beginKpPlot
  gp$endGpPlot <- gp$global.kp$endGpPlot


  #At this point the gp object is finished.

  #As a convenience (and to not create an empty plot by default) optionally plot a few elements

  #TODO: Plot the exon/intron structure
  if(!is.null(gene.plotter)) {
    gene.plotter(gp, ...)
  }

  #TODO: Plot the name of the gene/transcript...

  #Plot the main title
  if(!is.null(main)) {
    gpAddMainTitle(gp, main=main, ...)
  }



  #Everything finished. Return the GenoPlot object
  return(gp)


}
