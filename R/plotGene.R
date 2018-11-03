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






#If the function name is plotGene, regions should be called exons.
#In any case it should accept a list of GRanges (or GRangesList) with
#coding and non-coding exons or a TxDb and transcript name?
#Or maybe the trasformation from transcript should be done outside
#in an exernal function
plotGene <- function(exons, plot.type="only.exons", main=NULL, plot.params=NULL, gene.plotter=gpAddGeneStructure, ...) {
  #TODO: Check parameters

  #Get plot params given the plot.type
  if(is.null(plot.params)) {
    plot.params <- getDefaultPlotParams(plot.type=plot.type)
  }

  exons <- toGRanges(exons) #So it's possible to give them in any valid format

  total.gene.region <- toGRanges(seqnames(exons[1]), start(exons[1]), end(exons[length(exons)]))
  introns <-subtractRegions(total.gene.region, exons)


  #TODO: To implement a zoom, we should intersect the exons (and introns!) with the zoom region and continue as usual but keeping their mcols

  #Compute the region size values for the selected plot.type
  #Case 1: only exons visible (remove introns)
  if(plot.type=="only.exons") {
    total.plot.length.bases <-sum(width(exons)) + #The actual regions to plot
                              2*plot.params$inner.margin.bases*(length(exons)-1) + #The margin in bases between the regions
                              2*plot.params$outer.margin.bases #The margin in the outer ends of the first and last region
    internal.margin <- plot.params$margin.between.regions*(length(exons)-1)
    available.space <- 1 - plot.params$leftmargin - plot.params$rightmargin - internal.margin
    exons.space.per.base <- available.space/total.plot.length.bases
    regions <- exons #Plot only the exons
    regions <- regioneR::joinRegions(regions, min.dist = 1)
    regions$space.per.base <- exons.space.per.base
  } else if(plot.type=="compressed.introns.proportional") {
    stop("Unimplemented plot.type")
  } else if(plot.type=="compressed.introns.all.equal") {
    stop("Unimplemented plot.type")
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

  #Create the plot
  #Start creating an empty karyoplot (zoomed in on the correct chromosome) and get the KaryoPlot object
  gp$global.kp <- karyoploteR::plotKaryotype(zoom=exons[1], plot.type=2, labels.plotter = NULL, ideogram.plotter = NULL, plot.params=plot.params) #The zoom used here is not used to plot anything, but sets kp in the correct state

  #And prepare the list of KaryoPlot objecrts needed to plot
  #Now create a modified copy of this KaryoPlot adapted to each plotted region
  gp$regions.kp <- list()
  #Iteratively make the left margin bigger to plot each region at its position
  last.region.end <- gp$global.kp$plot.params$leftmargin
  for(num.reg in seq_len(length(regions))) {
    r <- regions[num.reg]

    if(num.reg==1) { #If the first region, add the out margin on the left and do not add the margin between regions to the last left
      r <- extendRegions(r, extend.start = plot.params$outer.margin.bases, extend.end = plot.params$inner.margin.bases)
      region.size <- width(r)*r$space.per.base #The space per base may be different for each region
      left.margin <- last.region.end
    } else {
      if(num.reg==length(regions)) { #If it's the last region, add the out margin to the right
        r <- extendRegions(r, extend.start = plot.params$inner.margin.bases, extend.end = plot.params$outer.margin.bases)
        region.size <- width(r)*r$space.per.base #The space per base may be different for each region
        left.margin <- last.region.end + plot.params$margin.between.regions
      } else { #if any other region, just add the regular intron-exon margins
        r <- extendRegions(r, extend.start = plot.params$inner.margin.bases, extend.end = plot.params$inner.margin.bases)
        region.size <- width(r)*r$space.per.base #The space per base may be different for each region
        left.margin <- last.region.end + plot.params$margin.between.regions
      }
    }
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
