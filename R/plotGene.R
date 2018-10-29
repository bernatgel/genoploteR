#plotGene



library(karyoploteR)


setGenomeToGRanges <- function(gr, genome) {
  if(!is.null(genome)) {
    if(is.character(genome)) genome <- characterToBSGenome(genome)
    if(!methods::is(genome, "BSgenome")) {
      warning("Invalid 'genome' argument. Ignoring.")
      genome <- NULL
    }
  }
  if(!is.null(genome)) {
    GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(genome)
    GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::seqlevels(genome)
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(genome)
  }
  return(gr)
}

coding.exons <- toGRanges(rep("chr1",3), c(50, 100, 200), c(70, 130, 220))
non.coding.exons <- toGRanges(rep("chr1", 3), c(20, 45, 221), c(25, 49, 270))
strand <- "+"



coding.exons$type <- "coding"
non.coding.exons$type <- "non.coding"

exons <- sort(setGenomeToGRanges(c(coding.exons, non.coding.exons), "hg19"))



#total.plot.region <- regioneR::extendRegions(total.gene.region, extend.start = out.margin, extend.end = out.margin)

#TODO: Move o its own file
getDefaultPlotParams <- function(plot.type) {
  #TODO: check parameters are valid

  if(plot.type=="only.exons") {
   pp <- karyoploteR::getDefaultPlotParams(plot.type=2)
   #TODO: Modify any default values to adjust to our needs

   #Add genoploteR specific values
   pp$outer.margin.bases <- 10
   pp$inner.margin.bases <- 5
   pp$margin.between.regions <- 0.02

   return(pp)
  }
  stop("getDefaultPlotParams: plot type not found: ", plot.type)
}


#If the function name is plotGene, regions should be called exons.
#In any case it should accept a list of GRanges (or GRagesList) with
#coding and non-coding exons or a TxDb and transcript name?
#Or maybe the trasformation from transcript should be done outside
#in an exernal function
plotGene <- function(exons, plot.type="only.exons", main=NULL, plot.params=NULL) {
  #Check parameters

  #Get plot params given the plot.type
  if(is.null(plot.params)) {
    plot.params <- getDefaultPlotParams(plot.type=plot.type)
  }

  total.gene.region <- toGRanges(seqnames(exons[1]), start(exons[1]), end(exons[length(exons)]))
  introns <-subtractRegions(total.gene.region, exons)


  #TODO: To implement a zoom, we should intersect the exons with the zoom region and continue as usual

  #Compute the region size values for the selected plot.type
  #Case 1: only exons visible (remove introns)
  if(plot.type=="only.exons") {
    total.plot.length.bases <-sum(width(exons)) + #The actual regions to plot
                              2*plot.params$inner.margin.bases*(length(exons)-1) + #The margin in bases between the regions
                              2*plot.params$outer.margin.bases #The margin in the outer ends of the first and last region
    internal.margin <- margin.between.regions*(length(exons)-1)
    available.space <- 1 - plot.params$leftmargin - plot.params$rightmargin - internal.margin
    exons.space.per.base <- available.space/total.plot.length.bases
    regions <- exons #Plot only the exons
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
  gp$global.kp <- plotKaryotype(zoom=exons[1], plot.type=2, labels.plotter = NULL, ideogram.plotter = NULL, plot.params=plot.params) #The zoom used here is not used to plot anything, but sets kp in the correct state

  #And prepare the list of KaryoPlot objecrts needed to plot
  #Now create a modified copy of this KaryoPlot adapted to each plotted region
  gp$regions.kp <- list()
  #Iteratively make the left margin bigger to plot each region at its position
  last.region.end <- kp$plot.params$leftmargin
  for(num.reg in seq_len(length(regions))) {
    r <- regions[num.reg]

    if(num.reg==1) { #If the first region, add the out margin on the left and do not add the margin between regions to the last left
      r <- extendRegions(r, extend.start = out.margin.bases, extend.end = intron.exon.margin.bases)
      region.size <- width(r)*space.per.base
      left.margin <- last.region.end
    } else {
      if(num.reg==length(regions)) { #If it's the last region, add the out margin to the right
        r <- extendRegions(r, extend.start = intron.exon.margin.bases, extend.end = out.margin.bases)
        region.size <- width(r)*space.per.base
        left.margin <- last.region.end + margin.between.regions
      } else { #if any other region, just add the regular intron-exon margins
        r <- extendRegions(r, extend.start = intron.exon.margin.bases, extend.end = intron.exon.margin.bases)
        region.size <- width(r)*space.per.base
        left.margin <- last.region.end + margin.between.regions
      }
    }
    right.margin <- 1 - left.margin - region.size
    last.region.end <- 1 - right.margin

    #DEBUG
    gp$global.kp$beginKpPlot()
    abline(v=left.margin)
    abline(v=1-right.margin)
    text(y=kp$chromosome.height/2, x=mean(c(left.margin, (1-right.margin))), labels = num.reg)
    gp$global.kp$endKpPlot()
    #END DEBUG

    #Build the region specific karyoplot
    reg.kp <- gp$global.kp
      #Set the plot.params
        reg.kp$plot.params <- kp$plot.params
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
    gp$regions.kp[[i]] <- reg.kp
  }


  #Set some more elements in the GenoPlot object
  gp$beginGpPlot <- gp$global.kp$beginKpPlot
  gp$endGpPlot <- gp$global.kp$endGpPlot


  #At this point the gp object is finished.

  #As a convenience (and to not create an empty plot by default) optionally plot a few elements

  #TODO: Plot the exon/intron structure

  #TODO: Plot the name of the gene/transcript...

  #Plot the main title
  if(!is.null(main)) {
    gpAddMainTitle(gp, main=main, ...)
  }



  #Everything finished. Return the GenoPlot object
  return(gp)


}



MUNTAR una funció tipus prepare parameters 2 que
preprocessi les dades (transfrmció de cdna a genomic? o fer-ho fora en un a una funcio auxiliar?)
Pero aquesta hauria d etraduir el data.panel="gene" a plot.panel="ideogram"

plotGene(exons = exons)


#data background
for(i in seq_len(length(regions))) {

  kpDataBackground(regions.kp[[i]], data.panel = 1, color = rainbow(length(regions))[i])
  kpAbline(regions.kp[[i]], v=start(regions[i]), h=0.5)
  #r <- regions[i]
  #height <- ifelse(r$type=="coding", 0.8, 0.4)
  #kpAbline(regions.kp[[i]], h=0.5, clipping=TRUE)
  #kpRect(regions.kp[[i]], data=extendRegions(regions[i], extend.end = 100), y0=0.5-height/2, y1=0.5+height/2)
}

