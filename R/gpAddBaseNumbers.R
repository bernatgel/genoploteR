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



#Funcio per a dibuixar un num i un palet en una pos determinada (però fent servir kpSegment i kpText)

#posicions dels nums allà on sigui, però la y determina el punt "top" de les etiquetes. Així, ideogram, y=0 les dibuixa o les volem normalment

#Pero data,panel i y es poden modificar lliurement (i r0 i r1 tambe)

#Side effect de fer servir kpSegment: si canviem de data.panel, la longitud dels ticks canvia. Ho volem? O fem la transformació fent servir ccf
#i crida directe a segments. (potser la primera opció és coherent amb tota la resta de funcions.)

#D'altra banda, funció per a determinar les posicions dels ticks grans  i petits.


#TODO: Reimplement here from scratch to:
# - Allow moving the ruler around
# - Getting at least one tick per exon (first base ideally?)
# - Possibility of setting the numbers to only to ticks per exon (first and last base)
# - Selecting between genomic and cDNA bases (genomic could actually use the genomic2cDNA function to compute?)
# - Decide whether to label number introns

#internal
toLabel <- function(n, add.units=TRUE, digits=2, shorten=TRUE, cDNA=FALSE) {
  if(add.units==TRUE) {
    units <- c("b", "Kb", "Mb")
  } else {
    units <- c("", "", "")
  }
  if(cDNA) n <- genomicTocDNA(gp, n)
  if(abs(n) < 1000 | shorten==FALSE) return(paste0(as.character(n), units[1]))
  if(abs(n) < 1000000) return(paste0(as.character(round(n/1000, digits=digits)), units[2])) #Kb
  return(paste0(as.character(round(n/1000000, digits=digits)), units[3])) #Mb
}



gpAddBaseNumbers <- function(genoplot, label.exon.ends=TRUE, label.every.x=FALSE,
                             add.horizontal.line=FALSE,
                             custom.labels=NULL, ticks.in.introns=TRUE,
                             position="below.gene", inverted=FALSE,
                             add.units=TRUE, digits=2, shorten.labels=TRUE, cDNA=FALSE, cex=0.6,
                             tick.dist=500, tick.len=0.2, tick.col=NULL,
                             minor.tick.dist=100, minor.tick.len=0.1, minor.tick.col=NULL,
                             r0=0, r1=1, data.panel=1, ...) {

    #TODO: Check arg values
    position <- match.arg(position, c("custom", "above.gene", "below.gene"))

    old.scipen <- options("scipen")
    options(scipen=999)
    on.exit(options(scipen=old.scipen), add=TRUE)

    #Define the x position where labels will be plotted
    if(is.null(custom.labels)) {
      lab.pos <- numeric()
      if(label.exon.ends) {
        lab.pos <- unique(sort(c(lab.pos, start(genoplot$exons), end(genoplot$exons))))
      }
      if(label.every.x) {

      }
    } else {
      lab.pos <- custom.labels
    }

    #If necessary, create the labels themselves
    if(is.null(names(lab.pos))) {
      names(lab.pos) <- mapply(lab.pos, FUN = toLabel, digits = digits, add.units = add.units, shorten=shorten.labels, cDNA=cDNA)
    }

    #build a GRanges with the labels position for later filtering by overlaps
    lab.pos.gr <- toGRanges(data.frame(genoplot$chromosome, lab.pos, lab.pos, labels=names(lab.pos), stringsAsFactors=FALSE))


    #Define the vertical positions
    if(position=="below.gene") {
      r1 <- 0
      r0 <- -1*tick.len
      data.panel <- "ideogram"
      inverted <- FALSE
    } else {
      if(position=="above.gene") {
        r0 <- 1
        r1 <- 1+tick.len
        data.panel <- "ideogram"
        inverted <- TRUE
      } else {
        #Don't change anyting
      }
    }

    for(i in seq_len(length(genoplot$regions.kp))) {
      kp <- genoplot$regions.kp[[i]]
      local.lab.pos <- subsetByOverlaps(lab.pos.gr, kp$plot.region)
      if(add.horizontal.line) {
        kpAbline(kp, h=ifelse(inverted, 0, 1), r0=r0, r1=r1, ymin=0, ymax=1, col="red", data.panel = "ideogram", ...)
      }
      kpSegments(kp, chr=genoplot$chromosome, x0=start(local.lab.pos), x1=end(local.lab.pos), y0=0, y1=1, ymin=0, ymax=1, col="red", clipping=FALSE, r0=r0, r1=r1, data.panel="ideogram", ...)
      if(inverted) {
        kpText(kp, chr=genoplot$chromosome, x=start(local.lab.pos), labels = local.lab.pos$labels, y=1, ymin=0, ymax=1, col="red", clipping=FALSE, r0=r0, r1=r1, data.panel="ideogram", pos=3, cex=cex, ...)
      } else {
        kpText(kp, chr=genoplot$chromosome, x=start(local.lab.pos), labels = local.lab.pos$labels, y=0, ymin=0, ymax=1, col="red", clipping=FALSE, r0=r0, r1=r1, data.panel="ideogram", pos=1, cex=cex, ...)
      }
    }


}

# gpAddBaseNumbers <- function(genoplot, tick.dist=500, tick.len=5, add.units=FALSE,
#                              digits=2, minor.ticks=TRUE,
#                              minor.tick.dist=100, minor.tick.len=2,  cex=0.5,
#                              tick.col=NULL, minor.tick.col=NULL, clipping=TRUE,  ...) {
#   for(i in seq_len(length(genoplot$regions.kp))) {
#     kpAddBaseNumbers(genoplot$regions.kp[[i]], tick.dist=tick.dist, tick.len=tick.len, add.units=add.units,
#                      digits=digits, minor.ticks=minor.ticks,
#                      minor.tick.dist=minor.tick.dist, minor.tick.len=minor.tick.len,  cex=cex,
#                      tick.col=tick.col, minor.tick.col=minor.tick.col, clipping=clipping,  ...)
#   }
# }
#


