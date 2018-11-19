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
                             add.horizontal.line=FALSE, add.minor.ticks=TRUE,
                             custom.labels=NULL, ticks.in.introns=TRUE,
                             position="below.gene", inverted=FALSE,
                             add.units=TRUE, digits=2, shorten.labels=TRUE, cDNA=FALSE,
                             cex=0.6, label.col=NULL,
                             tick.dist=500, tick.len=0.2, tick.col=NULL,
                             minor.tick.dist=100, minor.tick.len=0.5, minor.tick.col=NULL,
                             r0=0, r1=1, data.panel=1,
                             col="black", pos=NULL, ...) {

    #TODO: Check arg values
    position <- match.arg(position, c("custom", "above.gene", "below.gene"))

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
      lab.pos <- c(lab.pos, start(genoplot$exons), end(genoplot$exons))
    }
    if(label.every.x) {
      first.tick <- floor(start(gp$regions[1])/tick.dist)*tick.dist
      last.tick <- ceiling(end(gp$regions[length(gp$regions)])/tick.dist)*tick.dist
      lab.pos <- c(lab.pos, seq(from=first.tick, to=last.tick, by=tick.dist))
    }

    lab.pos <- sort(unique(lab.pos)) #This unique should not remove custom label names because they are the first in the vector

    #If no predefined names for the labels, create them
    if(is.null(names(lab.pos)) || any(names(lab.pos)=="")) {
      standard.names <- mapply(lab.pos, FUN = toLabel, digits = digits, add.units = add.units, shorten=shorten.labels, cDNA=cDNA)
      if(is.null(names(lab.pos))) {
        names(lab.pos) <- standard.names
      } else {
        names(lab.pos)[names(lab.pos)==""] <- standard.names[names(lab.pos)==""]
      }
    }


    #build a GRanges with the labels position for later filtering by overlaps
    lab.pos.gr <- toGRanges(data.frame(genoplot$chromosome, lab.pos, lab.pos, labels=names(lab.pos), stringsAsFactors=FALSE))

    #Minor ticks
    if(add.minor.ticks) {
      first.tick <- floor(start(gp$regions[1])/minor.tick.dist)*minor.tick.dist
      last.tick <- ceiling(end(gp$regions[length(gp$regions)])/minor.tick.dist)*minor.tick.dist
      minor.ticks <- seq(from=first.tick, to=last.tick, by=minor.tick.dist)
      minor.ticks <- toGRanges(genoplot$chromosome, minor.ticks, minor.ticks)
    }


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

      local.lab.pos <- subsetByOverlaps(lab.pos.gr, kp$plot.region) #We need to subset because we cannot activate the clipping (we are plotting out of data.panels)


      if(add.minor.ticks) {
        local.minor.ticks <- subsetByOverlaps(minor.ticks, kp$plot.region)
        if(inverted) {
          r0.min <- r0
          r1.min <- r0+(r1-r0)*minor.tick.len
        } else {
          r0.min <- r1-(r1-r0)*minor.tick.len
          r1.min <- r1
        }
        kpSegments(kp, chr=genoplot$chromosome, x0=start(local.minor.ticks), x1=end(local.minor.ticks), y0=0, y1=1, ymin=0, ymax=1, col=minor.tick.col, clipping=FALSE, r0=r0.min, r1=r1.min, data.panel=data.panel, ...)
      }

      kpSegments(kp, chr=genoplot$chromosome, x0=start(local.lab.pos), x1=end(local.lab.pos), y0=0, y1=1, ymin=0, ymax=1, col=tick.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, ...)
      if(add.horizontal.line) {
        kpAbline(kp, h=ifelse(inverted, 0, 1), r0=r0, r1=r1, ymin=0, ymax=1, col=tick.col, data.panel = data.panel, ...)
      }
      kpSegments(kp, chr=genoplot$chromosome, x0=start(local.lab.pos), x1=end(local.lab.pos), y0=0, y1=1, ymin=0, ymax=1, col=tick.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, ...)
      if(inverted) {
        kpText(kp, chr=genoplot$chromosome, x=start(local.lab.pos), labels = local.lab.pos$labels, y=1, ymin=0, ymax=1, col=label.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, pos=ifelse(is.null(pos), 3, pos), cex=cex, ...)
      } else {
        kpText(kp, chr=genoplot$chromosome, x=start(local.lab.pos), labels = local.lab.pos$labels, y=0, ymin=0, ymax=1, col=label.col, clipping=FALSE, r0=r0, r1=r1, data.panel=data.panel, pos=ifelse(is.null(pos), 1, pos), cex=cex, ...)
      }
    }


}
