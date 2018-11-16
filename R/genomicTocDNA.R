#' genomicTocDNA
#'
#' @note: TODO documentation
#'
#' @export genomicTocDNA

last <- function(x) {return(x[length(x)])}


genomicTocDNA <- function(genoplot, g.pos) {
  #Check it's in range (or return the correct value, out of range?)
  exon.cum.width <- c(0,cumsum(width(gp$exons)))
  ex.num <- last(which(start(genoplot$exons)<=g.pos))
  if(g.pos<=end(genoplot$exons)[ex.num]) { #If it's exonic
    return(g.pos - start(genoplot$exons[ex.num]) + exon.cum.width[ex.num] + 1)
  } else { #if it's intronic (or UTR 3')
    #TODO: Depending on the distance to the previous and next exons, define as XXX+yy or ZZZ-ww
    return(paste0(exon.cum.width[ex.num+1], "+", g.pos - end(genoplot$exons[ex.num]) ))
  }
}

