#' cDNAToGenomic
#'
#' @note: TODO documentation
#'
#' @export cDNAToGenomic

cDNAToGenomic <- function(genoplot, cdna.pos) {
  if(any(cdna.pos<1)) stop("cDNA position must be greater or equal to 1")
  exon.cum.width <- c(0, cumsum(width(genoplot$exons)))
  ex.num <- sapply(cdna.pos, function(p) {return(which(exon.cum.width>=p)[1])})-1
  #if(is.na(ex.num)) stop("The provided cDNA position is out the transcript range.")
  return(start(genoplot$exons)[ex.num] - exon.cum.width[ex.num] + cdna.pos - 1)
}

