#Scratch Pad
library(genoploteR)

coding.exons <- toGRanges(rep("chr1",3), c(50, 100, 200), c(70, 130, 220))
non.coding.exons <- toGRanges(rep("chr1", 3), c(20, 45, 221), c(25, 49, 270))
strand <- "+"



coding.exons$type <- "coding"
non.coding.exons$type <- "non.coding"

exons <- sort(toGRanges(c(coding.exons, non.coding.exons), "hg19"))


gp <- plotGene(exons=exons)
gpAddGeneStructure(gp)

plotGene(exons = exons)


#Example with pasilla data bam coverage
exons <- toGRanges(c("chr4:340600-340900", "chr4:341300-341700", "chr4:343500-345500", "chr4:348400-348800"))
exons$type <- "coding"
gp <- plotGene(exons)



#data background
for(i in seq_len(length(regions))) {

  kpDataBackground(regions.kp[[i]], data.panel = 1, color = rainbow(length(regions))[i])
  kpAbline(regions.kp[[i]], v=start(regions[i]), h=0.5)
  #r <- regions[i]
  #height <- ifelse(r$type=="coding", 0.8, 0.4)
  #kpAbline(regions.kp[[i]], h=0.5, clipping=TRUE)
  #kpRect(regions.kp[[i]], data=extendRegions(regions[i], extend.end = 100), y0=0.5-height/2, y1=0.5+height/2)
}


roxygen2::roxygenise()



sessionInfo()
