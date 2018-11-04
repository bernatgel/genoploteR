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

library(pasillaBamSubset) #A package with 2 example bam files
un1.bam.file <- untreated1_chr4() # get the name of the first bam
un3.bam.file <- untreated3_chr4() #and the name of the second

kp <- kpPlotBAMCoverage(kp, data = un1.bam.file) #Warning and does not plot. region too large.
kp <- kpPlotBAMCoverage(kp, data = un1.bam.file, max.valid.region.size=2000000)


exons <- toGRanges(c("chr4:340600-340900", "chr4:341300-341700", "chr4:343500-345000", "chr4:348400-348800"))
exons$type <- "coding"

pp <- getDefaultPlotParams(show.introns = TRUE)
pp$inner.margin.bases <- 50
pp$outer.margin.bases <- 200
gp <- plotGene(exons, show.introns = TRUE, intron.to.exon.ratio = 0.3, plot.params = pp)

gp <- plotGene(exons, show.introns = TRUE, compress.introns = FALSE, proportional.introns = TRUE, intron.length = 600)
gpPlotBAMCoverage(gp, data=un1.bam.file)
gpSegments(gp, chr="chr4", x0=340600, y0=1, x1=348800, y1=0)
gpAddBaseNumbers(gp)
gpAbline(gp, chr="chr4", v=c(340600, 340900, 340940))
width(gp$regions)



roxygen2::roxygenise()

library(TxDb.)

sessionInfo()
