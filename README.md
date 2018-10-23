[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)


# genoploteR - An R package to plot data on genes

## Description

This is the initial work on a package to plot data on genes. It is based on and heavily influenced by [karyoploteR](https://github.com/bernatgel/karyoploteR)). 

The idea is to create a package internally using the karyoploteR infrastructure to plot any type of data on genes. The user should provide somehow (manually, from a TxDb object...) the structure of a transcript or gene (exons, introns...) and be able to plot data on it using either low level functions such as points or lines or higher level functions to plot more specialized data (BAM coverage, mutations...).

One of the key features of the package is the ability to either hide or shrink the introns so the exonic data is clearly visible.

While the first versions will need the user to provide everything in genomic coordinates the aim is to be able to internally convert c. or p. notation to genomic coordinates.

The aim is to offload as much work as possible on karyoploteR and to have a range of plotting function available equivalent to it.

**Note:** Although the package is aimed at plotting data on genes it can be used to plot data on any regions as is not restricted to exons in any way.

## Expected use cases

There are a few use cases in our group, for example 

  * plotting the coverage level of a gene when sequencing a gene panel 
  * contextualizing the position of an unclassified variant on a transcript with respect to other known variants, protein domains...
  * plotting the distribution of pathogenic and non-pathogenic variants on a gene
  
  
but the package should be flexible enough to basically plot anything on a gene

