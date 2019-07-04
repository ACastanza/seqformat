Load the script with source("/seqformat.R") from your R session.

This script requires at minimum, a GSEA compatible chip file to convert to gene symbols matching the correct ENSEMBL release.

R Package Dependencies:

From R base:

tools

From CRAN:

dplyr
tidyr
readr (optional)
tibble (there is a call to a tibble command, but I don’t declare the library, so it must be a requirement of something else)

 

From Bioconductor:

GEOquery
DESeq2
tximport
GenomicFeatures
rhdf5

