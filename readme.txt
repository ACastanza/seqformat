The point of this script is that if you carefully follow the prompts using the information provided, you should be able to process any (bulk) RNA-seq dataset into normalized, properly formatted, data files ready to be used for GSEA (http://www.gsea-msigdb.org)

Load the script with source("/seqformat.R") from your R session.

This script requires at minimum, a GSEA compatible chip file to convert to gene symbols matching the correct ENSEMBL release. ENSEMBL97 based CHIP files targeting MSigDB7 are available for Human and Mouse Ensembl IDs, and Gene Symbols in the Additional Datafiles directory.

Warning: Orthology conversions after normalization are tentatively supported starting in the v2.0 releases. Do not attempt orthology conversion in prior builds as the results are not statistically valid.

R Package Dependencies:

From R base:

tools (Required)

From CRAN:

dplyr (Required. Essential!)
tidyr (required only if necessicary to split gene identifiers merged like ENSG000001_GeneSymbol, from a weird pipleine)
readr (optional)
tibble (there is a call to a tibble command, but I don’t declare the library, so it must be a requirement loaded from something else)



From Bioconductor:

GEOquery (required only if you want to pull datasets directly from GEO)
DESeq2 (required only for data that is not already normalized)
tximport (required only for transcript level features)
GenomicFeatures (required only if building tx2gene files from GTF/GFF3 for transcript data)
rhdf5 (required only if parsing kallisto abundance.h5 files)
biomaRt (required only if building CHIP file from Biomart rather than supplying one).
