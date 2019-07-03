cat("We're going to interactively process your RNA-seq data to output GSEA Compatible files.\n")
cat("Lets get started\n")
cat("Loading dplyr library for data formatting\n")
library("dplyr")
cat("Loaded dplyr for data structuring\n")
cat("\n")

if(exists("altnorm") == FALSE){altnorm <- FALSE}
## before running script set altnorm <- "userlog" or "usevst" or "useall" to output additional normalize using deseq2's rlog or vst functions

txlevel <- FALSE
FAIL <- FALSE
isnormalized <- FALSE
TYPE <- NULL
txtype <- 0
is_directory <- FALSE
getgeofiles <- FALSE

#Get Input Files and Prompt User for Necessary Information
#Set Working Directory

path <- readline(prompt=("Drop a directory into R window to use as the output folder or enter directory path: "))
setwd(path)
cat("Done\n")
cat("\n")

getgeofiles <- askYesNo("Attempt to get data directly from GEO \"supplementary files\"?")
  if(getgeofiles == TRUE){
  library("GEOquery")
  library("tools")
  geoid <- readline(prompt=("Enter the GEOID for the datafiles (eg: GSE38786): "))
  geooutfiles <- getGEOSuppFiles(geoid, makeDirectory = TRUE, baseDir = getwd(), fetch_files = TRUE, filter_regex = NULL)
  cat("Getting experiment information from Series Matrix...\n")
  gsefile <- getGEO(geoid)
  phenotypedata <- as.data.frame(pData(phenoData(gsefile[[1]]))[,c("geo_accession","title")], stringsAsFactors = FALSE)
  rownames(phenotypedata) <- 1:nrow(phenotypedata)
  geostructure <- pData(phenoData(gsefile[[1]]))
  cat("\n")
  print(phenotypedata)
  cat("\n")
  cat("Files acquired from GEO:\n")
  geoimporttable <- as.data.frame(rownames(geooutfiles), stringsAsFactors=FALSE)
  print(rownames(geooutfiles))
  cat("\n")
  geoselected <- readline(prompt=("Select which GEO file to use for downstream processing: "))
    geoselectednumber <- match(geoselected, cbind(rownames(geoimporttable),geoimporttable)[,1])
    if(is.na(geoselectednumber) == TRUE){
      geoselectednumber <- match(geoselected, cbind(rownames(geoimporttable),geoimporttable)[,2])
      }
  outfile <- rownames(geooutfiles)[geoselectednumber]
  if(file_ext(rownames(geooutfiles)[geoselectednumber]) == "tar"){
    untar(rownames(geooutfiles)[geoselectednumber], exdir=paste0(dirname(rownames(geooutfiles)[geoselectednumber]),"/",basename(tools::file_path_sans_ext(rownames(geooutfiles)[geoselectednumber]))))
    outfile <- paste0(dirname(rownames(geooutfiles)[geoselectednumber]),"/",basename(tools::file_path_sans_ext(rownames(geooutfiles)[geoselectednumber])))
    }
  genematrix <- outfile
  cat("\n")
  cat("Experiment Imported.\n")
  cat("\n")
  findnormal <- apply(geostructure,2,function(x){grepl("normalized|Normalized|NORMALIZED|normalised|Normalised|NORMALISED",x)})
  if(any(findnormal) == TRUE){
    cat("Series Matrix implies that this data might be ALREADY NORMALIZED:\n")
    print(unique(geostructure[findnormal]))
    isnormalized <- askYesNo("Is this data already normalized?")
        if(isnormalized == TRUE){
        NORM <- TRUE}
  } else if (any(findnormal) == FALSE){
      cat("Series Matrix implies that this data might require normalization.\n")
      cat("We'll prompt you about this later.\n")}
  findcounts <- apply(geostructure,2,function(x){grepl("counts|Counts",x)})
  if(any(findcounts) == TRUE){
    cat("Series Matrix implies that this data consists of COUNTS:\n")
    print(unique(geostructure[findcounts]))
      iscounts <- askYesNo("Do you agree that this is gene COUNTS data?")
          if(iscounts == TRUE){
          countsdetected <- TRUE
            if(isnormalized == TRUE){
            TYPE <- "NORMCOUNTS"}
            if(isnormalized == FALSE){
            TYPE <- "RAWCOUNTS"}}
    } else if (any(findcounts) == FALSE){
        cat("We couldn't detect the datatype\n")
        cat("We'll prompt you to manually select datatype next.\n")}
    findtxquant<- apply(geostructure,2,function(x){grepl("almon|ailfish|allisto",x)})
    if(any(findtxquant) == TRUE){
      cat("Series Matrix implies that this data might be transcript level quantifications:\n")
      print(unique(geostructure[findtxquant]))
      cat("\n")
      cat("Validate the presence of transcript level quantifications in data files, then follow prompts.\n")
      cat("\n")
      }
  }

istx <- askYesNo("Is your dataset transcript level abundance measurements from Salmon/Kallisto/Sailfish?")
  if(istx == FALSE){
  isrsem <- askYesNo("Is your dataset raw abundance measurements from RSEM?")
    if(isrsem == TRUE){
    if(getgeofiles == TRUE){dir <- dirname(genematrix)} #GEO
    txtype <- 4
    txlevel <- TRUE
    TYPE <- "RSEM"
    NORM <- FALSE
    library("tximport")}}
if(istx == TRUE){
TYPE <- "TXABUNDANCE"
txlevel <- TRUE
library("tximport")
NORM <- FALSE
iscounts <- FALSE
cat("\n")
cat("To process transcript alignments we need either an existing \"tx2gene\" file or a Gencode GTF/GFF3.\n")
tx2geneexists <- askYesNo("Do you have an existing tx2gene annotation file you wish to provide?")
if(tx2geneexists == TRUE){
tx2genepath <- readline(prompt=("Drop Your tx2gene file into R Window or Enter full path to file: "))
cat("\n")
library("tools")
tx2genetype <- file_ext(gsub(".gz$", "", tx2genepath))
  if(tx2genetype == "csv"){
  tx2gene <- read.table(tx2genepath, sep=",", header=T)
  cat("Identifiers:\n")
  print(as.data.frame(head(tx2gene[,1:2],3)))
  txhasversions <- askYesNo("Do your TRANSCRIPT IDs have decimal versions? (eg. ENST00000342771.9)?")
    if(txhasversions == TRUE){
    tx2gene[,1] <- gsub("\\..*","",tx2gene[,1])
    txhasversions <- FALSE}
  genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)?")
    if(genehasversions == TRUE){
    tx2gene[,2] <- gsub("\\..*","",tx2gene[,2])
    genehasversions <- FALSE}
  }
  if(tx2genetype == "txt" | tx2genetype == "tsv"| tx2genetype == "tabular"){
  tx2gene <- read.table(tx2genepath, sep="\t", header=T)
  print(head(tx2gene))
  txhasversions <- askYesNo("Do your TRANSCRIPT IDs have decimal versions? (eg. ENST00000342771.9)?")
    if(txhasversions == TRUE){
    tx2gene[,1] <- gsub("\\..*","",tx2gene[,1])
    txhasversions <- FALSE}
  genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)?")
    if(genehasversions == TRUE){
    tx2gene[,2] <- gsub("\\..*","",tx2gene[,2])
    genehasversions <- FALSE}
  }
tx2gene <- distinct(tx2gene)
print(head(tx2gene))
cat("This should be formatted as a two column file with TX IDs on the left and Gene IDs on the right and no version decimals.\n")
tx2genebuild <- FALSE
} else if (tx2geneexists == FALSE){
cat("Ok... \n")
cat("We can attempt to construct the tx2gene file automatically if you provide the GTF/GFF3 file used for your quantification. \n")
cat("This requires the GenomicFeatures Package from Bioconductor to be available.\n")
tx2genebuild <- askYesNo("(This has only been tested for Gencode Transcriptomes) Continue?")
}
  if(tx2genebuild == TRUE){
  library("GenomicFeatures")
  txgft <- readline(prompt=("Drop Your GTF/GFF3 file into R Window or Enter full path to file: "))
  TxDb <- makeTxDbFromGFF(file = txgft)
  k <- keys(TxDb, keytype = "TXNAME")
  tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
  print(head(tx2gene,3))
  txhasversions <- askYesNo("Do your TRANSCRIPT IDs have decimal versions? (eg. ENST00000342771.9)?")
    if(txhasversions == TRUE){
    tx2gene[,1] <- gsub("\\..*","",tx2gene[,1])}
  genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)?")
    if(genehasversions == TRUE){
    tx2gene[,2] <- gsub("\\..*","",tx2gene[,2])
    genehasversions <- FALSE}
  tx2geneexists <- TRUE
  cat("The Tx2Gene file should be formatted as a two column file with TX IDs on the left and Gene IDs on the right.\n")
  tx2gene <- distinct(tx2gene)}
  if(tx2geneexists == FALSE){
  cat("We can't continue with transcript level data without transcript to gene mappings...\n")
    FAIL <- TRUE
    if (FAIL == TRUE) {stop("Necessary files are not available so processing was terminated.")}}
cat("\n")
cat("tximport supports the following platforms:\n")
cat("\n")
cat("[1] Salmon\n")
cat("[2] Sailfish\n")
cat("[3] Kallisto\n")
cat("[4] RSEM\n")
#cat("[5] Stringtie\n")
#cat("[6] Generic Transcript Level Quantification Table\n")
#cat("WARNING: ONLY SALMON [1], Sailfish [2], and KALLISTO [3] ARE CURRENTLY SUPPORTED\n")
cat("\n")
txtype <- readline(prompt=("Select your platform by entering the cooresponding >>NUMBER<< (without brackets): "))
cat("\n")
if(getgeofiles == TRUE){dir <- dirname(genematrix)} #GEO
if (txtype==1){
  if(getgeofiles == FALSE){
  dir <- readline(prompt=("Drop a directory containing Salmon output into R Window or Enter Directory Path: "))
  }
files <- list.files(dir, recursive=TRUE, full.names = TRUE, pattern = ".sf|.sf.gz")
names(files) <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = ".sf|.sf.gz"))))
  if(length(files) > 0){
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
  tximportcounts <- as.data.frame(txi.salmon$counts)
  } else if(length(files) == 0){txtype <- 6}
} else if (txtype==2){
  if(getgeofiles == FALSE){
  dir <- readline(prompt=("Drop a directory containing Sailfish output into R Window or Enter Directory Path: "))
  }
files <- list.files(dir, recursive=TRUE, full.names = TRUE, pattern = "quant.sf|quant.sf.gz")
names(files) <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = "quant.sf|quant.sf.gz"))))
  if(length(files) > 0){
  txi.sailfish <- tximport(files, type = "sailfish", tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
  tximportcounts <- as.data.frame(txi.sailfish$counts)
  } else if(length(files) == 0){txtype <- 6}
} else if (txtype==3){
cat("Import Kallisto abundances from \"abundance.h5\" (requires rhdf5 library), or \"abundance.tsv*\".\n")
kallistotype <- readline(prompt=("Select kallisto datatype by entering \"h5\" or \"tsv\" (without quotes): "))
  if (kallistotype=="h5"){
  library("rhdf5")
    if(getgeofiles == FALSE){
    dir <- readline(prompt=("Drop a directory containing Kallisto output into R Window or Enter Directory Path: "))
    }
  files <- list.files(dir, recursive=TRUE, full.names = TRUE, pattern = "abundance.h5|abundance.h5.gz")
  if(length(files) > 0){
    names(files) <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = "abundance.h5|abundance.h5.gz"))))
    txi.kallisto.h5 <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
    tximportcounts <- as.data.frame(txi.kallisto.h5$abundance)
    } else if(length(files) == 0){txtype <- 6}}
  if (kallistotype=="tsv"){
    if(getgeofiles == FALSE){
    dir <- readline(prompt=("Drop a directory containing Kallisto output into R Window or Enter Directory Path: "))
    }
  files <- list.files(dir, recursive=TRUE, full.names = TRUE, pattern = "abundance.tsv|abundance.tsv.gz")
    if(length(files) > 0){
    names(files) <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = "abundance.tsv|abundance.tsv.gz"))))
    txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
    tximportcounts <- as.data.frame(txi.kallisto.tsv$abundance)
    } else if(length(files) == 0){txtype <- 6}}
} else if (txtype==5){
  cat("Stringtie suppot not yet implemented\n")
#  cat("We can import files produced by running the \"stringtie -eB -G transcripts.gff <source_file.bam>\" command per the TXImport tutorial.\n")
#  dir <- readline(prompt=("Drop a directory containing Stringtie output folder into R Window or Enter Directory Path: "))
#  samples <- list.files(dir)
#  files <- file.path(dir, samples)
#  tmp <- read_tsv(files[1])
#  print(head(tmp))
#  strtx <- readline(prompt=("Enter the column name that contains the stringtie transcript identifiers: "))
#  strgene <- readline(prompt=("Enter the column name that contains gene identifiers for which you have a GSEA compatible chip file (or gene symbols): "))
#  tx2gene <- tmp[, c(strtx, strgene)]
#  print(head(tx2gene))
#  readline(prompt=("If this looks correct, press any key to continue..."))
#  txi <- tximport(tmp, type = "stringtie", tx2gene = tx2gene)
  FAIL <- TRUE
}
if (txtype==6){
  cat("Failover Mode. This method isn't really supported. Use at your own risk.\n")
  cat("Could not detect transcript quantifications with tximport.\n")
  cat(" We're going to have to make a lot of guesses here, but we'll get through this together.\n")
  if(getgeofiles == FALSE){
  txmatrix <- readline(prompt=("Drop Transcript Expression Matrix into R Window or Enter File Path (supports csv or tab delimited txt): "))
  library("tools")
  txmatrixtype <- file_ext(gsub(".gz$", "", txmatrix))
    if(txmatrixtype == "csv"){
    full <- read.table(txmatrix, sep=",", header=T, row.names=NULL)
      if(colnames(full)[1] == "X"){colnames(full)[1] <- colnames(full["ID"])}
      if(colnames(full)[1] == "row.names"){colnames(full)[1] <- colnames(full["ID"])}}
    if(txmatrixtype == "txt" | txmatrix == "tabular" | txmatrix == "tsv" | txmatrix == "tab"){
    full <- read.table(txmatrix, sep="\t", header=T, row.names=NULL)
      if(colnames(full)[1] == "X"){colnames(full)[1] <- colnames(full["ID"])}
      if(colnames(full)[1] == "row.names"){colnames(full)[1] <- "ID"}}
  cat("Transcript expression Matrix Imported\n")
  txlevel <- FALSE
  cat("\n")
  print(head(full,3))
  cat("\n")
  cat("When you're prompted for a CHIP file, instead provide a table mapping Transcript IDs to Gene Symbols and Descriptions USING CHIP HEADERS. Good Luck.\n")
  readline(prompt=("Press any key to continue..."))}
  if(getgeofiles == TRUE){
        cat("Attempting to parse multiple individually quantified samples into single matrix...\n")
        inputsamples <- list.files(outfile)
        inputsamplestxt <- list.files(outfile, pattern=".txt|.tsv|.tabular|.tab")
        txtsamplenumber <- length(inputsamplestxt)
        inputsamplescsv <- list.files(outfile, pattern=".csv")
        csvsamplenumber <- length(inputsamplescsv)
        totalsamplenumber <- txtsamplenumber + csvsamplenumber
        if(txtsamplenumber > 0){
          print(inputsamplestxt)
          cat(paste(txtsamplenumber),"tabular formatted samples were detected in input directory.\n")
          txtsamplepaths <- file.path(outfile, inputsamplestxt)
          names(txtsamplepaths) <- paste0(inputsamplestxt)
          import_txt <- lapply(txtsamplepaths, read.table, header = TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
          for(i in 1:txtsamplenumber){
            do
            colnames(import_txt[[i]]) <- paste(inputsamplestxt[i], colnames(import_txt[[i]]), sep = "_")
          }}
        if(csvsamplenumber > 0){
          print(inputsamplescsv)
          cat(paste(csvsamplenumber),"csv formatted samples were detected in input directory.\n")
          csvsamplepaths <- file.path(outfile, inputsamplescsv)
          names(csvsamplepaths) <- paste0(inputsamplescsv)
          import_csv <- lapply(csvsamplepaths, read.table, header = TRUE, row.names=NULL, sep=",", stringsAsFactors=FALSE)
          for(i in 1:csvsamplenumber){
            do
            colnames(import_csv[[i]]) <- paste(inputsamplestxt[i], colnames(import_csv[[i]]), sep = "_")
          }}
        if((txtsamplenumber > 0) == (csvsamplenumber > 0)){
          import.list <- append(import_txt, import_csv) #combine csv and txt imports
          cat(paste(length(import.list)),"total samples were detected in input directory.\n")
          } else if(csvsamplenumber == 0){
            import.list <- import_txt  #process txt only
            } else if(txtsamplenumber == 0){
              import.list <- import_csv  #process csv only
              }
        cat("We need to make sure that all quantifications come from the same analysis pipeline.\n")
        keep <- as.numeric(readline(prompt=("HOW MANY of the extracted quantifications do you want to keep for analysis? ")))
        if(keep < length(names(import.list))){
        list <- c(1:length(names(import.list)))
        for(i in 1:keep){
          do
            cat("\n")
            cat("Sample: ",paste0(names(import.list[i])),"\n")
            value <- askYesNo("Include in analysis? ")
            list[i] <- value}
            import.list <- import.list[list]}

        importdata <- as.data.frame(colnames(import.list[[1]]), stringsAsFactors=FALSE, header=TRUE)
        colnames(importdata)[1] <- "EXPERIMENT"
        cat("\n")
        print(importdata)
        cat("\n")
        cat("Selected a sample to use for learning dataset format: \n")
        expids <- readline(prompt=("Enter the name or row number from the EXPERIMENT column above that defines your transcript identifiers: "))
        expidnumber <- match(expids, cbind(rownames(importdata),importdata)[,1])
          if(is.na(expidnumber) == TRUE){
          expidnumber <- match(expids, cbind(rownames(importdata),importdata)[,2])}
        mergedexp <- Reduce(function(x, y) merge(x, y, all=FALSE, by= expidnumber , all.x=TRUE, all.y=TRUE), import.list,accumulate=F)
        colnames(mergedexp)[expidnumber] <- "txID"

        mergedexp %>% select_if(is.numeric) -> mergedexp_2
        mergedexp_2 <- cbind(mergedexp[expidnumber], mergedexp_2)
        mergedexp_2 <- mergedexp_2[, -grep("ength$", colnames(mergedexp_2))]
        cat("Cleaned up identifiable extraneous non-numeric columns.\n")

        fullimportdata <- as.data.frame(colnames(mergedexp_2), stringsAsFactors=FALSE, header=TRUE)
        colnames(fullimportdata)[1] <- "EXPERIMENT"
        print(fullimportdata)
        cat("There are",paste(length(colnames(mergedexp_2))),"entries in the Experiment data matrix.\n")
        keepcols <- readline(prompt=("How many are tx quantifications? (this should be one per sample): "))
        removecols <- as.numeric(length(colnames(mergedexp_2))) - (1 + as.numeric(keepcols))
        if(removecols > 0){
          repeat{
        importdataloop <- fullimportdata
        colnames(importdataloop) <- c("EXPERIMENT")
        fullloop <- mergedexp_2
        print(importdataloop)
          for (i in 1:(as.numeric(removecols))){
            do
              cat("Discard everything except your Transcript IDs and your per sample quantifications. \n")
              expidsdrop <- readline(prompt=("Enter a name or row number from the EXPERIMENT column above that defines fields you want to DISCARD: "))
                expiddropnumber <- match(expidsdrop, cbind(rownames(importdataloop),importdataloop)[,1])
                if(is.na(expiddropnumber) == TRUE){expiddropnumber <- FALSE}
                if(expiddropnumber == expidsdrop){
                expidsdrop <- importdataloop[expiddropnumber,]}
                if(length(expidsdrop) == "1"){
                  cat("Dropping unused identifiers\n")
                  cat("\n")
                  fullloop <- fullloop[ , -which(names(fullloop) %in% c(expidsdrop))]
                  importdataloop <- as.data.frame(colnames(fullloop), stringsAsFactors=FALSE, header=TRUE)
                  colnames(importdataloop) <- c("EXPERIMENT")
                  print(importdataloop)
                  cat("\n")
                  expidsdrop <- NULL}
          }
        check <- askYesNo("Does this display only a single transcript identifier and the sample ids?")
        if(check == TRUE){
        coldata <- importdataloop
        tximportcounts <- fullloop
        expids <- expidnumber
        print(head(tximportcounts[expidnumber],3))
        break
        }}
      } else if (removecols == 0){
        coldata <- fullimportdata
        tximportcounts <- mergedexp_2
        expids <- expidnumber
        print(head(tximportcounts[expidnumber],3))
        }
      txhasversions <- askYesNo("Do your TRANCRIPT IDs have decimal versions? (eg. ENST00000342771.9)? ")
        if(txhasversions == TRUE){
            tximportcounts[,1] <- gsub("\\..*","",tximportcounts[,1])
            txhasversions <- FALSE}

        cat("\n")
        cat("Expression Matrix Imported\n")
        cat("\n")
        print(head(tximportcounts[expidnumber],3))
        cat("\n")
        cat("When you're prompted for a CHIP file, instead provide a table mapping Transcript IDs to Gene Symbols and Descriptions USING CHIP HEADERS. Good Luck.\n")
        readline(prompt=("Press any key to continue..."))
}
  }}


if (txtype==4){
TYPE <- "RSEM"
  if(getgeofiles == FALSE){
  dir <- readline(prompt=("Drop a directory containing RSEM output into R Window or Enter Directory Path: "))
  }
rsemcontents <- as.data.frame(paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = ".results.gz")))), stringsAsFactors=FALSE)
rsemgenes <- as.data.frame(paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = ".genes.results.gz")))), stringsAsFactors=FALSE)
rsemtxs <- as.data.frame(paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = ".isoforms.results.gz")))), stringsAsFactors=FALSE)
rsemgenetest <- rsemgenes[,1] %in% rsemcontents[,1]
rsemtxtest <- rsemtxs[,1] %in% rsemcontents[,1]
  if(all(rsemgenetest == TRUE) == TRUE){
  cat("Gene level RSEM quantifications are available. Using these directly.\n ")
  USEDRSEMTXLEVEL <- FALSE
  files <- list.files(dir, recursive=TRUE, full.names = TRUE, pattern = ".genes.results.gz")
  names(files) <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = ".genes.results.gz"))))
  txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  txi.rsemcounts <- as.data.frame(txi.rsem$counts)
  } else if(all(rsemtxtest == TRUE) == TRUE){
  cat("Gene level RSEM quantifications are NOT available.\n ")
  cat("Using RSEM isoform abundances. This is a little messy.\n ")
  USEDRSEMTXLEVEL <- TRUE
  files <- list.files(dir, recursive=TRUE, full.names = TRUE, pattern = ".isoforms.results.gz")
  names(files) <- paste0(tools::file_path_sans_ext(tools::file_path_sans_ext(list.files(dir, recursive=TRUE, full.names = FALSE, pattern = ".isoforms.results.gz"))))
  txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
  txi.rsemcounts <- as.data.frame(txi.rsem$counts)}
print("Done") #fixes a small bug in tximport rsem output
cat("\n")
cat("Identifiers:\n")
print(head(txi.rsemcounts[,1],3))
txi.rsemcounts <- tibble::rownames_to_column(txi.rsemcounts, "geneID")
cat("\n")
if(all(rsemgenetest == TRUE) == TRUE){genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)?")
} else if(USEDRSEMTXLEVEL == TRUE){genehasversions <- askYesNo("Do your TRANSCRIPT IDs have decimal versions? (eg. ENST00000342771.9)?")}
cat("\n")
if(genehasversions == TRUE){
txi.rsemcounts[,1] <- gsub("\\..*","",txi.rsemcounts[,1])
genehasversions <- FALSE}
txi.rsemcounts <- distinct(txi.rsemcounts)
txi.rsemcounts %>%
  group_by(geneID) %>%
  summarise_all(sum) %>%
  data.frame() -> txi.rsemcounts
tximportcounts<- txi.rsemcounts
if(USEDRSEMTXLEVEL == TRUE){colnames(tximportcounts)[1] <- "txID"
  cat("When you're prompted for a CHIP file, instead provide a table mapping Transcript IDs to Gene Symbols and Descriptions USING CHIP HEADERS. Good Luck.\n")
  readline(prompt=("Press any key to continue..."))}
}

if(txlevel == FALSE){
iscounts <- askYesNo("Is your dataset gene level COUNTS measurements from HTSeq-counts, FeatureCounts, or similar?")
  if(iscounts == TRUE){
  if(isnormalized == FALSE){
  isnormalized <- askYesNo("Are your counts already normalized?")
    if(isnormalized == TRUE){
    NORM <- TRUE
    TYPE <- "NORMCOUNTS"
    DESEQ2DONE <- FALSE}
    if(isnormalized == FALSE){
    NORM <- FALSE
    TYPE <- "RAWCOUNTS"}
    } else if (isnormalized == TRUE){
      NORM <- TRUE
      TYPE <- "NORMCOUNTS"
      DESEQ2DONE <- FALSE}}
else if (iscounts == FALSE){
cat("Your dataset is not supported at this time.\n")
FAIL <- TRUE
}}

if (FAIL == TRUE) {stop("An unsupported datatype was encountered and processing was terminated.")}

#Import Gene Expression Matrix ##GEO conditional Needed for genematrix prompt
if(txlevel == FALSE){
  if(getgeofiles == FALSE){
    genematrix <- readline(prompt=("Drop Gene Expression Matrix into R Window or Enter File Path (supports csv or tab delimited txt): "))
    }
  library("tools")
  genematrixtype <- file_ext(gsub(".gz$", "", genematrix))
    if(genematrixtype == "csv"){
    full <- read.table(genematrix, sep=",", header=T, row.names=NULL)
      if(colnames(full)[1] == "X"){colnames(full)[1] <- "ID"}
      if(colnames(full)[1] == "row.names"){colnames(full)[1] <- "ID"}
      cat("\n")
      cat("Expression Matrix Imported\n")
      cat("\n")
      print(head(full,3))
      cat("\n")}
    if(genematrixtype == "txt" | genematrixtype == "tabular" | genematrixtype == "tsv" | genematrixtype == "tab"){
    full <- read.table(genematrix, sep="\t", header=T, row.names=NULL)
      if(colnames(full)[1] == "X"){colnames(full)[1] <- "ID"}
      if(colnames(full)[1] == "row.names"){colnames(full)[1] <- "ID"}
      cat("\n")
      cat("Expression Matrix Imported\n")
      cat("\n")
      print(head(full,3))
      cat("\n")}
    if((genematrixtype == "") == TRUE){
      is_directory <- TRUE
        cat("Directory detected. This directory must contain ONLY samples to be imported\n")
        cat("Attempting to parse multiple individually quantified samples into single matrix...\n")
        inputsamples <- list.files(genematrix)
        inputsamplestxt <- list.files(genematrix, pattern=".txt|.tsv|.tabular|.tab")
        txtsamplenumber <- length(inputsamplestxt)
        inputsamplescsv <- list.files(genematrix, pattern=".csv")
        csvsamplenumber <- length(inputsamplescsv)
        totalsamplenumber <- txtsamplenumber + csvsamplenumber
        if(txtsamplenumber > 0){
          print(inputsamplestxt)
          cat(paste(txtsamplenumber),"tabular formatted samples were detected in input directory.\n")
          txtsamplepaths <- file.path(genematrix, inputsamplestxt)
          names(txtsamplepaths) <- paste0(inputsamplestxt)
          import_txt <- lapply(txtsamplepaths, read.table, header = TRUE, row.names=NULL, sep="\t", stringsAsFactors=FALSE)
          for(i in 1:txtsamplenumber){
            do
            colnames(import_txt[[i]]) <- paste(inputsamplestxt[i], colnames(import_txt[[i]]), sep = "_")
          }}
        if(csvsamplenumber > 0){
          print(inputsamplescsv)
          cat(paste(csvsamplenumber),"csv formatted samples were detected in input directory.\n")
          csvsamplepaths <- file.path(genematrix, inputsamplescsv)
          names(csvsamplepaths) <- paste0(inputsamplescsv)
          import_csv <- lapply(csvsamplepaths, read.table, header = TRUE, row.names=NULL, sep=",", stringsAsFactors=FALSE)
          for(i in 1:csvsamplenumber){
            do
            colnames(import_csv[[i]]) <- paste(inputsamplestxt[i], colnames(import_csv[[i]]), sep = "_")
          }}
        if((txtsamplenumber > 0) == (csvsamplenumber > 0)){
          import.list <- append(import_txt, import_csv) #combine csv and txt imports
          cat(paste(length(import.list)),"total samples were detected in input directory.\n")
          } else if(csvsamplenumber == 0){
            import.list <- import_txt  #process txt only
            } else if(txtsamplenumber == 0){
              import.list <- import_csv  #process csv only
              }
        cat("We need to make sure that all quantifications come from the same analysis pipeline.\n")
        keep <- as.numeric(readline(prompt=("HOW MANY of the extracted quantifications do you want to keep for analysis? ")))
        if(keep < length(names(import.list))){
        list <- c(1:length(names(import.list)))
        for(i in 1:keep){
          do
            cat("\n")
            cat("Sample: ",paste0(names(import.list[i])),"\n")
            value <- askYesNo("Include in analysis? ")
            list[i] <- value}
            import.list <- import.list[list]
          }
        #importlist<-
        #import.list <-
        importdata <- as.data.frame(colnames(import.list[[1]]), stringsAsFactors=FALSE, header=TRUE)
        colnames(importdata)[1] <- "EXPERIMENT"
        cat("\n")
        print(importdata)
        cat("\n")
        cat("Selected a sample to use for learning dataset format: \n")
        expids <- readline(prompt=("Enter the name or row number from the EXPERIMENT column above that defines your gene identifiers: "))
        expidnumber <- match(expids, cbind(rownames(importdata),importdata)[,1])
          if(is.na(expidnumber) == TRUE){
          expidnumber <- match(expids, cbind(rownames(importdata),importdata)[,2])}
        mergedexp <- Reduce(function(x, y) merge(x, y, all=FALSE, by= expidnumber , all.x=TRUE, all.y=TRUE), import.list,accumulate=F)
        colnames(mergedexp)[expidnumber] <- "geneID"

        mergedexp %>% select_if(is.numeric) -> mergedexp_2
        mergedexp_2 <- cbind(mergedexp[expidnumber], mergedexp_2)
        mergedexp_2 <- mergedexp_2[, -grep("ength$", colnames(mergedexp_2))]
        cat("Cleaned up identifiable extraneous columns.\n")

        fullimportdata <- as.data.frame(colnames(mergedexp_2), stringsAsFactors=FALSE, header=TRUE)
        colnames(fullimportdata)[1] <- "EXPERIMENT"
        print(fullimportdata)
        cat("There are",paste(length(colnames(mergedexp_2))),"entries in the Experiment data matrix.\n")
        keepcols <- readline(prompt=("How many are gene quantifications? (this should be one per sample): "))
        removecols <- as.numeric(length(colnames(mergedexp_2))) - (1 + as.numeric(keepcols))
        if(removecols > 0){
          repeat{
        importdataloop <- fullimportdata
        colnames(importdataloop) <- c("EXPERIMENT")
        fullloop <- mergedexp_2
        print(importdataloop)
          for (i in 1:(as.numeric(removecols))){
            do
              cat("Discard everything except your Gene IDs and your per sample quantifications. \n")
              expidsdrop <- readline(prompt=("Enter a name or row number from the EXPERIMENT column above that defines fields you want to DISCARD: "))
                expiddropnumber <- match(expidsdrop, cbind(rownames(importdataloop),importdataloop)[,1])
                if(is.na(expiddropnumber) == TRUE){expiddropnumber <- FALSE}
                if(expiddropnumber == expidsdrop){
                expidsdrop <- importdataloop[expiddropnumber,]}
                if(length(expidsdrop) == "1"){
                  cat("Dropping unused identifiers\n")
                  cat("\n")
                  fullloop <- fullloop[ , -which(names(fullloop) %in% c(expidsdrop))]
                  importdataloop <- as.data.frame(colnames(fullloop), stringsAsFactors=FALSE, header=TRUE)
                  colnames(importdataloop) <- c("EXPERIMENT")
                  print(importdataloop)
                  cat("\n")
                  expidsdrop <- NULL}
          }
        check <- askYesNo("Does this display only a single gene identifier and the sample ids?")
        if(check == TRUE){
        coldata <- importdataloop
        full2 <- fullloop
        expids <- expidnumber
        print(head(full2[expidnumber],3))
        break
        }}
      } else if (removecols == 0){
        coldata <- fullimportdata
        full2 <- mergedexp_2
        expids <- expidnumber
        print(head(full2[expidnumber],3))
        }
      genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)? ")
        if(genehasversions == TRUE){
            full2[,1] <- gsub("\\..*","",full2[,1])
            genehasversions <- FALSE}

        cat("\n")
        cat("Expression Matrix Imported\n")
        cat("\n")
        print(head(full2[expidnumber],3))
        cat("\n")
      }

} else if(txlevel == TRUE){
    full <- tximportcounts
    if(txtype != 6){
    cat("Using TXImport result...\n")}
    coldata <- as.data.frame(colnames(full), stringsAsFactors=FALSE, header=TRUE)
    coldata <- rbind(c("geneID"), coldata)
    colnames(coldata) <- c("EXPERIMENT")
    expids <- 0
    full2 <- full
    cat("\n")
    print(head(full2,3))
    cat("\n")}

if ((txlevel == FALSE | TYPE == "RSEM") == (is_directory == FALSE)){
cat("There are",paste(length(colnames(full))),"columns in your data table.\n")
samplesize <- readline(prompt=("How many of these are sequenced SAMPLES? "))
geneids <-  (length(colnames(full)) - as.numeric(samplesize))
repeat{
if(geneids == "1"){
#Prompt User for Original Experiment Namespace to Merge with CHIP On
cat("Gene Expression File Header:\n")
cat("\n")
coldata <- as.data.frame(colnames(full), stringsAsFactors=FALSE, header=TRUE)
colnames(coldata) <- c("EXPERIMENT")
print(as.data.frame(coldata))
cat("\n")
expids <- readline(prompt=("Enter the name or row number from the EXPERIMENT column above that defines your gene identifiers: "))
expidnumber <- match(expids, cbind(rownames(coldata),coldata)[,1])
  if(is.na(expidnumber) == TRUE){expidnumber <- FALSE}
  if(expidnumber == expids){
  expids <- coldata[expidnumber,]}
# FIX ENSG_SYMBOL MERGE
  if (txlevel == FALSE){
  if(all((grepl("_", full[,expids], fixed=TRUE))) == TRUE){
  library(tidyr)
  full <- separate(full, expids, c("geneID", NA), sep = "_", remove = TRUE, extra = "warn")
  expids <- "geneID"
  }}
full2 <- full
cat("Gene IDs:\n")
print(head(full2[,expids],3))
  if (txlevel == FALSE){genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)?")
    if(genehasversions == TRUE){
    full2[,1] <- gsub("\\..*","",full2[,1])
    genehasversions <- FALSE}}
coldata <- as.data.frame(colnames(full2), stringsAsFactors=FALSE, header=TRUE)
colnames(coldata) <- c("EXPERIMENT")
break }
else if(geneids > "1"){
fullloop <- full
#Prompt User for Original Experiment Namespace to Merge with CHIP On
cat("Gene Expression File Header:\n")
cat("\n")
coldata <- as.data.frame(colnames(full), stringsAsFactors=FALSE, header=TRUE)
colnames(coldata) <- c("EXPERIMENT")
print(as.data.frame(coldata))
cat("\n")
cat("We now need to pick one set of identifiers to use for downstream processing, and prune away the others one at a time.\n")
cat("For example, if you have both ENSEMBL Gene IDs and Gene Symbols columns KEEP the ENSEMBL IDs and DISCARD the Gene Symbols\n")
cat("\n")
expids <- readline(prompt=("Enter the name or row number from the EXPERIMENT column above that defines the gene identifiers you want to KEEP: "))
expidnumber <- match(expids, cbind(rownames(coldata),coldata)[,1])
  if(is.na(expidnumber) == TRUE){expidnumber <- FALSE}
  if(expidnumber == expids){
  expids <- coldata[expidnumber,]}
coldataloop <- coldata
for (i in 1:(as.numeric(geneids) - 1)){
  do
    expidsdrop <- readline(prompt=("Enter a name or row number from the EXPERIMENT column above that defines gene identifiers you want to DISCARD: "))
      expiddropnumber <- match(expidsdrop, cbind(rownames(coldataloop),coldataloop)[,1])
      if(is.na(expiddropnumber) == TRUE){expiddropnumber <- FALSE}
      if(expiddropnumber == expidsdrop){
      expidsdrop <- coldataloop[expiddropnumber,]}
    if(length(expidsdrop) == "1"){
    cat("Dropping unused identifiers\n")
    cat("\n")
    fullloop <- fullloop[ , -which(names(fullloop) %in% c(expidsdrop))]
    coldataloop <- as.data.frame(colnames(fullloop), stringsAsFactors=FALSE, header=TRUE)
    colnames(coldataloop) <- c("EXPERIMENT")
    print(coldataloop)
    cat("\n")
    expidsdrop <- NULL}}
#check <- askYesNo("Does this look right?")}
#if(check == TRUE){
#
#
#}
}
check <- askYesNo("Does this display only a single gene identifier and the sample ids?")
if(check == TRUE){
coldata <- coldataloop
full2 <- fullloop
print(head(full2,3))
if(TYPE == "NORMCOUNTS" | TYPE == "RAWCOUNTS" | TYPE == "TXABUNDANCE"){
genehasversions <- askYesNo("Do your GENE IDs have decimal versions? (eg. ENSG00000158321.16)?")
} else if(TYPE == "RSEM"){genehasversions<- FALSE}
if(genehasversions == TRUE){
full2[,1] <- gsub("\\..*","",full2[,1])
genehasversions <- FALSE}
break
}
}
}

#Import CHIP File for Processing
CHIPpath <- readline(prompt=("Drop appropriate CHIP File matching the namespace of identifiers into R Window or Enter File Path: "))
rawchip <- read.table(CHIPpath, sep="\t", comment.char = "", quote="", stringsAsFactors=FALSE, fill = TRUE, header=T)
chip <- rawchip[c("Probe.Set.ID","Gene.Symbol")]
colnames(chip) <- c("Raw_IDs", "NAME")
fullchip <- unique(rawchip[c("Gene.Symbol","Gene.Title")])
colnames(fullchip) <- c("NAME", "Description")

cat("Done\n")
cat("\n")

if(median(nchar(colnames(full2))) > 25){
cat("Your sample names are pretty long.\n")
replacenames <- askYesNo("Would you like to replace them with something friendlier?")
  if(replacenames == TRUE){
    replacenameslength <- length(colnames(full2))
    for (i in 1:(replacenameslength - 1)){
      do
      print(colnames(full2[i +1]))
      newname <- readline(prompt=("What would you like to call this sample? "))
      colnames(full2)[i + 1] <- newname
      }
    coldata <- as.data.frame(colnames(full2), stringsAsFactors=FALSE, header=TRUE)
    colnames(coldata) <- c("EXPERIMENT")
    cat("\n")
    cat("Formatted experiment:\n")
    print(coldata)
    cat("\n")
  }}

#Prompt User for DESeq2 formatted "coldata" Experiment Design File
cat("We can build the GSEA CLS file from a DESeq2 \"coldata\" file.\n")
colavailable <- askYesNo("Do you have an existing coldata file?")
if(colavailable == TRUE){
design <- readline(prompt=("Drop DESeq2 formatted coldata File into R Window or Enter File Path, see DESeq2 Tutorial for Formatting (https://bioconductor.org/packages/3.8/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#input-data): "))
coldata <- read.table(design, sep="\t", header=T, row.names =1)
colnames(coldata) <- c("condition")
} else if(colavailable == FALSE){
cat("Okay. We'll build annotations from scratch...\n")
#construct coldata directly
coldata <- cbind(coldata,c(""), stringsAsFactors=FALSE)
idrow <- which(grepl(expids, coldata))
coldata <- as.data.frame(coldata[-c(idrow), ], stringsAsFactors=FALSE)
rownames(coldata) <- 1:nrow(coldata)
colnames(coldata) <- c("name", "condition")
looplen <- length(rownames(coldata))
#levels(coldata$condition)= c("EX", "CT", "ex", "ct", "Experimental", "Control", "experimental", "control")
cat("\n")
cat("We now need to generate an experimental labels file to build the GSEA CLS file.\n")
cat("\n")
phnumber <- readline(prompt=("How many phenotypes are in your experiment? (numbers only) "))
Phenotypes <- 0
for (i in 1:phnumber){
  do
cat("What single word do you want to name phenotype",paste0(i),"? ")
Phenotypes[i] <- readline()
}
cat("\n")
print(as.data.frame(Phenotypes))
cat("\n")
Phenotypeframe <- as.data.frame(Phenotypes, stringsAsFactors=FALSE)
repeat{
for (i in 1:looplen){
  do
  cat("Select a phenotype for sample >",paste(coldata[i,1]),": ")
  value <- readline()
  Phenotypesnumber <- match(value, cbind(rownames(Phenotypeframe),Phenotypeframe)[,1])
  if(is.na(Phenotypesnumber) == TRUE){coldata[i,2] <- value
  } else if(Phenotypesnumber == value){
  coldata[i,2] <- Phenotypeframe[value,]}
}

cat("\n")
print(as.data.frame(coldata))
cat("\n")
cat("Make absolutely sure the conditions are listed correctly before continuing!\n")
cat("We detected",paste(length(unique(coldata$condition))),"phenotypes:",paste(unique(coldata$condition)),"\n")
check2 <- askYesNo("Does the sample to phenotype assignment apear correct (make absolutely sure before continuing)?")
if(check2 == TRUE){
rownames(coldata) <- coldata$"name"
coldata <- coldata[c(2)]
break
}}}

#Merge Gene Expression Matrix with CHIP File and Process Identifiers
mappedexp <- merge(x = chip, y= full2, by.x="Raw_IDs", by.y=expids, all=FALSE)
size <- length(colnames(mappedexp))
mappedexp <- mappedexp[c(2:size)]

cat("Summing Counts that mapped to the same gene after applying mapping chip\n")

#SUM Counts for Identifiers Mapping to the Same Gene
mappedexp <- distinct(mappedexp)
mappedexp %>%
  group_by(NAME) %>%
  summarise_all(sum) %>%
  data.frame() -> mappedexp_sum

#Set Gene Names as Index Column
mappedexp_sum2 <- mappedexp_sum[,-1]
rownames(mappedexp_sum2) <- mappedexp_sum[,1]
rownames(coldata) <- colnames(mappedexp_sum2)
cat("Done\n")
outprefix <- readline(prompt=("Enter a prefix to label output files: "))
cat("\n")

if(NORM == FALSE){
cat("Your dataset is flagged as needing normalization to be compatible with GSEA.\n")
donormalize <- askYesNo("We can now normalize your dataset with DESEq2. Continue?")
if(donormalize == FALSE){
  cat("Perofrming standard filtering of low count genes without additional normalization.\n")
  cat("You probably don't want to do this. We don't think this dataset is properly normalized.\n")
  mappedexp_sum3 <- subset(mappedexp_sum2, rowSums(mappedexp_sum2[])>=5 )
  protoGCT <- merge(x = fullchip, y= mappedexp_sum3, by.x="NAME", by.y=0, all=FALSE)
  bound <- rbind(colnames(protoGCT), protoGCT)
  bound <- rbind(NA, bound)
  bound <- rbind(NA, bound)
  bound[1,1] <- "#1.2"
  numberofsamples <- length(colnames(bound))-2
  numberofgenes <- length(bound$NAME)-3
  bound[2,1] <- numberofgenes
  bound[2,2] <- numberofsamples
  cat("Writing final .GCT file for GSEA\n")
  write.table(bound, paste0(outprefix,"_Formatted.gct"), sep="\t", quote=F, row.names=FALSE, col.names=FALSE, na="")}
if(donormalize == TRUE) {
  cat("Loading DESEq2 Library...\n")
  library(DESeq2)
  cat("Begin DESeq2 Normalization...\n")
  mappedexp_sum2 <- round(mappedexp_sum2)
  dds <- DESeqDataSetFromMatrix(countData = mappedexp_sum2,
                                colData = coldata,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds)
  DESEQ2DONE <- TRUE
}}
if(DESEQ2DONE == TRUE){
cat("Writing Size Factor Normalizated GCT\n")
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized=TRUE)
norm <- tibble::rownames_to_column(as.data.frame(norm), "NAME")
protoGCT <- merge(x = fullchip, y= norm, by.x="NAME", by.y="NAME", all.y=TRUE)
bound <- rbind(colnames(protoGCT), protoGCT)
bound <- rbind(NA, bound)
bound <- rbind(NA, bound)
bound[1,1] <- "#1.2"
numberofsamples <- length(colnames(bound))-2
numberofgenes <- length(bound$NAME)-3
bound[2,1] <- numberofgenes
bound[2,2] <- numberofsamples
write.table(bound, paste0(outprefix,"_Default_Normalized_Counts.gct"), sep="\t", quote=F, row.names=FALSE, col.names=FALSE, na="")
}

if(NORM == TRUE){
cat("You've indicated that you dataset is already normalized.\n")
cat("Perofrming standard filtering of low count genes without additional normalization.\n")
mappedexp_sum3 <- subset(mappedexp_sum2, rowSums(mappedexp_sum2[])>=5 )
protoGCT <- merge(x = fullchip, y= mappedexp_sum3, by.x="NAME", by.y=0, all=FALSE)
bound <- rbind(colnames(protoGCT), protoGCT)
bound <- rbind(NA, bound)
bound <- rbind(NA, bound)
bound[1,1] <- "#1.2"
numberofsamples <- length(colnames(bound))-2
numberofgenes <- length(bound$NAME)-3
bound[2,1] <- numberofgenes
bound[2,2] <- numberofsamples
cat("Writing final .GCT file for GSEA\n")
write.table(bound, paste0(outprefix,"_Formatted.gct"), sep="\t", quote=F, row.names=FALSE, col.names=FALSE, na="")
}

if(altnorm == "usevst" | altnorm == "useall"){cat("Writing Variance Stabilizing Transformation Normalizated GCT\n")
vsd <- vst(dds, blind=FALSE)
vstnorm <- as.data.frame(assay(vsd))
vstnorm <- tibble::rownames_to_column(as.data.frame(vstnorm), "NAME")
vstprotoGCT <- merge(x = fullchip, y= vstnorm, by.x="NAME", by.y="NAME", all.y=TRUE)
vstbound <- rbind(colnames(vstprotoGCT), vstprotoGCT)
vstbound <- rbind(NA, vstbound)
vstbound <- rbind(NA, vstbound)
vstbound[1,1] <- "#1.2"
vstnumberofsamples <- length(colnames(vstbound))-2
vstnumberofgenes <- length(vstbound$NAME)-3
vstbound[2,1] <- vstnumberofgenes
vstbound[2,2] <- vstnumberofsamples
write.table(vstbound, paste0(outprefix,"_vst_normalized_Counts.gct"), sep="\t", quote=F, row.names=FALSE, col.names=FALSE, na="")}

if(altnorm == "userlog" | altnorm == "useall"){cat("Writing RLOG Normalizated GCT\n")
rld <- rlog(dds, blind=FALSE)
rlognorm <- as.data.frame(assay(rld))
rlognorm <- tibble::rownames_to_column(as.data.frame(rlognorm), "NAME")
rlogprotoGCT <- merge(x = fullchip, y= rlognorm, by.x="NAME", by.y="NAME", all.y=TRUE)
rlogbound <- rbind(colnames(rlogprotoGCT), rlogprotoGCT)
rlogbound <- rbind(NA, rlogbound)
rlogbound <- rbind(NA, rlogbound)
rlogbound[1,1] <- "#1.2"
rlognumberofsamples <- length(colnames(rlogbound))-2
rlognumberofgenes <- length(rlogbound$NAME)-3
rlogbound[2,1] <- rlognumberofgenes
rlogbound[2,2] <- rlognumberofsamples
write.table(rlogbound, paste0(outprefix,"_rlog_normalized_Counts.gct"), sep="\t", quote=F, row.names=FALSE, col.names=FALSE, na="")}

uniqueclsclasses <- unique(coldata$condition)
clsclasses <- coldata$condition
clsspan <- length(rownames(coldata))
clslen <- length(uniqueclsclasses)
cls = data.frame(matrix(vector(), 3, clsspan,
                dimnames=list(c())),
                stringsAsFactors=F)
cls[1,1] <- clsspan
cls[1,2] <- clslen
cls[1,3] <- "1"
cls[2,1] <- "# "
cls[2,2:(1 + clslen)] <- uniqueclsclasses
cls[3,] <- clsclasses
cat("Writing phenotype definitions .CLS file for GSEA.\n")
write.table(cls, paste0(outprefix,"_phenotypes.cls"), sep="\t", quote=F, row.names=FALSE, col.names=FALSE, na="")
cat("All Done!\n")
