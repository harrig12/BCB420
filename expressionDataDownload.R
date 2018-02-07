#get GSE84712 data
#modified from expressionAnalysisSampleCode.R

if (! require(Biobase, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biobase")
  library(Biobase)
}

if (! require(GEOquery, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GEOquery")
  library(GEOquery)
}

# Work with RNAseq seems much easier by comparison: we already get unique gene
# Gene.IDs.

# download the normalized, annotated dataset from GEO
#
fName <- "GSE84712_TPM_all_samples.txt"
library(readr)
GSE84712 <- read_tsv(fName)
str(GSE84712)
save(GSE84712, file="/GSE84712.RData")

load("GSE84712.RData")
load("zu/inst/extdata/HUGOsymbols.RData")

nrow(GSE84712) # 19097 rows
any(is.na(GSE84712$Gene.ID))       # FALSE
any(GSE84712$Gene.ID == "")        # FALSE
any(grepl("/", GSE84712$Gene.ID))  # FALSE
any(duplicated(GSE84712$Gene.ID))  # FALSE

sum(HUGOsymbols %in% GSE84712$Gene.ID) # 17628: 92 %

missingRNAseqGene.IDs <- HUGOsymbols[!(HUGOsymbols %in% GSE84712$Gene.ID)]

# how many of these can be recovered, because they are aliases?

extraRNAseqGene.IDs <- GSE84712$Gene.ID[!(GSE84712$Gene.ID %in% HUGOsymbols)]
# what are these?
cat(sprintf("# [%d] %s\n", 1:20, head(extraRNAseqGene.IDs, 20)))
# [1] A1BG-AS1       -  antisense
# [2] A2MP1          -  pseudogene
# [3] AA06           -  ???
# [4] AAA1           -  synonym for NPSR1-AS1
# [5] AACSP1
# [6] AATK-AS1
# [7] ABCA11P
# [8] ABCA17P
# [9] ABCC13         - pseudogene
# [10] ABCC6P1
# [11] ABCC6P2
# [12] ABHD14A-ACY1  - readthrough (NMD candidate)
# [13] ABP1          - NEW SYMBOL: AOC1
# [14] ACN9          - NEW SYMBOL: SDHAF3
# [15] ACPL2         - NEW SYMBOL: PXYLP1
# [16] ACPT          - NEW SYMBOL: ACP4
# [17] ACRC          - NEW SYMBOL: GCNA
# [18] ACTR3BP2
# [19] ACTR3BP5
# [20] ADAM21P1

load("zu/inst/extdata/synMap.RData")
str(synMap)
# mappable:
x <- extraRNAseqGene.IDs[extraRNAseqGene.IDs %in% synMap$synonyms] # 1105

# recover 1105 Gene.IDs will raise coverage
# to 17628 + 1105 = 18733 ... 98 %

# value distributions:

boxplot(GSE84712[ , 2:7],
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE84712",
        outline = FALSE,
        col = c(rep(myArrayCols[5], 3), rep(myArrayCols[17], 3)))

# etc ...


# [END]


