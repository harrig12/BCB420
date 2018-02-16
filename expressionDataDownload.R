#get GSE84712 data and map to HUGO symbols
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
# fName <- "GSE84712_TPM_all_samples.txt"
# library(readr)
# GSE84712 <- read_tsv(fName)
# str(GSE84712)
# save(GSE84712, file="/GSE84712.RData")

setwd("..")
load("GSE84712.RData")
load("zu/inst/extdata/HUGOsymbols.RData")
load("zu/inst/extdata/synMap.RData")

nrow(GSE84712) # 19097 rows
any(is.na(GSE84712$Gene.ID))       # FALSE
any(GSE84712$Gene.ID == "")        # FALSE
any(grepl("/", GSE84712$Gene.ID))  # FALSE
any(duplicated(GSE84712$Gene.ID))  # FALSE

#start by summarizing the data we're looking at
(colnames(GSE84712))
#controls are columns 3-29, experiment conditions are columns 30-81
GSE84712$controlMean <- rowMeans(GSE84712[,3:29])
GSE84712$treatedMean <- rowMeans(GSE84712[,30:81])
head(GSE84712[,c("controlMean", "treatedMean")], 20)

sum(HUGOsymbols %in% GSE84712$Gene.ID) # 17628: 92 %

missingRNAseqGene.IDs <- HUGOsymbols[!(HUGOsymbols %in% GSE84712$Gene.ID)]

# how many of these can be recovered, because they are aliases?
extraRNAseqGene.IDs <- GSE84712$Gene.ID[!(GSE84712$Gene.ID %in% HUGOsymbols)]

#synMap acts as a lookup table
str(synMap)

# mappable
mappable <- extraRNAseqGene.IDs[extraRNAseqGene.IDs %in% synMap$synonyms]
length(mappable) # 1105

# recovering 1105 Gene.IDs will raise coverage
# to 17628 + 1105 = 18733 ... 98 %

#which couldn't be mapped?
unmappable <- extraRNAseqGene.IDs[!(extraRNAseqGene.IDs %in% mappable)]
write(unmappable, file = "symUnmapped.txt", ncolumns = 1)

#update symMap with HUGO symbols that could be found
source("HUGOmapperParsing.R")
mappable <- extraRNAseqGene.IDs[extraRNAseqGene.IDs %in% synMap$synonyms]

#now that we have more mappings, there are some duplicates.
#(Two different synonyms mapping to a single HUGO symbol for example,
#probably due to the reference genome used for alignment)
wDuplicates <- GSE84712
wDuplicates$reMappedID <- wDuplicates$Gene.ID

#if a gene is mappable, rename it in wDuplicates
for (x in mappable){
  wDuplicates$reMappedID[wDuplicates$Gene.ID == x] <- synMap$symbols[synMap$synonyms == x]
}

#find which HUGO symbols are represented more than once in the experiment
#first find the duplicated geneIDs
dupHUGOs <- wDuplicates$reMappedID[which(duplicated(wDuplicates$reMappedID))]
(length(dupHUGOs)) #31

#then find original versions of those geneIDs
doubleMappedGenes <- wDuplicates[, c("Gene.ID","reMappedID", "controlMean", "treatedMean")][which(wDuplicates$reMappedID %in% doubleMapped),]
(nrow(doubleMappedGenes)) # 60, two of the doubleMapped genes were actually triple mapped

#Show original and new mappings, manually decide which expression values should be used
#as the values for that HUGO symbol. In class we decided we will keep only probe with
#highest mean/median with robust statistics. This translates similarly for RNAseq.
#print the average expression value beside each and choose which to keep or not.
(doubleMappedGenes)
doublesClean <- read.table("HUGOduplicateCleaning.txt", stringsAsFactors = F)
discardDouble <- doubleMappedGenes$Gene.ID[!doublesClean$retainClean]

#store the cleaned version of the experiment; all rows except those with
#the duplicates we eliminated
CGSE84712 <- GSE84712[!(GSE84712$Gene.ID %in% discardDouble),]
(nrow(CGSE84712)) #30 rows removed

#write mapped experiment to file
save(CGSE84712, file="MappedGSE84712.RData")

# value distributions:

cyclicPalette <- colorRampPalette(c("#00AAFF",
                                    "#DDDD00",
                                    "#FFAA00",
                                    "#00AAFF",
                                    "#DDDD00",
                                    "#FFAA00",
                                    "#00AAFF"))

boxplot(CGSE84712[ , 3:29],
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE84712 controls",
        outline = FALSE,
        col = cyclicPalette(27))

boxplot(CGSE84712[ , 30:55],
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE84712 Lead30",
        outline = FALSE,
        col = cyclicPalette(25))

boxplot(CGSE84712[ , 56:81],
        boxwex = 0.6,
        notch = TRUE,
        main = "GSE84712 Lead3",
        outline = FALSE,
        col = cyclicPalette(25))


# [END]


