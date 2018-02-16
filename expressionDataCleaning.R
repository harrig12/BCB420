# expressionDataCleaning.R
#
# Purpose: clean and normalize MappedGSE84712 data
# Version: 1
# Date: 02 2018
# Author: Cait Harrigan
#
# Input:MappedGSE84712.RData
# Output:CleanGSE84712
# Dependencies:
#
# Version history:
#
# ToDo:
# Notes:
#
# ==============================================================================


# ====  PARAMETERS  ============================================================
# Define and explain all parameters. No "magic numbers" in your code below.



# ====  PACKAGES  ==============================================================
# Load all required packages. Example for the idiom:

if (! require(stringr, quietly=TRUE)) {
  install.packages("stringr")
  library(stringr)
}
# Package information:
#  library(help = stringr)       # basic information
#  browseVignettes("stringr")    # available vignettes
#  data(package = "stringr")     # available datasets

if (! require(edgeR, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("edgeR")
  library(edgeR)
}

# ====  FUNCTIONS  =============================================================

# Define functions or source external files


# ====  PROCESS  ===============================================================
# Enter the step-by-step process of your script here. Strive to write your
# code so that you can simply run this entire file and re-create all
# intermediate results.

load("MappedGSE84712.RData")

head(CGSE84712)
CGSE84712 <- as.data.frame(CGSE84712)

#find how many rows have little (or no) expression values recorded
(nrow(CGSE84712[CGSE84712$controlMean == 0 & CGSE84712$treatedMean == 0,])) # 737

#may want to remove these (see journal for rational)
noLowGSE84712 <- CGSE84712[!(CGSE84712$controlMean == 0 & CGSE84712$treatedMean == 0),]
lowExpression <- CGSE84712[(CGSE84712$controlMean == 0 & CGSE84712$treatedMean == 0),]

#have a peek
plot(noLowGSE84712[,"controlMean"], noLowGSE84712[,"treatedMean"],
     xlab = "control mean", ylab = "treated mean")

#pretty linear. What are those highly expressed genes?
(noLowGSE84712[noLowGSE84712$controlMean>10000,"Gene.ID"])
#tubulin and cytochrome C oxidase... not surprising

#plot again, without them (to get a better look at the rest)
plot(noLowGSE84712[noLowGSE84712$controlMean<10000,"controlMean"],
     noLowGSE84712[noLowGSE84712$controlMean<10000,"treatedMean"],
     xlab = "control mean", ylab = "treated mean")

#only need the averages as factors for downstream
noLowGSE84712 <- noLowGSE84712[,c("Gene.ID", "controlMean", "treatedMean")]

#use edgeR for transcript normalization
noLowDGE <- DGEList(noLowGSE84712[, c("controlMean", "treatedMean")], group = factor(c(1,2)))
noLowDGE <- calcNormFactors(noLowDGE)

#consider the lows-included case
#incLowsDGE <- DGEList(CGSE84712[, c("controlMean", "treatedMean")], group = factor(c(1,2)))
#incLowsDGE <- calcNormFactors(incLowsDGE)
#look at all samples
#allDGE <- DGEList(CGSE84712[3:81])
#allDGE <- calcNormFactors(allDGE)

#scale with norm factors
noLowGSE84712$controlMean <- noLowGSE84712$controlMean*noLowDGE$samples$norm.factors[1]
noLowGSE84712$treatedMean <- noLowGSE84712$treatedMean*noLowDGE$samples$norm.factors[2]

save(noLowGSE84712, file="normGSE84712.RData")

# ====  TESTS  =================================================================
# Enter your tests here...  (Note: these are not tests for functions in the
# package!)


# [END]
