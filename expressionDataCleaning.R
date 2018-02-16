# expressionDataCleaning.R
#
# Purpose: clean MappedGSE84712 data
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

myFunction <- function(a, b=1) {
	# Purpose:
	#     Describe ...
	# Parameters:
	#     a: ...
	#     b: ...
	# Value:
	#     result: ...

	# code ...

	return(result)
}



# ====  PROCESS  ===============================================================
# Enter the step-by-step process of your script here. Strive to write your
# code so that you can simply run this entire file and re-create all
# intermediate results.

load("MappedGSE84712.RData")

head(CGSE84712)

#find how many rows have little (or no) expression values recorded
(nrow(CGSE84712[CGSE84712$controlMean == 0 & CGSE84712$treatedMean == 0,])) # 737
#remove these (see journal for rational)
CGSE84712 <- CGSE84712[!(CGSE84712$controlMean == 0 & CGSE84712$treatedMean == 0),]

#have a peek
plot(CGSE84712$controlMean, CGSE84712$treatedMean,
     xlab = "control mean", ylab = "treated mean")

#pretty linear. What are thoes highly expressed genes?
(CGSE84712[CGSE84712$controlMean>10000,"Gene.ID"])
#tubulin and cytochrome C oxidase... not surprising

#plot again, without them (to get a better look at the rest)
plot(CGSE84712[CGSE84712$controlMean<10000,"controlMean"],
     CGSE84712[CGSE84712$controlMean<10000,"treatedMean"],
     xlab = "control mean", ylab = "treated mean")

# ====  TESTS  =================================================================
# Enter your tests here...  (Note: these are not tests for functions in the
# package!)


# [END]
