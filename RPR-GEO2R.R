# RPR_GEO2R.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR_GEO2R unit.
#
# Version:  1.1
#
# Date:     2017 09  -  2018 01
# Author:   Boris Steipe <boris.steipe@utoronto.ca>
#
# Versions:
#           1.1    Add section on GPL annotations
#           1.0    Updates for BCH441 2017
#           0.1    First code copied from 2016 material.
#
#
# TODO:
#
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...
#
# Note: to submit tasks for credit for this unit, report on the sections
# that have "Task ..." section headers, and report on the lines that are
# identified with #TASK>  comments.
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                Line
#TOC> --------------------------------------------------------------------
#TOC>   1        Preparations                                           53
#TOC>   2        Loading a GEO Dataset                                  84
#TOC>   3        Column wise analysis - time points                    154
#TOC>   3.1      Task - Comparison of experiments                      160
#TOC>   3.2      Grouped Samples                                       207
#TOC>   4        Row-wise Analysis: Expression Profiles                242
#TOC>   4.1      Task - Read a table of features                       277
#TOC>   4.2      Selected Expression profiles                          325
#TOC>   5        Differential Expression                               366
#TOC>   5.1      Final task: Gene descriptions                         490
#TOC>   6        Improving on Discovery by Differential Expression     495
#TOC>   7        Annotation data                                       577
#TOC>
#TOC> ==========================================================================


# =    1  Preparations  ========================================================

# To load and analyze GEO datasets we use a number of Bioconductor packages:


if (! require(Biobase, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biobase")
  library(Biobase)
}
# Package information:
#  library(help = Biobase)       # basic information
#  browseVignettes("Biobase")    # available vignettes
#  data(package = "Biobase")     # available datasets


if (! require(GEOquery, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GEOquery")
  library(GEOquery)
}
# Package information:
#  library(help = GEOquery)       # basic information
#  browseVignettes("GEOquery")    # available vignettes
#  data(package = "GEOquery")      # available datasets


# =    2  Loading a GEO Dataset  ===============================================

# The R code below is adapted from the GEO2R scripts produced by GEO for GSE3635
# for the experiment conducted in the BIN-EXPR-GEO unit

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Wed Jan 11 17:39:46 EST 2017


# Load series and platform data from GEO. The GEO server is a bit flakey and
# I have experienced outages over several hours. If the command below does
# not work for you, skip ahead to the fallback procedure.

GSE3635 <- getGEO("GSE3635", GSEMatrix =TRUE, getGPL=FALSE)
# Note: GEO2R scripts call the expression data set
#       "gset" throughout ... in this script I give
#       it the name "GSE3635" for clarity.

# Subset, if the data contains multiple experiments with different platforms
# since we should not be mixing data from different platforms. (That's
# unfortunately very hard to do correctly). The "platform" is the type of chip
# that was used, and what I selected online when this script was produced, so
# GEO2R knows about the correct platform ID. Technically, subsetting to just one
# platform is not necessary here because we _know_ that GSE3635 contains only
# data produced with the GPL1914 platform. But that's the way GEO2R scripts are
# setup be default since they have to be able to handle a variety of cases.
if (length(GSE3635) > 1) {
  idx <- grep("GPL1914", attr(GSE3635, "names"))
} else {
  idx <- 1
}

GSE3635 <- GSE3635[[idx]]

# FALLBACK
# ... in case the GEO server is not working, load the "GSE3635" object from
# the data directory:
# load(file="./data/GSE3635.RData")

# Checkpoint ...
if (! exists("GSE3635")) {
  stop("PANIC: GSE3635 was not loaded. Can't continue.")
}


# GSE3635 is an "Expression Set" - cf.
# https://bioconductor.org/packages/release/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf

# What does this contain?
help("ExpressionSet-class")

# Print it
GSE3635

# Access contents via methods:
featureNames(GSE3635)[1:20]   # Rows. These are the proteins being expressed. There are 6228.
sampleNames(GSE3635)[1:10]    # Columns. These are the samples at different timepoints

# Access contents by subsetting:
( tmp <- GSE3635[12:17, 1:6] )

# Access data
exprs(tmp)   # exprs() gives us the actual expression values.


#What are the data: in each cell?
#The expression value for a certain gene at a certain timepoint
#in each column?
#The expression values for all genes at a certain timepoint
#in each row?
#The expression values for a certain gene over all timepoints

# =    3  Column wise analysis - time points  ==================================

#  Each column represents one experiment.
#What are these experiments? These are the expression values found for different genes
#at a given timepoint, as measured by log ratios of present cDNA.


# ==   3.1  Task - Comparison of experiments  ==================================

# Get an overview of the distribution of data values in individual columns
summary(exprs(GSE3635)[ , 1])
summary(exprs(GSE3635)[ , 4])
summary(exprs(GSE3635)[ , 7])

# as a boxplot
cyclicPalette <- colorRampPalette(c("#00AAFF",
                                    "#DDDD00",
                                    "#FFAA00",
                                    "#00AAFF",
                                    "#DDDD00",
                                    "#FFAA00",
                                    "#00AAFF"))
tCols <- cyclicPalette(13)
boxplot(exprs(GSE3635), col = tCols)


#TASK>     Study this boxplot. What's going on? Are these expression values?
#These expression values are log ratios of the observed expression, relative to some chosen
#control column (cDAN of RNA from an asynchronous population of W303a cells)
#TASK>     What do the numbers that exprs() returns from the dataset mean? They are the relative
#expression level of the genes examined (rows) compared to the control (W303a) at a given
#timepoint (columns)


# Lets plot the distributions of values in a more fine-grained manner:
hT0  <- hist(exprs(GSE3635)[ ,  1], breaks = 100)
hT1  <- hist(exprs(GSE3635)[ ,  2], breaks = 100)
hT2  <- hist(exprs(GSE3635)[ ,  3], breaks = 100)
hT3  <- hist(exprs(GSE3635)[ ,  4], breaks = 100)
hT4  <- hist(exprs(GSE3635)[ ,  5], breaks = 100)
hT5  <- hist(exprs(GSE3635)[ ,  6], breaks = 100)
hT6  <- hist(exprs(GSE3635)[ ,  7], breaks = 100)
hT7  <- hist(exprs(GSE3635)[ ,  8], breaks = 100)
hT8  <- hist(exprs(GSE3635)[ ,  9], breaks = 100)
hT9  <- hist(exprs(GSE3635)[ , 10], breaks = 100)
hT10 <- hist(exprs(GSE3635)[ ,  11], breaks = 100)
hT11 <- hist(exprs(GSE3635)[ ,  12], breaks = 100)
hT12 <- hist(exprs(GSE3635)[ , 13], breaks = 100)


plot(  hT0$mids,  hT0$counts,  type = "l", col =  tCols[1],
       xlim = c(-0.5, 0.5), ylim = c(0,1700), main = "All Timepoints")
points(hT1$mids,  hT1$counts,  type = "l", col =  tCols[2])
points(hT2$mids,  hT2$counts,  type = "l", col =  tCols[3])
points(hT3$mids,  hT3$counts,  type = "l", col =  tCols[4])
points(hT4$mids,  hT4$counts,  type = "l", col =  tCols[5])
points(hT5$mids,  hT5$counts,  type = "l", col =  tCols[6])
points(hT6$mids,  hT6$counts,  type = "l", col =  tCols[7])
points(hT7$mids,  hT7$counts,  type = "l", col =  tCols[8])
points(hT8$mids,  hT8$counts,  type = "l", col =  tCols[9])
points(hT9$mids,  hT9$counts,  type = "l", col =  tCols[10])
points(hT10$mids,  hT10$counts,  type = "l", col =  tCols[11])
points(hT11$mids,  hT11$counts,  type = "l", col =  tCols[12])
points(hT12$mids, hT12$counts, type = "l", col =  tCols[13])
legend("topright",
       legend = c("hT0","hT1","hT2", "hT3","hT4","hT5", "hT6","hT7","hT8", "hT9","hT10","hT11", "hT12"),
       col = tCols,
       lwd = 1)


#TASK>   Study this plot. What does it tell you? Is there systematic, global
#TASK>   change in the values over time? Within a cycle? Over the course of the
#TASK>   experiment?

#I might have expected hT0, hT6 and hT12 to be more similar in counts compared with
#hT3 and hT9, but this was not the case. This means that really the total gene expression
#Isn't consistantly higher or lower at the start/end versus middle of the cell cycle.
#There does appear to be a separation of values around 1200 and around 1500, but this might
#not appear as significant over many experiments. From the boxplot, I suspect this variation
#is not well correlated will time, or within a cell cycle. When I plotted all time points, there
#still seemed to be a separation of values, however this still didn't seem to have much to
#do with time.



# ==   3.2  Grouped Samples  ===================================================

#  This is the GEO2R code that produces the histogram you saw on the NCBI
#  Website if you went through the BIN-EXPR-GEO unit.

#  Group names for all samples in a series
gsms <- "0123450123450" # Each digit identifies one of the 13 columns
sml <- c()
for (i in 1:nchar(gsms)) {
  sml[i] <- substr(gsms,i,i)
}
sml <- paste("G", sml, sep="")  # set group names

# order samples by group
ex <- exprs(GSE3635)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("t0","t10","t20","t30","t40","t50") # these are the labels we
                                                # assigned in the BIN-EXPR-GEO
                                                # unit

# Set parameters and draw the plot. I changed this from the original GEO2R
# code which overwrote the global palette(). That's evil! Utility code
# should _never_ mess with global parameters!
GEOcols <- c("#dfeaf4", "#f4dfdf", "#f2cb98", "#dcdaa5",
             "#dff4e4", "#f4dff4",  "#AABBCC")
dev.new(width = 4 + dim(GSE3635)[[2]] / 5, height = 6) # plot into a new window
par(mar = c(2 + round(max(nchar(sampleNames(GSE3635))) / 2), 4, 2, 1))
title <- paste ("GSE3635", '/', annotation(GSE3635),
                " grouped samples", sep ='')
boxplot(ex, boxwex = 0.6, notch = TRUE, main = title, outline=FALSE,
        las = 2, col = GEOcols[fl])
legend("topleft", labels, fill = GEOcols, bty = "n")


# =    4  Row-wise Analysis: Expression Profiles  ==============================

# What we did above was column-wise analysis and it gave us an idea about how
# our experiments relate to each other. To analyze individual genes, we need to
# do row-wise analyis.

#TASK>  Try to answer the following questions:
#TASK>     Are all rows genes?
#No, there are also loading controls and other control rows (some are genes, but not yeast genes)
featureNames(GSE3635)[grep("^[^Y]", featureNames(GSE3635))]
#TASK>     What identifiers are being used?
#          (cf. https://sites.google.com/view/yeastgenome-help/community-help/nomenclature-conventions)
#Nuclear ORFs
#TASK>     Are all rows/genes unique?
#Yes
length(unique(featureNames(GSE3635))) == length(featureNames(GSE3635))
#TASK>     Are all yeast genes accounted for?
#No, for example RNA-coding genes and mitochondrial genes are not present.
nNonYeastNuclear <- length(featureNames(GSE3635)[grep("^[^Y]", featureNames(GSE3635))])
nYeastNuclear <- length(featureNames(GSE3635)[grep("^Y", featureNames(GSE3635))])
length(featureNames(GSE3635)) == nNonYeastNuclear + nYeastNuclear

#     These are crucially important questions but you can't answer the last
#     question  because you have no information other than the gene identifiers.
#     However, any biological interpretation relies absolutely on understanding
#     the semantics of the data. Merely manipulating abstract identifiers and
#     numbers, that would surely be Cargo Cult bioinformatics.

#     To answer these questions about the data semantics, I have provided the
#     file "SGD_features.tab" in the data directory of this project. I have
#     downloaded it from SGD
#     (http://www.yeastgenome.org/download-data/curation), for a description of
#     its contents see here:
file.show("./data/SGD_features.README.txt")
#     Note: the file as downloaded from SGD actually crashed RStudio due to an
#           unbalanced quotation mark which caused R to try and read the whole
#           of the subsequent file into a single string. This was caused by an
#           alias gene name (B"). I have removed this abomination
#           by editing the file. The version in the ./data directory can be
#           read without issues.

# Leets peek into the file:
readLines("./data/SGD_features.tab", n = 5)

# ==   4.1  Task - Read a table of features  ===================================

# This data file is rather typical of datasets that you will encounter "in the
# wild". To proceed, you need to write code to read it into an R-object. Develop
# the code in your script file according to the following specification:
#
#    -  read "./data/SGD_features.tab" into a data frame
#          called "SGD_features"
#    -  remove unneeded columns - keep the following information:
#         -   Primary SGDID
#         -   Feature type
#         -   Feature qualifier
#         -   Feature name - (the systematic name !)
#         -   Standard gene name
#         -   Description
#    -  give the data frame meaningful column names:
#    colnames(SGD_features) <- c("SGDID",
#                                "type",
#                                "qual",
#                                "sysName",
#                                "name",
#                                "description")
#
#    -  remove all rows that don't have a systematic name. (You'll have to check
#          what's in cells that don't have a systematic name)
#    -  check that the systematic names are unique (Hint: use the duplicated()
#          function.)
#    -  assign the systematic names as row names
#    -  confirm: are all rows of the expression data set represented in
#                  the feature table? Hint: use setdiff() to print all that
#                  are not.
#                  Example:  A <- c("duck", "crow", "gull", "tern")
#                            B <- c("gull", "rook", "tern", "kite", "myna")
#                            setdiff(A, B)
#                            setdiff(B, A)

#       If some of the features in the expression set are not listed in the
#       systematic names, you have to be aware of that, when you try to get
#       more information on them. I presume they are missing because revisions
#       of the yeast genome after these experiments were done showed that these
#       genes did not actually exist.

#    -  confirm: how many / which genes in the feature table do not
#                have expression data?

#  How should we handle rows/columns that are missing or not unique?

source("../dataParsingGEO2R.R")
head(SGD_features)

# ==   4.2  Selected Expression profiles  ======================================


# Here is an expression profile for Mbp1.

gName <- "MBP1"
(iFeature <- which(SGD_features$name == gName))
(iExprs   <- which(featureNames(GSE3635) == SGD_features$sysName[iFeature]))
plot(seq(0, 120, by = 10),
     exprs(GSE3635)[iExprs, ],
     main = paste("Expression profile for", gName),
     xlab = "time (min)",
     ylab = "expression",
     type = "b",
     col= "maroon")
abline(h =  0, col = "#00000055")
abline(v = 60, col = "#00000055")

# Print the description
SGD_features$description[iFeature]

# Here is a list of gene names that may be involved in the cell cycle switch,
# and some genes that are controls (cf. BIN-SYS-Concepts):

# Turning it on
# Cdc14, Mbp1, Swi6, Swi4, Whi5, Cdc28, Cln1, Cln2, Cln3

# Turning it off
# Rad53, Cdc28, Clb1, Clb2, Clb6, Nrm1

# Housekeeping genes
# Act1, and Alg9

#TASK>  Plot expression profiles for these genes and study them. What do you
#TASK>  expect the profiles to look like, given the role of these genes? What
#TASK>  do you find? (Hint: you will make your life much easier if you define
#TASK>  a function that plots and prints descriptions with a gene name as input.
#TASK>  Also: are the gene names in the feature table upper case, or lower case?
#TASK>  Also: note that the absolute values of change are quite different.
#TASK>  Also: note that some values may be outliers i.e. failed experiments.)

#expect the profiles to spike and dip in a consistent way such that the genes that
#"turn on" are inversely spiking and dipping to the genes that "turn off"

geneList = c("Cdc14", "Mbp1", "Swi6", "Swi4", "Whi5", "Cdc28", "Cln1", "Cln2", "Cln3",
             "Rad53", "Cdc28", "Clb1", "Clb2", "Clb6", "Nrm1",
             "Act1", "Alg9")

plotAndPrint <- function(gene){
  #plot expression
  (iFeature <- which(SGD_features$name == gene))
  (iExprs   <- which(featureNames(GSE3635) == SGD_features$sysName[iFeature]))
  plot(seq(0, 120, by = 10),
      exprs(GSE3635)[iExprs, ],
       main = paste("Expression profile for", gene),
       xlab = "time (min)",
       ylab = "expression",
       type = "b",
       col= "maroon")
  abline(h =  0, col = "#00000055")
  abline(v = 60, col = "#00000055")

  # and print description
  print(SGD_features$description[iFeature])
}

#all gene names in the table are capitalized
for (i in (1:length(geneList))){
  plotAndPrint(toupper(geneList[i]))
}

#prediction was more or less correct, the "on" genes showed a dip around 60,
#while the "off" genes showed a hump. (except for Cln3, which showed a hump
# and Rad53 which showed a dip).

# =    5  Differential Expression  =============================================

# GEO2R discovers the top differentially expressed  expressed genes by
# using functions in the Bioconductor limma package.

if (! require(limma, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("limma")
  library(limma)
}
# Package information:
#  library(help = limma)       # basic information
#  browseVignettes("limma")    # available vignettes
#  data(package = "limma")     # available datasets

# The GEO2R limma code is virtually uncommented, and has not been written for
# clarity, but for being easily produced by a code-generator that is triggered
# by the parameters on the GEO Web-site. Most of it is actually dispensable for
# our purposes.

# In principle, the code goes through three steps:
# 1. Prepare the data
# 2. Define groups that are to be contrasted to define "Differential"
# 3. Find genes whose expression levels are significantly different across
#     the groups
# 4. Format results.

# We are recapitulating the experiment in which we assigned the 0, 10, 60 and
# 70 minute samples to one group, the 30, 40, 90 and 100 minute samples to
# another group, and calculated differential expression values between these
# two groups.

setA <- c(1, 2, 7, 8)    # columns for set A
setB <- c(4, 5, 10, 11)  # columns for set B

# We remove columns we don't need, and for simplicity put the experiments for
# group A in the first four columns, and the group B in the next four.

mySet <- GSE3635[ , c(setA, setB)]

# limma needs the column descriptions as factors
mySet$description <- as.factor(c(rep("A", 4), rep("B", 4)))

# Next we build the "design Matrix" for the statistical test:
myDesign <- model.matrix(~ description + 0, mySet)
colnames(myDesign) <- levels(mySet$description)
myDesign

# Now we can calculate the fit of all rows to a linear model that depends
# on the two groups as specified in the design:
myFit <- lmFit(mySet, myDesign)

# Next we calculate the contrasts, given the fit ...
myCont.matrix <- makeContrasts(A - B, levels = myDesign)
myFit2 <- contrasts.fit(myFit, myCont.matrix)

# ... compute appropriate probabilites from a modified t-test
#     (empirical Bayes) ...
myFit2 <- eBayes(myFit2, 0.01)

# ... add the gene names to the fit - object ...
myFit2$genes <- featureNames(mySet)

# ... and pick the top N differentially expressed genes while controlling
#     for multiple testing with a False Discovery Rate (fdr) correction. GEO2R
#     gave us only the top 250 genes, but we might as well do 1000, just so we
#     can be reasonable sure that our gens of interest are included.
N <- 1000
myTable <- topTable(myFit2, adjust.method = "fdr", sort.by = "B", number = N)

str(myTable)
# The gene names are now in the $ID column

# These are the top 10
write.table(myTable[1:10 , c("ID","P.Value","B")],
            file = stdout(),
            row.names = FALSE,
            sep="\t")

# Let's see what we got: let's plot the full expression profiles (all 13
# columns) for the top ten genes from exprs(GSE3635).

plot(seq(0, 120, by = 10),
     rep(0, 13),
     type = "n",
     ylim = c(-1, 1),
     xlab = "time",
     ylab = "log-ratio expression")
rect( 0,  -2,  15, 2, col = "#dfeaf4", border = NA)   # setA
rect( 55, -2,  75, 2, col = "#dfeaf4", border = NA)   # setA
rect( 25, -2,  45, 2, col = "#f4dfdf", border = NA)   # setB
rect( 85, -2, 105, 2, col = "#f4dfdf", border = NA)   # setB
abline(h = 0, col = "#00000055")

for (i in 1:10) {
  thisID <- myTable$ID[i]
  points(seq(0, 120, by = 10), exprs(GSE3635)[thisID, ], type = "b")
}

# Our guess that we might discover interesting genes be selecting groups A and B
# like we did was not bad. But limma knows nothing about the biology and though
# the expression profiles look good, there is no guarantee that these are the
# most biologically relevant genes. Significantly different in expression
# according to the groups we define is not necessarily the same as a cyclically
# varying gene, nor does it necessarily find the genes whose expression levels
# are _most_ different, i.e. if the variance of a highly, differentially
# expressed gene within a group is large, it may not be very significant. Also,
# we are not exploiting the fact that these values are time series.
# Nevertheless, we find genes for which we see a change in expression levels
# along two cell-cycles.

# Let's superimpose some "real" cell-cycle genes:
myControls <- c("Cdc14", "Mbp1", "Swi6", "Swi4", "Whi5", "Cln1", "Cln2", "Cln3")
for (name in toupper(myControls)) {
  thisID <- SGD_features$sysName[which(SGD_features$name == name)]
  points(seq(0, 120, by=10), exprs(GSE3635)[thisID, ], type="b", col="#AA0000")
}

# Indeed, the discovered gene profiles look much "cleaner" than the real cycle
# genes and this just means that differential expression in the way that we
# have performed it is an approximation to the biology.

# ==   5.1  Final task: Gene descriptions  =====================================

#found "YML117W-A" not in my systematic names, looked up its alias
#YML116W-A and found it there so renamed it
SGD_features$sysName[which("YML116W-A" == SGD_features$sysName)] = "YML117W-A"

#    Print the descriptions of the top ten differentially expressed genes.
for (i in 1:10){
  print(SGD_features$name[SGD_features$sysName == myTable$ID[i]])
  print(SGD_features$description[SGD_features$sysName == myTable$ID[i]])
}



# =    6  Improving on Discovery by Differential Expression  ===================


# There are many ways to improve on purely statistical methods, if we have
# better ideas about the biology. I would just like to demonstrate one
# possibility here: calculate correlation values to a sample gene - such as Cln2
# for which we know that it is involved in the cell cycle. Let's plot it first:

gName <- "CLN2"
(iFeature <- which(SGD_features$name == gName))
(iExprs   <- which(featureNames(GSE3635) == SGD_features$sysName[iFeature]))
Cln2Profile <- exprs(GSE3635)[iExprs, ]
plot(seq(0, 120, by = 10),
     Cln2Profile,
     ylim = c(-1, 1),
     main = paste("Expression profile for", gName),
     xlab = "time (min)",
     ylab = "expression",
     type = "b",
     col= "maroon")
abline(h =  0, col = "#00000055")
abline(v = 60, col = "#00000055")

# Set up a vector of correlation values


myCorrelations <- numeric(nrow(exprs(GSE3635)))
names(myCorrelations) <- featureNames(GSE3635)
for (i in 1:length(myCorrelations)) {
  myCorrelations[i] <- cor(Cln2Profile, exprs(GSE3635)[i, ])
}

myTopC <- order(myCorrelations, decreasing = TRUE)[1:10]  # top ten

# Number 1
(ID <- featureNames(GSE3635)[myTopC[1]])

# Get information
SGD_features[which(SGD_features$sysName == ID), ]
# Of course: the highest correlation is Cln2 itself. This is our positive
# control for the experiment.

# Let's plot the rest
for (i in 2:length(myTopC)) {
  ID <- featureNames(GSE3635)[myTopC[i]]
  points(seq(0, 120, by = 10),
       exprs(GSE3635)[ID, ],
       type = "b",
       col= "chartreuse")
  print(SGD_features[which(SGD_features$sysName == ID), c("name", "description")])
}

# Note that all of these genes are highly correlated with a known cell cycle
# gene, but because the absolute values of expression differences are not very
# large for some of them, they might not be picked up by any algorithm that is
# focussed on large differential expression changes. Do small relative changes
# mean small biological effects? Certainly not!

# And we haven't even looked at the anticorrelated genes yet...
myBottomC <- order(myCorrelations, decreasing = FALSE)[1:10]  # bottom ten
for (i in 1:length(myBottomC)) {
  ID <- featureNames(GSE3635)[myBottomC[i]]
  points(seq(0, 120, by = 10),
         exprs(GSE3635)[ID, ],
         type = "b",
         col= "coral")
  print(SGD_features[which(SGD_features$sysName == ID), c("name", "description")])
}
# ... which are very interesting in their own right.

# What I hope you appreciate from this example is:
#   -  convenient general purpose methods exist that analyze expression data
#        with sophisticated statistical methods;
#   -  the results of these methods depend on whether the statistics model
#        the biology well;
#   -  you will draw Cargo Cult conclusions if you can't cast your biology
#        of interest in terms of a model, "canned" solutions are unlikely
#        to give you all the answers you need;
#   -  being able to write your own code gives you the freedom to experiment
#        and explore. There is a learning curve - but the payoffs are
#        significant.

# =    7  Annotation data  =====================================================
#
# Loading feature data "by hand" as we've done above, is usually not necessary
# since GEO provides rich annotations in the GPL platform files, which are
# associated with its Gene Expression Sets files. In the code above,
# we used getGEO("GSE3635", GSEMatrix = TRUE, getGPL = FALSE), and the GPL
# annotations were not loaded. We could use getGPL = TRUE instead ...

GSE3635annot <- getGEO("GSE3635", GSEMatrix = TRUE, getGPL = TRUE)
GSE3635annot <- GSE3635annot[[1]]

# ... and the feature data is then available in the GSE3635@featureData@data
#     slot:
str(GSE3635annot@featureData@data)
GSE3635annot@featureData@data[ 1:20 , ]

# ... from where you can access the columns you need, e.g. spotIDs and gene
#     symbols:

myAnnot <- GSE3635annot@featureData@data[ , c("SPOT_ID", "Gene")]
str(myAnnot)

# ... Note that this is a data frame, but - oy veh - the gene symbols are
#     factors. Really, we need to fix this! To convert a factor into a string,
#     we need to cast it to character.

myAnnot$Gene <- as.character(myAnnot$Gene)

# ... whereupon we can find things we might be looking for ...
myAnnot[which(myAnnot$Gene == "MBP1"), ]

# ... or identify rows that might give us trouble, such as probes that
#     hybridize to more than one gene.


# Alternatively, we could have identified the GPL file for this set:
GSE3635@annotation   # "GPL1914"

# ... and downloaded it directly from NCBI:
GPL1914 <- getGEO("GPL1914")
str(GPL1914)

# ... from which we can get the data - which is however NOT necessarily

# matched to the rows of our expression dataset. Note that here to: the majority
# of data elements are factors and will likely have to be converted before
# use.



# [END]
