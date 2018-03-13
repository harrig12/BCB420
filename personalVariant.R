# personalVariant.R
#
# Functions for downloading and exploring personal genome canada project data.
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Cait Harrigan <cait.harrigan@mail.utoronto.ca>
#
# Purpose: Explore variant density in and out of protein-coding genes
#
# Versions: 1.1 Stream with chr20 parsing scripts
#           1.0 Initial exploration
#
# Notes: Associated with personal genomes unit (2.1) in my BCB420 journal
# http://steipe.biochemistry.utoronto.ca/abc/students/index.php/User:Caitlin_Harrigan/Journal420
#
# ==============================================================================

if (!require(vcfR, quietly=TRUE)) {
  install.packages("vcfR")
  library(vcfR)
}

if (!require(readr, quietly=TRUE)) {
  install.packages("readr")
  library(readr)
}

if (!require(intervals, quietly=TRUE)) {
  install.packages("intervals")
  library(intervals)
}

CHRLEN20 <- 64444167
CH20CENTROMERE <- c(26369569, 29369569)

# load in vcf information
load("../data/PGP/chr20_001.Rdata")
CH20_001 <- CH20
rm(CH20)

#lets get to know our chromosome
show(CH20_001)
head(CH20_001)
head(variant.table(CH20_001))

#looks like PGP didn't document MQ (mapping quality) or DP (depth)
#we can still look at the quality measure
chromoqc(CH20_001)

#read in ensemble data (download details documented in journal)
tmp <- read_tsv("../data/HGNC_chr20.txt")

#check for uniqness in the HGNC symbol
sum(duplicated(tmp$`HGNC symbol`))  #0

#subset as appropriate
Chr20GeneData <- data.frame(sym = tmp$`HGNC symbol`,
                            start = tmp$`Gene start (bp)`,
                            end = tmp$`Gene end (bp)`,
                            stringsAsFactors = FALSE)

#check all elements have coordinates on the chromosome (0 <= coord <= CHRLEN20)

#genes
(min(Chr20GeneData$end) >= 0)  #True
(max(Chr20GeneData$end) <= CHRLEN20)  #True

#variants
(min(CH20_001@var.info$POS) >= 0)  #True
(max(CH20_001@var.info$POS) <= CHRLEN20)  #True

#quick plot
geneIntervals <- Intervals(Chr20GeneData[, c("start", "end")])
rownames(geneIntervals) <- Chr20GeneData$sym

plot(geneIntervals,
     use_points = F,
     xlab = "Coordinate",
     ylab = "n Overlaps",
     xlim = c(0, CHRLEN20),
     ylim = c(0,2),
     lwd = 5,
     use_names = F,
     main = "Chromosome 20 participant 001 genes",
     axes = T)
axis(2, at = c(0,1,2))

#draw centromere
rect(CH20CENTROMERE[1], -0.05, CH20CENTROMERE[2], 2.05, col = "#ff00007F", border = NA)
legend(x = 4e+07, y = 0.5, c("HGNC gene", "Centromere"), col=c(1,"#ff00007F"), lwd=3)

#find start and end coordinates of non-genes (interval complement of gene intervals)
nonGeneIntervals <- interval_complement(geneIntervals)

#adjuset the starts and ends for chromosome coordinates
nonGeneIntervals[[1]] <- 0
nonGeneIntervals[length(nonGeneIntervals)/2][[2]] <- CHRLEN20

#make colummn to flag genes vs non-genes
Chr20GeneData$gFlag <- rep(1, nrow(Chr20GeneData))

#add non-genes to chr20 data
nonGeneNames <- paste0("notAGene", 1:nrow(nonGeneIntervals))
nonGenes <- as.data.frame(cbind(nonGeneNames, as.matrix(nonGeneIntervals), 0))
names(nonGenes) <- names(Chr20GeneData)
Chr20GeneData <- rbind(Chr20GeneData, nonGenes)

#init nVariants column
Chr20GeneData$nVariants <- rep(0, nrow(Chr20GeneData))

#count the number of variants per gene
#takes a few seconds
for (var_i in 1:length(CH20_001@var.info$POS)){
  pos <- CH20_001@var.info$POS[var_i]
  inGenes <- names(which(geneIntervals[,1]<=pos & geneIntervals[,2]>=pos))
  for (g in inGenes){ #should not be more than 3 for any gene
    Chr20GeneData$nVariants[Chr20GeneData$sym == g] <- Chr20GeneData$nVariants[Chr20GeneData$sym == g] + 1
  }
}

#count the number of variants per non-gene
#takes a few seconds
for (var_i in 1:length(CH20_001@var.info$POS)){
  pos <- CH20_001@var.info$POS[var_i]
  inNonGenes <- paste0("notAGene", (which(nonGeneIntervals[,1]<=pos & nonGeneIntervals[,2]>=pos)))
  for (g in inNonGenes){ #should not be more than 1 for any non-gene
    Chr20GeneData$nVariants[Chr20GeneData$sym == g] <- Chr20GeneData$nVariants[Chr20GeneData$sym == g] + 1
  }
}

#plot raw n variants per gene
plot(1:nrow(Chr20GeneData), Chr20GeneData$nVariants, pch=18)

#find where to plot the centromere
cenLoc <- length(which(Chr20GeneData$end<CH20CENTROMERE[1]))
abline(v = cenLoc, col =2)

# [END]
