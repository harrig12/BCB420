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
# Versions: 1.2 Plot 8 participants using makeChrElements.R
#           1.1 Stream with chr20 parsing scripts
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
load("data/PGP/chr20_001.Rdata")
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
tmp <- read_tsv("data/HGNC_chr20.txt")

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

#adjust the starts and ends for chromosome coordinates
nonGeneIntervals[[1]] <- 0
nonGeneIntervals[length(nonGeneIntervals)/2][[2]] <- CHRLEN20

#make colummn to flag genes vs non-genes
Chr20GeneData$gFlag <- rep(1, nrow(Chr20GeneData))

#add non-genes to non-gene intervals and chr20 data  
nonGeneNames <- paste0("notAGene", 1:nrow(nonGeneIntervals))
rownames(nonGeneIntervals) <- nonGeneNames
nonGenes <- cbind(nonGeneNames, as.data.frame(nonGeneIntervals), 0)
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

#sort elements of ch20 by coordinate
Chr20GeneData <- Chr20GeneData[order(Chr20GeneData$start, Chr20GeneData$end),]

#plot raw n variants per gene
plot(1:nrow(Chr20GeneData),
     Chr20GeneData$nVariants,
     col = ifelse(Chr20GeneData$gFlag==1, 1, 3),
     pch=18,
     xlab = "Gene Position",
     ylab = "n Variants",
     main = "Raw Variant Counts Per Gene on Chomosome 20 Participant 1")

#find where to plot the centromere
cenLoc <- length(which(Chr20GeneData$end<CH20CENTROMERE[1]))
abline(v = cenLoc, col =2)
legend("topright", c("Gene", "Non-Gene", "Centromere"),
       pch=c(18, 18, NA),
       lty=c(NA, NA, 1),
       col=c(1,3,2))

#normalize variant counts to gene length
#use two different methods for indexing columns Chr20GeneData because of a quirk of sweep()
#index with $end to rename
Chr20GeneData$gLen <- sweep(Chr20GeneData["end"], 1,
                            Chr20GeneData$start, FUN="-")$end
Chr20GeneData$normCounts <- sweep(Chr20GeneData["nVariants"], 1,
                                  Chr20GeneData$gLen, FUN="/")$nVariants

#plot the counts normalized for gene size
plot(1:nrow(Chr20GeneData),
     Chr20GeneData$normCounts,
     col = ifelse(Chr20GeneData$gFlag==1, 1, 3),
     pch=18,
     xlab = "Position",
     ylab = "Variants per Basepair",
     main = "Length Normalized Variant Counts \n Per Element on Chomosome 20 Participant 1")

#find where to plot the centromere
cenLoc <- length(which(Chr20GeneData$end<CH20CENTROMERE[1]))
abline(v = cenLoc, col =2)
legend("topright", c("Gene", "Non-Gene", "Centromere"),
       pch=c(18, 18, NA),
       lty=c(NA, NA, 1),
       col=c(1,3,2))

#the most highly varying non-genes seem to be near the centromere
#there does not appear to be a particular clustering of highly varying genes.

#plot distribution of variants per bp in genes
hist(Chr20GeneData$normCounts[Chr20GeneData$gFlag==1],
     col = '#00000090',
     main = 'Variants in Gene and Non-Gene elements \nChromosome 20 Participant 001',
     xlab = 'Variants per bp')

#plot distribution of variants per bp in non-genes
hist(Chr20GeneData$normCounts[Chr20GeneData$gFlag==0], add=T,
     col = '#00FF0090')

legend("topright", c("Gene", "Non-Gene"),
       lty=c(1, 1),
       lwd=c(3,3),
       col=c('#00000090', '#00FF0090'))

#these distributions looks quite similar, let's check what the qq plot looks like
qqplot(Chr20GeneData$normCounts[Chr20GeneData$gFlag==0], 
       Chr20GeneData$normCounts[Chr20GeneData$gFlag==1],
       main = 'Gene and Non-gene Variants per bp qq-plot',
       xlab = 'Non-gene quantiles',
       ylab = 'Gene quantiles')
       
abline(0, 1)

#the distributions still seem fairly similar, particularly near the the lower end. 
#there is some divergence and for chromosome segments with more variants, they seem 
#to more frequently represent a gene. This is also supported by the normalized
#counts plots.

boxplot(normCounts~gFlag, 
        Chr20GeneData,
        main = 'Gene Length Normalized Variant Counts \nChromosome 20 Participant 001',
        ylab = 'Variants per bp',
        col = c('#00FF0090', '#00000090'),
        axes = F)
axis(2)
axis(1, at = c(1:2), tick = F, labels = c('Non-Gene', 'Gene'))



#plot the first 8 participants!

source("makeChrElements.R")

#setup 
for (nP in 1:8){
  print(paste0('working on participant ', nP))
  
  #load the chr dataframe
  load(paste0('data/PGP/chr20_00', nP, '.Rdata'))
  
  #set up the elements
  assign(paste0('CHR20_00', nP), makeChrEle(CH20))
  rm(CH20)
}

#plot
#prepare labels for box plot

boxLabels <- NULL

for (nP in 1:8){
  boxLabels <- c(boxLabels, 
                 paste0('Participant ', nP)) 
}

#draw box plot
for (nP in 1:8){
  
  participantN <- paste0('CHR20_00', nP)
  
  #draw the first boxplot
  if (nP == 1){
    boxplot(normCounts~gFlag, 
            get(participantN), 
            col = c('#00FF0090', '#00000090'),
            xlim = c(0,17),
            ylim = c(0, 0.03),
            axes = F,
            main = 'Gene Length Normalized Variant Counts \nChromosome 20, 8 Participants from PGP',
            ylab = 'Variants per bp'
    )
    axis(1, tick = F, at = (1:8)*2 - 0.5, labels = boxLabels, las = 2, hadj = 0.8)
    axis(2)
    legend("topright", c("Gene", "Non-Gene"),
           lty=c(1, 1),
           lwd=c(3,3),
           col=c('#00000090', '#00FF0090'))
  }
  
  #add the others
  else{
    boxplot(normCounts~gFlag, 
            get(participantN), 
            col = c('#00FF0090', '#00000090'),
            add = T,
            at = 1:2 + 2*(nP-1),
            axes = F
    )  
  }
  
} 

#Distributions appear similar for the 8 participants

#Are the elements with high variant counts clustered?

thresh <- 0.005

highVarIntervals <- Intervals(cbind(Chr20GeneData$start[Chr20GeneData$normCounts > thresh],
                                    Chr20GeneData$end[Chr20GeneData$normCounts > thresh]))

rownames(highVarIntervals) <- Chr20GeneData$sym[Chr20GeneData$normCounts > thresh]


#colour by gene or non-gene 
#caution, the intervals are ordered differently than in the dataframe!
highVarCol <-  (Chr20GeneData[which(Chr20GeneData$sym %in% rownames(highVarIntervals)), 
                             'gFlag'] * -2) + 3

plot(highVarIntervals,  
     y = c(0,0),
     use_points = F,
     xlab = "Coordinate",
     xlim = c(0, CHRLEN20),
     lwd = 5,
     use_names = F,
     main = "Chromosome 20 participant 001 \nElements with > 0.005 Variants per bp",
     col = highVarCol, 
     axes = T)

#draw centromere
rect(CH20CENTROMERE[1], -0.05, CH20CENTROMERE[2], 0.05, col = "#ff00007F", border = NA)

legend("topright", c("Gene", "Non-Gene", "Centromere"),
       pch=c(18, 18, 15),
       col=c(1,3,'#FF0000BB'))

#Are these highly varying genes particularly associated with disease states?
#Data stolen from BioHacks repo. See github.com/hyginn/BCBBH-2018 for data prep information

#Check GWAStraits 
GWAStraits <- read_tsv('data/CHr20GWAStraits.tsv')
GWAStraits[which(GWAStraits$sym %in% rownames(highVarIntervals)),]

#nothing interesting, it seems

#Check HPAprognostic
HPAprog <- read_tsv('data/Chr20GeneData.tsv', )
HPAprog <- HPAprog[,c("sym", "HPAprognostic")]

HPAprog[which(HPAprog$sym %in% rownames(highVarIntervals)),]

#two cancers make an appearance; liver and renal.

#expand the search - reduce the threshold!



# [END]

