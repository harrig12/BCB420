# BWT.R
#
# Functions for finding a query in a DNA srting, making use of BWT.
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Cait Harrigan <cait.harrigan@mail.utoronto.ca>
#
# Purpose: Explore BWT 
#
# Versions: 1.0 BWT exploration
#
# Notes: Associated with population genomes unit (2.2) in my BCB420 journal
# http://steipe.biochemistry.utoronto.ca/abc/students/index.php/User:Caitlin_Harrigan/Journal420
#
# ==============================================================================

if (!require(Biostrings, quietly=TRUE)) {
  biocLite("Biostrings")
  library(Biostrings)
}

if (!require(seqinr, quietly=TRUE)) {
  biocLite("seqinr")
  library(seqinr)
}

seqFile = 'data/seqBanana.fa'

#load some sequence of interest (if fasta file)
mySeq <- read.fasta(seqFile)
mySeq <- unlist(mySeq)

#load some sequence of interest (if txt file)
# mySeq <- readLines(seqFile)

#compute the suffix array (method from https://support.bioconductor.org/p/24049/ 
#and in consultation of http://stanford.edu/class/cs262/notes/lecture4.pdf )
mySA <- order(rank(mySeq, ties.method = 'last'))

#We can construct BWT(X) from the suffix array S. At each position, take
#the letter to the left of the one pointed to by S.
# (from http://stanford.edu/class/cs262/notes/lecture4.pdf )

#find BWT
s2BWT <- function(s, suffixArray){ 
  # s is a character array
  # suffixArray gives the indices of the suffix array of s in order
  stringBWT <- character()
  for (i in suffixArray){
    i <- i-1
    if (i == 0){i <- max(suffixArray)}
    stringBWT <- c(stringBWT, s[i])
  }
  return(stringBWT)
}

#We can reconstruct the original string from its BWT by repeatedly sorting a
#list of strings then prepending the BWT.
# (from http://stanford.edu/class/cs262/notes/lecture4.pdf)

#build original string from BWT
BWT2s <- function(sBWT){
  
  stringO <- character()
  dfBWT <- as.data.frame(sBWT)
  
  for (i in 1:nrow(dfBWT)){
    stringO <- paste0(dfBWT[,1], sort(stringO))
    dfBWT <- cbind(sBWT, stringO)
  }

  return(dfBWT[which(dfBWT[,1] == "$") - 1, 2])
}

#did we get the initial string back?
stringI <- BWT2s(s2BWT(mySeq, mySA))

(c2s(mySeq) == stringI) #True

#search a substring
#to fins a substring, we look for the lowest and highest index of an occurance
#of a prefix W

#build the dataframe from which we found the string from the BWT
#very similar to BWT2s
permuteBWT <- function(sBWT){
  
  stringO <- character()
  dfBWT <- as.data.frame(sBWT)
  
  for (i in 1:nrow(dfBWT)){
    stringO <- paste0(dfBWT[,1], sort(stringO))
    dfBWT <- cbind(sBWT, stringO)
  }
  
  return(dfBWT[,2])
}

getLowAndHighW <- function(BWTpermutations, W){
  #reduce each string to a prefix
  for (indPrefix in 1:length(BWTpermutations)){
    BWTpermutations[indPrefix] <- c2s(s2c(BWTpermutations[indPrefix])[1:nchar(W)])
  }
  #find the indices of occurances of W in the prefix array
  indW <- which(BWTpermutations == W)
  
  #return the first and last occurances
  return (c(indW[1], indW[length(indW)]))
}

# [END]