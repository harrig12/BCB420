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

if (!require(stringr, quietly=TRUE)) {
  install.packages("stringr")
  library(stringr)
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

#========= slow BWT method ==================================================

#We can reconstruct the original string from its BWT by repeatedly sorting a
#list of strings then prepending the BWT.
# (from http://stanford.edu/class/cs262/notes/lecture4.pdf )

#build original string from BWT
slowBWT2s <- function(sBWT){
  #sBWT is a BW-transformed string
  
  stringO <- character()
  dfBWT <- as.data.frame(sBWT)
  
  #prepend the BWT to sorted "in progress" string until original 
  #string is reformed
  for (i in 1:nrow(dfBWT)){
    stringO <- paste0(dfBWT[,1], sort(stringO))
    dfBWT <- cbind(sBWT, stringO)
  }

  #the original string is that one row above the one begining
  #with the last character in the string - $ should be unique in the string
  return(dfBWT[which(dfBWT[,1] == "$") - 1, 2])
}

#did we get the initial string back?
stringI <- slowBWT2s(s2BWT(mySeq, mySA))

(c2s(mySeq) == stringI) #True


#========= fast BWT method ==================================================

lexiScore <- function(s, char){
  #s is a character vector
  #char is the character to score
  
  #find the number of characters lexicographically smaller than each 
  #character in s
  
  #subtract 1 for indexing from 1 in R
  return (min(which(sort(s)==char)) - 1)
}

getIthOccurence <- function(s){
  #s is a character vector
  
  #find the ith occurence of each character in s
  
  #keep track of how many of each character we've seen
  occCounter <- numeric(length(unique(s)))
  names(occCounter) <- unique(s)
  
  #place the occurence into the ith position corresponding to
  #the character in s
  ith <- numeric()
  
  i = 1
  for (c in s){
    occCounter[c] = occCounter[c] + 1
    ith[i] = occCounter[c]
    i = i + 1
  }
  
  return(ith)
  
}

computeLF <- function(cvBWT){
  #cvBWT is a BW-transformed character vector of the initial string
  
  #find lexicographical score of each character in sBWT
  #this is "C(a)" in the referenced material
  lexiScores <- numeric()
  for (c in cvBWT){
    lexiScores <- c(lexiScores, lexiScore(cvBWT, c))
  }
  ithOccurence <- getIthOccurence(cvBWT)
  return(lexiScores + ithOccurence)
}

# taken directly from RECONSTRUCT() in http://stanford.edu/class/cs262/notes/lecture4.pdf
fastBWT2s <- function(cvBWT){
  #cvBWT is a BW-transformed character vector of the initial string
  
  LF <- computeLF(cvBWT)
  S <- character()
  r <- 1
  ch <- cvBWT[r]
  while (ch != '$'){
    S <- c(ch, S)
    r <- LF[r]
    ch <- cvBWT[r]
  }
  return (S)
}


#search a substring
#to fins a substring, we look for the lowest and highest index of an occurence
#of a prefix W

#build the dataframe from which we found the string from the BWT
#very similar to slowBWT2s
rotateBWT <- function(sBWT){
  
  stringO <- character()
  dfBWT <- as.data.frame(sBWT)
  
  for (i in 1:nrow(dfBWT)){
    stringO <- paste0(dfBWT[,1], sort(stringO))
    dfBWT <- cbind(sBWT, stringO)
  }
  
  return(dfBWT[,2])
}

getLowAndHighW <- function(BWTrotations, W){
  #reduce each string to a prefix
  for (indPrefix in 1:length(BWTrotations)){
    BWTrotations[indPrefix] <- c2s(s2c(BWTrotations[indPrefix])[1:nchar(W)])
  }
  #find the indices of occurences of W in the prefix array
  indW <- which(BWTrotations == W)
  
  #return the first and last occurences
  return (c(indW[1], indW[length(indW)]))
}



# [END]