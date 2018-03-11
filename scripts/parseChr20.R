# parseChr20.R
#
# Purpose: Script to create and save a chromR object from PGP downlads for chromosome 20
#
# License: GPL-3 (https://www.gnu.org/licenses/gpl-3.0.en.html)
#
# Author: Cait Harrigan <cait.harrigan@mail.utoronto.ca>
#
# Notes: Associated with personal genomes unit (2.3) in my BCB420 journal
# http://steipe.biochemistry.utoronto.ca/abc/students/index.php/User:Caitlin_Harrigan/Journal420
#
# ==============================================================================
#!/usr/bin/Rscript
args <- commandArgs(TRUE)

# test for command line arguments
# (modified from https://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/)
if (length(args)!=1) {
  stop("Please supply only a .vcf or vcf.gz file. Script stopped.", call.=FALSE)
} else if (length(args)==1) {

  if (!require(vcfR, quietly=TRUE)) {
    install.packages("vcfR")
    library(vcfR)
  }

  if (!require(seqinr, quietly=TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("seqinr")
    library(seqinr)
  }

  # read in vcf information
  parData <- read.vcfR(args[1])

  # subset for just chromosome 20
  sel <- getCHROM(parData) == "chr20"

  #parse file name for participant id
  parID <- c2s(s2c(args[1])[16:18])
  parName <- paste0("chr20_", parID)

  #create chromosome object
  CH20 <- create.chromR(parData[sel], parName)

  #save chromsome object
  dir.create(file.path(getwd(), "data/PGP/"), showWarnings = FALSE)
  setwd(file.path(getwd(), "data/PGP/"))

  save(CH20, file = paste0(parName, ".Rdata"))
}

#[END]
