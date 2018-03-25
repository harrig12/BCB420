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
  install.packages("Biostrings")
  library(Biostrings)
}

#load some gene of interest
mySeq <- readDNAStringSet('data/EGFrefseq.fa')

#Select some query... the 10th line of the fasta file will do fine. (the 701-770th bases)
myQuery <- 'CGATGGTGTGGGAGTGAAGGCTCTGTTGGAGACATCAGAGAAAATAACAGCTGTGTCATTGGATGTGCTT' 


# [END]