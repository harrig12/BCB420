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
# Versions:
#           1.0 Initial exploration
#
# Notes: Associated with personal genomes unit (2.3) in my BCB420 journal
# http://steipe.biochemistry.utoronto.ca/abc/students/index.php/User:Caitlin_Harrigan/Journal420
#
# ==============================================================================

if (!require(vcfR, quietly=TRUE)) {
  install.packages("vcfR")
  library(vcfR)
}

# load in vcf information
par001 <- read.vcfR("../data/PGP/PGPC_0001_S1.flt.vcf.gz")

# subset for just chromosome 20
sel <- getCHROM(par001) == "chr20"

#create chromosome object
CH20_001 <- create.chromR(par001[sel], "Ch20_001")

#lets get to know our chromosome
show(CH20_001)
head(CH20_001)
head(variant.table(CH20_001))

#looks like PGP didn't document MQ (mapping quality) or DP (depth)
#we can still look at the quality measure
chromoqc(CH20_001)




# [END]
