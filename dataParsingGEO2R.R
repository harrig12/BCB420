#Wild data parsing
#As specified in RPR-GEO2R.R

#    -  read "./data/SGD_features.tab" into a data frame
#          called "SGD_features"

# Data is one entry per line, and tab delimited.
SGD_features <- read.delim("./data/SGD_features.tab", header = F, stringsAsFactors = F)

#    -  remove unneeded columns - keep the following information:
#         -   Primary SGDID
#         -   Feature type
#         -   Feature qualifier
#         -   Feature name - (the systematic name !)
#         -   Standard gene name
#         -   Description
head(SGD_features)
SGD_features <- SGD_features[c(1,2,3,4,5,16)]

#    -  give the data frame meaningful column names:
colnames(SGD_features) <- c("SGDID",
                            "type",
                            "qual",
                            "sysName",
                            "name",
                            "description")
#
#    -  remove all rows that don't have a systematic name. (You'll have to check
#          what's in cells that don't have a systematic name)

SGD_features <- SGD_features[SGD_features$sysName!="",]

#    -  check that the systematic names are unique (Hint: use the duplicated()
#          function.)

sum(duplicated(SGD_features$sysName)) == 0

#    -  assign the systematic names as row names

rownames(SGD_features) <- SGD_features$sysName

#    -  confirm: are all rows of the expression data set represented in
#                  the feature table? Hint: use setdiff() to print all that
#                  are not.
#                  Example:  A <- c("duck", "crow", "gull", "tern")
#                            B <- c("gull", "rook", "tern", "kite", "myna")
#                            setdiff(A, B)
#                            setdiff(B, A)
exprDataFeatureNames <- featureNames(GSE3635)

setdiff(exprDataFeatureNames, rownames(SGD_features))

#       If some of the features in the expression set are not listed in the
#       systematic names, you have to be aware of that, when you try to get
#       more information on them. I presume they are missing because revisions
#       of the yeast genome after these experiments were done showed that these
#       genes did not actually exist.

#    -  confirm: how many / which genes in the feature table do not
#                have expression data?

length(setdiff(rownames(SGD_features), exprDataFeatureNames))

#  How should we handle rows/columns that are missing or not unique?
# See if we can account for them in synonyms or aliases and remove the data that doesn't make sense
# or is the least well supported.

# [END]
