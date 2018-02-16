#script to handle HUGO multisymbol mapper output

#write these to a file to use other mapping techniques to match them to HUGO symbols
#those with approved symbols not in HUGOsymbols should be removed
#those which are withdrawn should be removed
#those wich are unmatched should be removed and manually searched
#those with synonyms or previous symbols should be placed into the synMap

#read in mappings and update synMap
HUGOmap <- read.delim("symHUGOmultiMapper.txt", header = T, sep="\t",stringsAsFactors = F)
unmatched <- HUGOmap[(HUGOmap$Match.type == "Unmatched"), c("Input")]
mappable <- HUGOmap[(HUGOmap$Match.type == "Previous symbol"| HUGOmap$Match.type =="Synonyms"),
                     c("Input", "Approved.symbol")]
colnames(mappable) <- c("symbols", "synonyms")
synMap <- rbind(synMap, mappable)

#manual filtering of unmatched:
#certain genes has no HUGO equivalents
manualMatched <- read.delim("manMatchGSE84712.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(manualMatched) <- c("symbols", "synonyms")
synMap <- rbind(synMap, manualMatched)

#some of the symbols added are not in HUGOsymbols, so remove these from synMap
synMap <- synMap[(synMap$symbols %in% HUGOsymbols),]

#[END]
