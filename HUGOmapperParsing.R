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
#certain genes has no HUGO equivalents, these tended to be antibody receptors



