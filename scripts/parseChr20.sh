#!/bin/bash
for fn in data/PGP/*.vcf.gz
	do
	echo "parsing ch20 for ${fn}"
	Rscript scripts/parseChr20.R ${fn}
	done
