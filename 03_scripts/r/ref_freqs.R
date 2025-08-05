
REF = commandArgs(TRUE)[1]
suppressPackageStartupMessages(library(Biostrings))

g = readDNAStringSet(REF); f = colSums(trinucleotideFrequency(g))

f2 = sapply(which(substr(names(f), 2, 2) %in% c("C", "T")), function(i) {
    rc = as.character(reverseComplement(DNAString(names(f)[i]))); f[i] + f[rc];
});

write.table(as.matrix(f2), file="ref_freqs.txt", sep="\t", col.names=F)
