install.packages("BiocManager")
library(BiocManager)
install(c("genbankr","Biostrings","ggtree","annotate")) 
install.packages(c("reshape2","rentrez"))
library(genbankr)
library(annotate)

BbID<-GBAccession("KC984225.1")#using NCBI I seached for 16S ribosomal RNA, and canis lupus familiaris [organism], i find this accession number
BbGbk<-readGenBank(BbID)##use the accession number to browse through the gene bank 

##Pairwise alignment
BbGbkBLAST<-blastSequences(paste(BbGbk@sequence),as = 'data.frame',hitListSize = 40, timeout = 600)## pairwise alignment search using our gene 

### maybe we can exclude canis lupus familiaris to generate a better tree 
