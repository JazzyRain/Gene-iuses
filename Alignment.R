install.packages("BiocManager")
library(BiocManager)
install(c("genbankr","Biostrings","ggtree","annotate")) 
install.packages(c("reshape2","rentrez"))
library(genbankr)
library(annotate)
library(ape)

BbID<-GBAccession("KC984225.1")#using NCBI I seached for 16S ribosomal RNA, and canis lupus familiaris [organism], i find this accession number
BbGbk<-readGenBank(BbID)##use the accession number to browse through the gene bank 

##Pairwise alignment
BbGbkBLAST<-blastSequences(paste(BbGbk@sequence),as = 'data.frame',hitListSize = 40, timeout = 600)## pairwise alignment search using our gene 
##multiple alignment
BbHitsDF<-data.frame(ID=BbGbkBLAST$Hit_accession,Seq=BbGbkBLAST$Hsp_hseq,stringsAsFactors = FALSE)


library(rtracklayer)
test_path <- system.file("tests", package = "rtracklayer")
test_bed <- file.path(test_path, "test.bed")

test <- import(test_bed, format = "bed")
test