install.packages("BiocManager")
library(BiocManager)
install(c("genbankr","Biostrings","ggtree","annotate")) 
install.packages(c("reshape2","rentrez"))
library(genbankr)
library(annotate)
library(ape)
library(ggtree)


BbID<-GBAccession("Q9ZZ64")#using NCBI I seached for CO1, and canis lupus familiaris [organism], i find this accession number
BbGbk<-readGenBank(BbID)##use the accession number to browse through the gene bank 

##Pairwise alignment
BbGbkBLAST<-blastSequences(paste(BbGbk@sequence),as = 'data.frame',hitListSize = 40, timeout = 600)## pairwise alignment search using our gene 
##multiple alignment
BbHitsDF<-data.frame(ID=BbGbkBLAST$Hit_accession,Seq=BbGbkBLAST$Hsp_hseq,stringsAsFactors = FALSE)


##muscle+ alignemnt inspection
BbHitSeqs<-read.GenBank(BbGbkBLAST$Hit_accession[1:3])
BbHitsDNA<-sapply(BbHitsDF$Seq,strsplit,split="")## convert the object to one with seperate column for each basepair using strsplit()and sapply() to split the DA for each row
names(BbHitsDNA)<-paste(1:nrow(BbHitsDF),BbHitsDF$ID,sep="_")## give each sequence unique names.
BbHitsDNA<-as.DNAbin(BbHitsDNA)## convert it into a DNAbin object 
BbAlign<-muscle(BbHitsDNA,quiet=F)## run muscle on DNAbin object, allowing alignment 
checkAlignment(BbAlign,what=3)## to inspect alignments and determine Shannon indexes relatives to the sequence position, providing cues to diversity

#Distance Matreix 
BbDM<-dist.dna(BbAlign, model="K80") #use distance model K80, this estimates a pairwise distance matrix from the sequence data
BbDMmat<-as.matrix(BbDM)
## tree building 
BbTree<-nj(BbDM)
ggtree(BbTree,layout="rectangular") + geom_tiplab()


