install.packages("BiocManager")
library(BiocManager)
install(c("genbankr","Biostrings","ggtree","annotate")) 
install.packages(c("reshape2","rentrez"))
library(genbankr)
library(annotate)
library(ape)
library(ggtree)

BbID<-GBAccession("NC_002008.4")#using NCBI I seached for CO1, and canis lupus familiaris [organism], i find this accession number
BbGbk<-readGenBank(BbID)##use the accession number to browse through the gene bank 

##Pairwise alignment
BbGbkBLAST<-blastSequences(paste(BbGbk@sequence),as = 'data.frame',hitListSize = 40, timeout = 600)## pairwise alignment search using our gene 
##multiple alignment
BbHitsDF<-data.frame(ID=BbGbkBLAST$Hit_accession,Seq=BbGbkBLAST$Hsp_hseq,stringsAsFactors = FALSE)

##muscle+ alignemnt inspection+ generate distance matrix 
BbHitSeqs<-read.GenBank(BbGbkBLAST$Hit_accession[1:20])
BbHitsDNA<-sapply(BbHitsDF$Seq,strsplit,split="")## convert the object to one with seperate column for each basepair using strsplit()and sapply() to split the DA for each row
names(BbHitsDNA)<-paste(1:nrow(BbHitsDF),BbHitsDF$ID,sep="_")## give each sequence unique names.
BbHitsDNA<-as.DNAbin(BbHitsDNA)## convert it into a DNAbin object 
BbAlign<-muscle(BbHitsDNA,quiet=F)## run muscle on DNAbin object, allowing alignment 
checkAlignment(BbAlign,what=3)## to inspect alignments and determine Shannon indexes relatives to the sequence position, providing cues to diversity
BbDM<-dist.dna(BbAlign, model="K81") #use distance model K80, this estimates a pairwise distance matrix from the sequence data

## tree building 
BbTree<-nj(BbDM)#neighbour joining 
ggtree(BbTree,layout="rectangular") + geom_tiplab()##generate the phylogenetic tree using ggtree(), arranged in rectangular form with tip labels 

#problems: 
###I got a stick, and I don't think it is suppose to be a stick. Does anyone know how to fix it?
  ## this is the tree result i get from BLAST: 
    #https://blast.ncbi.nlm.nih.gov/blast/treeview/treeView.cgi?request=page&blastRID=4V81BY0501R&queryID=sp|Q9ZZ64|&entrezLim=&ex=&exl=&exh=&ns=100
###Also we are suppose to show the phylogenetic relationship of species, instead we are currently showing geneIDs, this should be fixed. 
