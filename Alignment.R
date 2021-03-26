library(BiocManager)
library(genbankr)
library(annotate)
library(ape)
library(ggtree)
library(rentrez)

mtID<-GBAccession("NC_002008.4") #using NCBI I seached for CO1, and canis lupus familiaris [organism], i find this accession number
mtRead<-readGenBank(BbID) #use the accession number to browse through the gene bank 
mtRead@genes[3] #since COX1, our reference gene of interest, is 3rd on the list, we can index it to see the exact range of the gene. The IRange given is 5349-6893. So let's snip to be just this range!

mt<-entrez_fetch(dbfrom="gene", id="NC_002008.4", db="nuccore", rettype="fasta")
mtSeq<-gsub("^>.*genome\\n([ATCG].*)","\\1",mt)
mtSeq<-gsub("\\n", "",mtSeq)

for(i in 1:nchar(mtSeq)){
  if(i<5349){
    COX1<-sub("\\w", "", mtSeq)
  } else if(i>6893){
    COX1<-sub("\\w$", "", COX1)
  }
}

##Pairwise alignment
COX1BLAST<-blastSequences(COX1, as='data.frame', hitListSize=40, timeout=600) #pairwise alignment search using our gene 
##multiple alignment
COXHitsDF<-data.frame(ID=COX1BLAST$Hit_accession, Seq=COX1BLAST$Hsp_hseq, stringsAsFactors=FALSE)


#muscle+ alignemnt inspection+ generate distance matrix 
COXHitsSeqs<-read.GenBank(COX1BLAST$Hit_accession[1:40])
COXHitsDNA<-sapply(COXHitsDF$Seq, strsplit, split="") #convert the object to one with seperate column for each basepair using strsplit()and sapply() to split the DA for each row
names(COXHitsDNA)<-paste(1:nrow(COXHitsDF), COXHitsDF$ID, sep="_") #give each sequence unique names.
COXHitsDNA<-as.DNAbin(COXHitsDNA) #convert it into a DNAbin object 
COXAlign<-muscle(COXHitsDNA, quiet=F) #run muscle on DNAbin object, allowing alignment 
checkAlignment(COXAlign, what=3) #to inspect alignments and determine Shannon indexes relatives to the sequence position, providing cues to diversity
COXDM<-dist.dna(COXAlign, model="K81") #use distance model K80, this estimates a pairwise distance matrix from the sequence data

# tree building 
COXTree<-nj(COXDM) #neighbour joining 
ggtree(COXTree, layout="rectangular")+geom_tiplab() #generate the phylogenetic tree using ggtree(), arranged in rectangular form with tip labels 
