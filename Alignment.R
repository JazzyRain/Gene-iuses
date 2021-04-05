library(BiocManager)
library(genbankr)
library(annotate)
library(ape)
library(ggtree)
library(rentrez)

mtID<-GBAccession("NC_002008.4") #using NCBI I seached for CO1, and canis lupus familiaris [organism], i find this accession number
mtRead<-readGenBank(mtID) #use the accession number to browse through the gene bank 
mtRead@genes[3] #since COX1, our reference gene of interest, is 3rd on the list, we can index it to see the exact range of the gene. The IRange given is 5349-6893. So let's snip to be just this range!

mt<-entrez_fetch(dbfrom="gene", id="NC_002008.4", db="nuccore", rettype="fasta") #download the COX1 reference genome from Genbank with the accession NC_002008.4"
mtSeq<-gsub("^>.*genome\\n([ATCG].*)","\\1",mt)#find the actual DNA sequence by capturing only the parts containing ATCG (after the word "genome", the last word in the name)
mtSeq<-gsub("\\n", "",mtSeq)#removing all the new line characters so that it is not cut in pieces
for(i in 1:nchar(mtSeq)){#for loop that runs through every character in the sequence
  if(i<5349){#for every time the number is smaller than 5349 (the position where the COX1 begins)
    COX1<-sub("\\w", "", mtSeq)## remove bo from beginning 
  } else if(i>6893){#for every time the number is larger than 6893 (the position where the COX1 ends)
    COX1<-sub("\\w$", "", COX1)##remove bp from the end 
  }
}

##Pairwise alignment
COX1BLAST<-blastSequences(COX1, as='data.frame', hitListSize=800, timeout=3000) #pairwise alignment search using our gene 
##multiple alignment
COXHitsDF<-data.frame(ID=COX1BLAST$Hit_def, Seq=COX1BLAST$Hsp_hseq, stringsAsFactors=FALSE)

COXHitsDF$ID

COXHitsDF$ID <- gsub("isolate.*|haplotype.*|mitochondrial.*|genome.*|mitochondrion.*", "", COXHitsDF$ID)## filter out different genes from the same species by removing the part that indicate the genetic difference, keeping only the species name
COXHitsDF <- COXHitsDF %>% distinct(ID, .keep_all = TRUE)## keep only distinct species and remove identical 


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