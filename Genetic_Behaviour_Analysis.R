library(BiocManager)
library(dplyr)
library(genbankr)
library(rentrez)
library(annotate)
library(ape)
library(reshape2)
library(ggplot2)
library(ggtree)


##CADPS2 Alignment
#pulling CADPS2 sequence from genbank
CAD_Lab <- GBAccession("XM_038687133.1")
#reading the genbank object
CAD_LabGBk <- readGenBank(CAD_Lab)

#blast searching the CADPS2 sequence
BLASTSearch_CAD <- blastSequences(paste(CAD_LabGBk@sequence), as = 'data.frame', hitListSize = 40, timeout = 600)

#making a simple dataframe with just ID numbers and sequences
CADHitsDF <- data.frame(ID=BLASTSearch_CAD$Hit_def, Seq=BLASTSearch_CAD$Hsp_hseq, stringsAsFactors = FALSE)

#changing ID to just species name with regular expressions
CADHitsDF$ID <- gsub("PREDICTED: (.+) calcium.+", "\\1", CADHitsDF$ID)
CADHitsDF$ID <- gsub("PREDICTED: (.+) Ca.+", "\\1", CADHitsDF$ID)

#making sure theres no replicate hits from the same species
CADHitsDF <- CADHitsDF %>% distinct(ID, .keep_all = TRUE)

#turning the dataframe into a DNABin object by splitting the sequences so each nucleotide is its own string
CADHits <- sapply(CADHitsDF$Seq, strsplit, split="")

#adding species as the name
names(CADHits) <- CADHitsDF$ID

#then converting to a DNABin object
CADHits <- as.DNAbin(CADHits)

#Using Multiple Sequence Comparison by Log-Expectation (MUSCLE) to align the Canis lupus DRD4 sequence with the top blast searches
CAD_align <- muscle(CADHits)

#generating a distance matrix from the alignment
CAD_DM <- dist.dna(CAD_align, model = "K80")

#making a neighbour joining tree from the distance matrix
CAD_Tree <- njs(CAD_DM)

#visualizing the tree
CAD_treeim <- ggtree(CAD_Tree, layout = "rectangular")+
  geom_tiplab()

#exporting the tree as a pdf
ggsave(filename = "CADtree.pdf", CAD_treeim, device = "pdf", width = 20, height = 4, units = "in" , limitsize = FALSE)

#creating a tree with no branch lengths
CAD_nodist <- ggtree(CAD_Tree, layout = "rectangular", branch.length = "none")+
  geom_tiplab()+
  xlim(NA, 24)

#exporting the no branch length tree as a pdf
ggsave(filename = "CADtree_nodist.pdf", CAD_nodist, device = "pdf", width = 5, height = 4, units = "in" , limitsize = FALSE)

##DRD4 Alignment

#Pulling DRD4 Sequence from GenBank and reading it
DRD4_Lab <- GBAccession("XM_038424041.1")
DRD4_LabGBk <- readGenBank(DRD4_Lab)

#Blast search to find close sequences to Canis lupus familiaris DRD4 Sequence
BLASTSearch_DRD4 <- blastSequences(paste(DRD4_LabGBk@sequence), as = 'data.frame', hitListSize = 20, timeout = 600)

#making a simple dataframe with just ID numbers and sequences
DRD4HitsDF <- data.frame(ID=BLASTSearch_DRD4$Hit_def, Seq=BLASTSearch_DRD4$Hsp_hseq, stringsAsFactors = FALSE)

#changing ID to just species name with regular expressions
DRD4HitsDF$ID <- gsub("PREDICTED: (.+) dopamine.+", "\\1", DRD4HitsDF$ID)
DRD4HitsDF$ID <- gsub("(.+)dopamine.+", "\\1", DRD4HitsDF$ID)

#making sure theres no replicate hits from the same species
DRD4HitsDF <- DRD4HitsDF %>% distinct(ID, .keep_all = TRUE)

#turning the dataframe into a DNABin object by splitting the sequences so each nucleotide is its own string
DRD4Hits <- sapply(DRD4HitsDF$Seq, strsplit, split="")

#adding species as the name
names(DRD4Hits) <- DRD4HitsDF$ID

#then converting to a DNABin object
DRD4Hits <- as.DNAbin(DRD4Hits)

#Using Multiple Sequence Comparison by Log-Expectation (MUSCLE) to align the Canis lupus DRD4 sequence with the top blast searches
DRD4_align <- muscle(DRD4Hits)

#using the alignment to create a DNA distance matrix
DRD4_DM <- dist.dna(DRD4_align, model = "K80")

#using the distance matrix to create a neighbour joining tree
DRD4_Tree <- njs(DRD4_DM)

#generating the tree image as an object
Dr_treeim <- ggtree(DRD4_Tree, layout = "rectangular")+
  geom_tiplab()

#exporting the tree image as a pdf
ggsave(filename = "DRD4tree.pdf", Dr_treeim, device = "pdf", width = 20, height = 4, units = "in" , limitsize = FALSE)

#generating a tree that doesnt make the branch lengths a factor
DR_nodist <- ggtree(DRD4_Tree, branch.length = "none")+
  geom_tiplab()+
  xlim(NA, 24)

#exporting the tree without distance
ggsave(filename = "DRD4tree_nodist.pdf", DR_nodist, device = "pdf", width = 8, height = 4, units = "in" , limitsize = FALSE)

##TH Alignment

#Pulling TH Sequence from GenBank
THLab <- GBAccession("AB097058.1")

#reading the genbank object
THLabGBk <- readGenBank(THLab)

#running a blast search with the TH sequence
BLASTSearch_TH <- blastSequences(paste(THLabGBk@sequence), as = 'data.frame', hitListSize = 40, timeout = 600)

#making a simple dataframe with just ID names and sequences
THHitsDF <- data.frame(ID=BLASTSearch_TH$Hit_def, Seq=BLASTSearch_TH$Hsp_hseq, stringsAsFactors = FALSE)

#changing ID to just species name with regular expressions
THHitsDF$ID <- gsub("PREDICTED: (.+) tyrosine.+", "\\1", THHitsDF$ID)
THHitsDF$ID[1] <- "Canis lupus familiaris"

#making sure theres no replicate hits from the same species
THHitsDF <- THHitsDF %>% distinct(ID, .keep_all = TRUE)

#turning the dataframe into a DNABin object by splitting the sequences so each nucleotide is its own string
THHits <- sapply(THHitsDF$Seq, strsplit, split="")

#adding species as the name
names(THHits) <- THHitsDF$ID

#then converting to a DNABin object
THHits <- as.DNAbin(THHits)

#Using Multiple Sequence Comparison by Log-Expectation (MUSCLE) to align the Canis lupus DRD4 sequence with the top blast searches
TH_align <- muscle(THHits)

#making a distance matrix based on the alignment
TH_DM <- dist.dna(TH_align)

#making a neighbour joining tree based on distance matrix
TH_Tree <- njs(TH_DM)

#visualizing the tree
TH_treeim <- ggtree(TH_Tree)+
  geom_tiplab()
#saving the tree image
ggsave(filename = "THtree.pdf", TH_treeim, device = "pdf", width = 20, height = 4, units = "in" , limitsize = FALSE)

#creating a tree without branch dists
TH_nodist <- ggtree(TH_Tree, branch.length = "none")+
  geom_tiplab()+
  xlim(NA, 24)
#saving no branch dist tree
ggsave(filename = "THtree_nodist.pdf", TH_nodist, device = "pdf", width = 8, height = 4, units = "in" , limitsize = FALSE)
