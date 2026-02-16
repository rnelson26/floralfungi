# Compile MiSeq library and sample data into a single wood phyloseq object
# see phyloseq tutorial at https://joey711.github.io/phyloseq/
  ## steps: clustering OTUs, filtering out host sequences (against P.trich. genome), filtering out contaminants (based on negative controls), filtering out rare sequences, combining OTU-by-sample matrix into the phyloseq object
  ## input: seqTab.rds, sample metadata (csv file)
  ## output: phyloseq object for analysis (rds file), OTU table (csv), taxon table (csv)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

#BiocManager::install("DECIPHER")

library(phyloseq)
library(magrittr)
library(dada2)
library(Biostrings)
library(tidyverse)
library(foreach)
library(DECIPHER)

sessionInfo()

# read in sample metadata (row.names refers to column one, where the sample IDs are)
meta <- read.csv("data/RAN_McLaughlin_flower_mycobiomes_metadata_round2.csv"
                 ,as.is=T,row.names = 1 ,sep=",")

# read in denoised sequence table
seqTab <- readRDS("output/dada/seqTab2.rds")
seqTab.df <- data.frame(seqTab) %>% rownames_to_column()
seqTab.df %>% write.csv("output/dada/seqTab2.csv")
#extract sequences from OTU table 
seqs <- getSequences(seqTab) %>% DNAStringSet #1756 sequences
#rename sequences
ASV.names <- paste0("ASV.",1:ncol(seqTab))
colnames(seqTab) <- ASV.names
names(seqs) <- ASV.names
#Create data frame for summary of processing
compile.summary <- data.frame(SampID=rownames(seqTab),
                              denoised=rowSums(seqTab))
# create scratch directory
dir.create("output/scratch")
# write temp file
writeXStringSet(seqs,file="output/scratch/tmp.fasta",width=1000)

# remove sequences that are shorter than 50 bp
seqs <- readDNAStringSet("output/scratch/tmp.fasta") %>% .[.@ranges@width > 50] # 1755 sequences
seqs %>% writeXStringSet(., file="output/scratch/ITS.tmp.fasta", width=1000)

#remove ASVs from seqTab that are too short
seqTab %<>% .[,names(seqs)] # 1755 sequences/ASVs in seqTab
#collapse any identical sequences 
colnames(seqTab) <- seqs %>% as.character %>% unname
seqTab %<>% collapseNoMismatch # 1728 sequences/ASVs remaining after collapsing 
# updating the sequence list so that it matches the number of ASVs in the trimmed ASV table
final.seqs <- dada2::getSequences(seqTab) %>% DNAStringSet()

#assign taxonomy with IdTaxa from the DECIPHER package (recommended by DADA2 developer, faster than assignTaxonomy) - A Murali et al. (2018) "IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences." Microbiome, doi:10.1186/s40168-018-0521-5.
load("data/UNITE_vAll_April2024.RData")
ids <- IdTaxa(final.seqs, trainingSet, strand = "top", processors = NULL, verbose = T)

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
## Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid.ITS <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid.ITS) <- ranks 
taxid.df <- data.frame(taxid.ITS) %>% rownames_to_column()
#rownames(taxid.ITS) <- getSequences(seqTab)

#set final ASV names
final.names <- paste0("ASV.",1:ncol(seqTab)) # different length from taxa

#final.seqs <- dada2::getSequences(seqTab) %>% DNAStringSet()
names(final.seqs) <- final.names
colnames(seqTab) <- final.names
rownames(taxid.ITS) <- final.names # number of ASVs don't match - removed duplicate/short sequences


#make phyloseq object
phy2 <- phyloseq(otu_table(seqTab,taxa_are_rows = F), 
                sample_data(meta),
                tax_table(taxid.ITS),
                refseq(final.seqs))
phy2.asv <- otu_table(phy2) %>% data.frame()

# remove non-fungal taxa from phyloseq object
min.reads <- 0
keep.samples <- sample_sums(phy2) > min.reads
phy2 %<>% subset_taxa(., kingdom=="Fungi") %>% prune_samples(keep.samples,.) %>% filter_taxa(function(x) {sum(x) > 0}, TRUE)

phy <- readRDS("output/compiled/phy.rds")
all.phy <- merge_phyloseq(phy, phy2)

#Save phyloseq object
dir.create("output/compiled")
saveRDS(phy2,"output/compiled/phy2.rds")
saveRDS(all.phy, "output/compiled/allphy.rds")

#Save components for possible manual inspection
otu_table(phy2) %>% write.csv("output/compiled/ASV.table2.csv")
tax_table(phy2) %>% write.csv("output/compiled/taxonomy.table2.csv")
otu_table(all.phy) %>% write.csv("output/compiled/ASV.table.all.csv")
tax_table(all.phy) %>% write.csv("output/compiled/taxonomy.table.all.csv")
