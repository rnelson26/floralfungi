# Analysis script
  ## based on the analysis_model.R script from Devin
# Steps: remove host/non-fungal sequences, remove putative contaminants, remove extraction blanks and controls, look at sequencing depth across samples, filter out samples with poor coverage based on seq depth, do preliminary ordinations and perMANOVA (of all pops, and subset without NZ seeds)
  ## input: phyloseq object
  ## output: ordination plots, PerMANOVA table

library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(ggthemes)
library(FUNGuildR)
library(decontam)
library(microViz)

sessionInfo()

#### Prepare data for analysis ####
(all.phy <- readRDS("output/compiled/allphy.rds"))
all.meta <- read.csv("data/RAN_McLaughlin_flower_mycobiomes_metadata_all.csv", row.names = 1)
sample_data(all.phy) <- all.meta

#### remove contaminants with decontam, compare to manual contaminant counts ####
# checking with the manual method first
control.phy <- subset_samples(all.phy, control_YN =="yes")
all.samps.phy <- subset_samples(all.phy, control_YN =="no")

control.phy.tab <- tax_table(control.phy) %>% data.frame %>% 
  transmute(
    OTU=taxa_names(control.phy),
    totSeqs=taxa_sums(all.phy %>% prune_taxa(taxa_names(control.phy),.)),
    totNegSeqs = taxa_sums(control.phy),
    pctSeqsNeg = round((100*totNegSeqs/totSeqs),4),
    genus=genus
  )

# decontam
df <- as.data.frame(sample_data(all.phy))
df$LibrarySize <- sample_sums(all.phy)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=control_YN)) + geom_point()

#identify negative controls
sample_data(all.phy)$is.neg <- sample_data(all.phy)$control_YN == "yes"

#run decontam with threshold 0.3 - this threshold was based on comparing results of decontam with threshold at 0.5 vs. manual method, most contaminants have a threshold below 0.3
contamdf.prev_ITS <- isContaminant(all.phy, method="prevalence", neg="is.neg", threshold=0.3)
table(contamdf.prev_ITS$contaminant) # 26 contaminants at treshold of 0.3
head(which(contamdf.prev_ITS$contaminant))

# show contaminant ASVs with associated taxonomy
contamdf.prev_ITS.list <- contamdf.prev_ITS %>% 
  filter(contaminant == TRUE)
contamdf.prev_ITS.list.taxo <-  merge(contamdf.prev_ITS.list, 
                                       as.data.frame(tax_table(all.phy)), 
                                       by="row.names")
contamdf.prev_ITS.list.taxo %>% write.csv("output/compiled/post_decontam_results.csv")

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(all.phy, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$control == "yes", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$control == "no", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev_ITS$contaminant)
ggplot(df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

contamdf.prev_ITS[contamdf.prev_ITS$contaminant == TRUE, ]  #26 true contaminants
# Prune TRUE contaminants in phyloseq objects
phy_v2 <- prune_taxa(!contamdf.prev_ITS$contaminant, all.phy)
phy_v2 #3262 taxa and 159 samples
phy_v2 %<>% subset_samples(., control_YN == "no")
phy_v2 %<>% prune_taxa(taxa_sums(.)>0, .) # 3248 taxa in 155 samples
phy_v2_asv <- as.data.frame(otu_table(phy_v2))

#### Look at sequencing depth ####
plotSeqDepth <- function(phy_v2,cutoff){
  pctPass <- phy_v2 %>% sample_data %>% data.frame %>% 
    mutate(SeqDepth=sample_sums(phy_v2),
           pass=SeqDepth>cutoff) %>% 
    #group_by(DataSet) %>%
    summarise(nSamps=n(),pctPass=round(100*sum(pass)/n(),1))
  phy_v2 %>% sample_data %>% data.frame %>%  
    mutate(SeqDepth=sample_sums(phy_v2)) %>% 
    .[order(.$SeqDepth),] %>%
    mutate(Index=factor(seq(nrow(.)))) %>%
    ggplot(aes(x=Index, y=SeqDepth,color=SeqDepth<cutoff)) + 
    geom_point(position=position_jitter(width=50),alpha=0.5,shape=19) +
    geom_text(aes(x=-Inf,y=Inf,label=paste0(nSamps," samples")),data=pctPass,inherit.aes = F, hjust=-0.1, vjust=1.5) +
    geom_text(aes(x=-Inf,y=Inf,label=paste0(pctPass, "% > ",cutoff, " seqs")),data=pctPass,inherit.aes = F, hjust=-0.1, vjust=3) +
    scale_color_calc(name=paste0("Sequencing depth > ",cutoff,"\n")) +
    #facet_wrap(~DataSet,scales="free") +
    theme_classic() + 
    theme(legend.position = "none",axis.text.x = element_blank(),axis.ticks.x = element_blank()) 
}
plotSeqDepth(phy_v2,2000) # 96.2 % of samples above this cutoff (6 samples fail)
minDepth <- 2000 # this minDepth is based off the cutoffs tested with the plot

# Remove samples below sequencing depth cutoff
(ITS.phy <- prune_samples(sample_sums(phy_v2)>minDepth,phy_v2)) # 3262 taxa in 152 samples

# Remove low abundance ASVs by prevalence - Only keep ASVs present in at least 2.5% of samples
ITS.abund.phy <- ITS.phy %>% 
  tax_filter(min_prevalence = 2.5, prev_detection_threshold = 10) # 596 taxa in 153 samples
ITS.abund.asv <- otu_table(ITS.abund.phy) %>% data.frame()

# convert objects to relative abundance
ITS.rel.phy <- transform_sample_counts(ITS.phy, function(x){(x/sum(x))*100}) # we will use this for alpha diversity analyses
ITS.abund.rel.phy <- transform_sample_counts(ITS.abund.phy, function(x){(x/sum(x))*100}) # we will use this for ordination, taxonomic visualizations

# save cleaned phy object and relevant information
dir.create("output/analysis")
saveRDS(ITS.rel.phy, "output/analysis/clean.ITS.rel.phy.rds")
saveRDS(ITS.abund.rel.phy, "output/analysis/clean.ITS.abund.rel.phy.rds")
otu_table(ITS.rel.phy) %>% write.csv("output/analysis/ITS.rel.ASV.table.csv")
tax_table(ITS.rel.phy) %>% write.csv("output/analysis/ITS.rel.ASVtaxonomy.csv")
otu_table(ITS.abund.rel.phy) %>% write.csv("output/analysis/ITS.abund.rel.ASV.table.csv")
tax_table(ITS.abund.rel.phy) %>% write.csv("output/analysis/ITS.abund.rel.ASVtaxonomy.csv")
