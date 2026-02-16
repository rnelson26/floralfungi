# Denoising 16S samples using DADA2
#for a more detailed tutorial see https://benjjneb.github.io/dada2/tutorial.html
# might not need most of this script, since Novogene trimmed, filtered and merged the reads already. Just need to compile into a sequence table. 

#### load packages, data, set file paths ####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("dada2", force = T)

library(tidyverse)
library(magrittr)
library(dada2)
library(ShortRead)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(stringr)
library(readxl)
library(tibble)

#path for output
out.path <- "output/dada"

#Identify paths to trimmed fwd and rev reads
in.path <- file.path("data/NovaSeq_files_round2/result_X202SC25014024-Z01-F001/00.CleanData_trimmed/")
path.seq <- in.path %>% list.files(pattern="effective.fastq.gz",full.names = T) %>% sort 
#path.rev <- in.path %>% list.files(pattern=".R2.fq",full.names = T) %>% sort 

#make file names / paths for trimmed and filtered files
filt.path <- file.path(out.path, "filt")
filt.seq <- path.seq %>% gsub(in.path,filt.path,.)
#filt.rev <- path.rev %>% gsub(in.path,filt.path,.)

#preview read quality of 12 randomly selected samples before trimming
qual.samps <- sample(1:length(path.seq),6)
plotQualityProfile(path.seq[qual.samps]) 
#plotQualityProfile(path.rev[qual.samps])

#### Filter and trim, dereplicate identical sequences ####
trim <- filterAndTrim(path.seq, filt.seq,
                      maxN=0, maxEE=c(2), truncQ=2,
                      multithread=TRUE)

#Update list of trimmed file paths to exclude samples with no reads passing filters
filt.seq <- list.files(filt.path, pattern = "effective.fastq.gz",full.names = T)
#filt.rev <- list.files(filt.path, pattern = ".R2.fq",full.names = T)

#Check quality of trimmed and filtered reads
qual.samps <- sample(1:length(filt.seq),12)
plotQualityProfile(filt.seq[qual.samps])
#plotQualityProfile(filt.rev[qual.samps]) 

#dereplicate identical sequences
derep.seq <- derepFastq(filt.seq, verbose=F)
#derep.rev <- derepFastq(filt.rev, verbose=F) 

#Trim names of derep objects to just sample names
names(derep.seq) %<>% gsub("effective.fastq.gz","",.)
#names(derep.rev) %<>% gsub(".R2.fq","",.)

#### Learn errors, denoise, merge sequences into table ####
err.seq <- learnErrors(filt.seq, multithread=TRUE, nbases = 3e08) 
plotErrors(err.seq, nominalQ=TRUE)
#err.rev <- learnErrors(filt.rev, multithread=TRUE, nbases = 3e08)
#plotErrors(err.rev, nominalQ=TRUE)

#Denoise
dada.seq <- dada(derep.seq, err=err.seq, multithread=TRUE)
#dada.rev <- dada(derep.rev, err=err.rev, multithread=TRUE)

#Merge reads - need to fix this since the reads have been trimmed and merged
merged <- mergePairs(dada.fwd, derep.fwd, dada.rev, derep.rev, trimOverhang = T)

#Make sequence table
seqtab <- dada.seq %>% makeSequenceTable

#### Remove chimeras, make summary report ####
seqtab.nonChimeras <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

#Make summary report of what happened
getN <- function(x) sum(getUniques(x))
trim.summary <- trim %>% data.frame %>% rownames_to_column("Sample") 
trim.summary$Sample %<>% str_sub(., 1, 12)
track <- cbind(sapply(dada.seq, getN), 
               #sapply(dada.rev, getN), 
               #sapply(merged, getN),
               rowSums(seqtab.nonChimeras)) %>% 
  data.frame %>% rownames_to_column("Sample") 
track %<>% left_join(trim.summary,.)
colnames(track) <- c("sample", "input", "filtered", "denoised", "ChimeraFiltered")
file.path(out.path,"dadaSummary2.csv") %>% write.csv(track,.,row.names = F)

#### Save output ####
file.path(out.path,"seqTab2.rds") %>% saveRDS(seqtab.nonChimeras,.)

