# Collating UpSet plot data

#### loading packages ####
library(tidyverse)
library(phyloseq)
library(UpSetR)
library(ComplexUpset)
library(microbiomeutilities)
library(microViz)
library(ggpubr)

#### upset plot by host species ####
# merging samples by plant species
variable1 = as.character(get_variable(ITS.rel.phy, "Plant_species"))
sample_data(ITS.rel.phy)$Plant_species <- mapply(paste0, variable1, collapse = "_")
merge_samples(ITS.rel.phy, "Plant_species")
merged_ITS_all <- merge_samples(ITS.rel.phy, "Plant_species")
#Fix continuous variable to discrete value
sample_data(merged_ITS_all)$Plant_species <- factor(sample_names(merged_ITS_all),
                                                    levels = c("Goldfield", "Clover","Vetch"))
#Extract count datasets from phyloseq
count_ITS_all <- t(otu_table(merged_ITS_all))
#Transform in binary matrix
count_ITS_all[count_ITS_all> 0] <- 1
#Convert to dataframe
df_ITS_all <- as.data.frame(count_ITS_all)
#Read abundance per ASV
rel.abu <- data.frame(taxa_sums(ITS.rel.phy))
rel.abu$log <- log10(rel.abu$taxa_sums.ITS.rel.phy.)
#Merge ASV dataframe and ASV abundance
df_ITS_all_upset <-merge.data.frame(df_ITS_all, rel.abu, by="row.names",
                                    all.x=TRUE)
#Append taxonomy to upset dataframe
df_ITS_taxa <- tax_table(ITS.rel.phy) %>% data.frame() %>% rownames_to_column()
names(df_ITS_taxa)[names(df_ITS_taxa)=="rowname"] <- "Row.names"

df_ITS_all_upset2 <- full_join(df_ITS_all_upset, df_ITS_taxa, by="Row.names")

#Basic UpSet plot will all comparisons
Upset_ITS_all <- UpSetR::upset(df_ITS_all_upset2,
                               nsets = 3, #set to number of columns to compare, here 3 Plant species
                               nintersects = NA,
                               order.by = "freq",
                               mainbar.y.label = "Number of ASVs",
                               sets.x.label="total plant ASV richness",
                               decreasing ="TRUE")
Upset_ITS_all

#### upset plot of just clover vs. vetch (overall) ####
# subsetting the dataset into a new phyloseq object of just clover and vetch samples
TV.rel.phy <- subset_samples(ITS.rel.phy, Plant_species == "Clover"|Plant_species == "Vetch")
TV.rel.phy <- prune_taxa(taxa_sums(TV.rel.phy)>0, TV.rel.phy)
# setting plant species as the comparison variable again, merging samples by species
variable2 = as.character(get_variable(TV.rel.phy, "Plant_species"))
sample_data(TV.rel.phy)$Plant_species <- mapply(paste0, variable2, collapse = "_")
merge_samples(TV.rel.phy, "Plant_species")
merged_TV_all <- merge_samples(TV.rel.phy, "Plant_species")
#Fix continuous variable to discrete value
sample_data(merged_TV_all)$Plant_species <- factor(sample_names(merged_TV_all),
                                                   levels = c("Clover","Vetch"))
#Extract count datasets from phyloseq
count_TV_all <- t(otu_table(merged_TV_all))
#Transform in binary matrix
count_TV_all[count_TV_all> 0] <- 1
#Convert to dataframe
df_TV_all <- as.data.frame(count_TV_all)
#Read abundance per ASV
TV.rel.abu <- data.frame(taxa_sums(TV.rel.phy))
TV.rel.abu$log <- log10(TV.rel.abu$taxa_sums.TV.rel.phy.)
#Merge ASV dataframe and ASV abundance
df_TV_all_upset <-merge.data.frame(df_TV_all, TV.rel.abu, by="row.names",
                                   all.x=TRUE)
#Append taxonomy to upset dataframe
df_TV_taxa <- tax_table(TV.rel.phy) %>% data.frame() %>% rownames_to_column()
names(df_TV_taxa)[names(df_TV_taxa)=="rowname"] <- "Row.names"

df_TV_all_upset2 <- full_join(df_TV_all_upset, df_TV_taxa, by="Row.names")

#Basic UpSet plot will all comparisons
Upset_TV_all <- UpSetR::upset(df_TV_all_upset2,
                              nsets = 2, #set to number of columns to compare, here 2 Plant species
                              nintersects = NA,
                              order.by = "freq",
                              mainbar.y.label = "Number of ASVs",
                              sets.x.label="Total plant\nASV richness",
                              decreasing ="TRUE")
Upset_TV_all

#### comparing vetch to clover at each distance bin within invasion levels ####
# subsetting the data
TV.low.rel.phy <- subset_samples(TV.rel.phy, Site_invasion_level == "Low")
TV.low.rel.phy <- prune_taxa(taxa_sums(TV.low.rel.phy)>0, TV.low.rel.phy)
TV.med.rel.phy <- subset_samples(TV.rel.phy, Site_invasion_level == "Medium")
TV.med.rel.phy <- prune_taxa(taxa_sums(TV.med.rel.phy)>0, TV.med.rel.phy)
TV.high.rel.phy <- subset_samples(TV.rel.phy, Site_invasion_level == "High")
TV.high.rel.phy <- prune_taxa(taxa_sums(TV.high.rel.phy)>0, TV.high.rel.phy)

# extracting the metadata from each object for binning distance levels
TV.low.rel.meta <- sample_data(TV.low.rel.phy) %>% data.frame()
TV.med.rel.meta <- sample_data(TV.med.rel.phy) %>% data.frame()
TV.high.rel.meta <- sample_data(TV.high.rel.phy) %>% data.frame()

TV.low.rel.meta[is.na(TV.low.rel.meta)] <- "Vetch"
TV.med.rel.meta[is.na(TV.med.rel.meta)] <- "Vetch"
TV.high.rel.meta[is.na(TV.high.rel.meta)] <- "Vetch"


# mutating distance from vetch for each object
TV.low.rel.meta2 <- TV.low.rel.meta %>%
  mutate(Distance = case_when(
    dist_from_vetch_m == "0" ~ "0-1",  
    dist_from_vetch_m == "1" ~ "0-1", 
    dist_from_vetch_m == "50" ~ "50",  
    dist_from_vetch_m == "Vetch" ~ "Vetch",
    TRUE ~ "5-10"
  ))

TV.med.rel.meta2 <- TV.med.rel.meta %>%
  mutate(Distance = case_when(
    dist_from_vetch_m == "0" ~ "0-1",  
    dist_from_vetch_m == "1" ~ "0-1", 
    dist_from_vetch_m == "50" ~ "50",  
    dist_from_vetch_m == "Vetch" ~ "Vetch",
    TRUE ~ "5-10"
  ))

TV.high.rel.meta2 <- TV.high.rel.meta %>%
  mutate(Distance = case_when(
    dist_from_vetch_m == "0" ~ "0-1",  
    dist_from_vetch_m == "1" ~ "0-1", 
    dist_from_vetch_m == "50" ~ "50",  
    dist_from_vetch_m == "Vetch" ~ "Vetch",
    TRUE ~ "5-10"
  ))

# setting "Distance" as a factor in each object
TV.low.rel.meta2$Distance <- factor(TV.low.rel.meta2$Distance, 
                                    levels = c("Vetch","0-1","5-10","50"))
TV.med.rel.meta2$Distance <- factor(TV.med.rel.meta2$Distance, 
                                    levels = c("Vetch","0-1","5-10","50"))
TV.high.rel.meta2$Distance <- factor(TV.high.rel.meta2$Distance, 
                                    levels = c("Vetch","0-1","5-10","50"))

# replacing phyloseq metadata with new data frames
sample_data(TV.low.rel.phy) <- TV.low.rel.meta2
sample_data(TV.med.rel.phy) <- TV.med.rel.meta2
sample_data(TV.high.rel.phy) <- TV.high.rel.meta2

# setting plant species as the comparison variable again, merging samples by species
variable3 = as.character(get_variable(TV.low.rel.phy, "Distance"))
sample_data(TV.low.rel.phy)$Distance <- mapply(paste0, variable3, collapse = "_")
merged_TV_low_all <- merge_samples(TV.low.rel.phy, "Distance")

variable4 = as.character(get_variable(TV.med.rel.phy, "Distance"))
sample_data(TV.med.rel.phy)$Distance <- mapply(paste0, variable4, collapse = "_")
merged_TV_med_all <- merge_samples(TV.med.rel.phy, "Distance")

variable5 = as.character(get_variable(TV.high.rel.phy, "Distance"))
sample_data(TV.high.rel.phy)$Distance <- mapply(paste0, variable5, collapse = "_")
merged_TV_high_all <- merge_samples(TV.high.rel.phy, "Distance")

#Extract count datasets from phyloseq
count_TV_low_all <- t(otu_table(merged_TV_low_all))
count_TV_med_all <- t(otu_table(merged_TV_med_all))
count_TV_high_all <- t(otu_table(merged_TV_high_all))

#Transform in binary matrix
count_TV_low_all[count_TV_low_all> 0] <- 1
count_TV_med_all[count_TV_med_all> 0] <- 1
count_TV_high_all[count_TV_high_all> 0] <- 1

#Convert to dataframe
df_TV_low_all <- as.data.frame(count_TV_low_all)
df_TV_med_all <- as.data.frame(count_TV_med_all)
df_TV_high_all <- as.data.frame(count_TV_high_all)

#Append taxonomy to upset dataframe
df_TV_low_taxa <- tax_table(TV.low.rel.phy) %>% data.frame() %>% rownames_to_column()
df_TV_low_all %<>% rownames_to_column()
df_TV_low_upset <- full_join(df_TV_low_all, df_TV_low_taxa, by="rowname")

df_TV_med_taxa <- tax_table(TV.med.rel.phy) %>% data.frame() %>% rownames_to_column()
df_TV_med_all %<>% rownames_to_column()
df_TV_med_upset <- full_join(df_TV_med_all, df_TV_med_taxa, by="rowname")

df_TV_high_taxa <- tax_table(TV.high.rel.phy) %>% data.frame() %>% rownames_to_column()
df_TV_high_all %<>% rownames_to_column()
df_TV_high_upset <- full_join(df_TV_high_all, df_TV_high_taxa, by="rowname")

# make the upset plots for each distance bin within each data subset, write out percentages
df_TV_low_upset %>%
  select(rowname, 'Vetch', '0-1') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 112 ASVs are shared (26.9%), 304 are unique to clover (73.1%)

df_TV_low_upset %>%
  select(rowname, 'Vetch', '5-10') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 118 ASVs are shared (21%), 444 are unique to clover (79%)

df_TV_low_upset %>%
  select(rowname, 'Vetch', '50') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 94 ASVs are shared (24%), 297 are unique to clover (76%)

df_TV_med_upset %>%
  select(rowname, 'Vetch', '0-1') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 59 ASVs are shared (17.5%), 278 are unique to clover (82.5%)

df_TV_med_upset %>%
  select(rowname, 'Vetch', '5-10') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 69 ASVs are shared (15%), 392 are unique to clover (85%)

#df_TV_med_upset %>%
#  select(rowname, 'Vetch', '50') %>%
#  UpSetR::upset(., nsets = 3, order.by = "freq") 
# skipping this distance because there are only two clover samples

df_TV_high_upset %>%
  select(rowname, 'Vetch', '0-1') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 125 ASVs are shared (24%), 395 are unique to clover (76%)

df_TV_high_upset %>%
  select(rowname, 'Vetch', '5-10') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 106 ASVs are shared (22%), 375 are unique to clover (78%)

df_TV_high_upset %>%
  select(rowname, 'Vetch', '50') %>%
  UpSetR::upset(., nsets = 3, order.by = "freq") # 118 ASVs are shared (28%), 303 are unique to clover (72%)

#### Making a plot summarizing the UpSet results ####
invasion <- c("Low invasion","Low invasion","Low invasion",
              "Low invasion","Low invasion","Low invasion",
              "Medium invasion","Medium invasion","Medium invasion","Medium invasion",
              "High invasion","High invasion","High invasion",
              "High invasion","High invasion","High invasion")
distance <- c('0-1','0-1','5-10','5-10','50','50',
              '0-1','0-1','5-10','5-10',
              '0-1','0-1','5-10','5-10','50','50')
ASV_type <- c("Shared","Unique","Shared","Unique","Shared","Unique",
              "Shared","Unique","Shared","Unique",
              "Shared","Unique","Shared","Unique","Shared","Unique")
Proportion <- c(26.9,73.1,21,79,24,76,
                17.5,82.5,15,85,
                24,76,22,78,28,72)
Upset_vals <- data.frame(invasion,distance,ASV_type,Proportion)
Upset_vals$invasion <- factor(Upset_vals$invasion,
                              levels = c("Low invasion","Medium invasion","High invasion"))

ASV.labs <- c("Shared with Vetch","Unique to Clover") 
names(ASV.labs) <- c("Shared","Unique")

Upset_summary_plot <- Upset_vals %>%
  filter(ASV_type == "Unique") %>%
  ggplot(aes(x=distance,y=Proportion))+
  geom_col(position = "stack")+
  facet_grid(~invasion, scales = "free") +
  theme_light()+
  coord_cartesian(ylim = c(40, 90))+
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=18),
        strip.text = element_text(size=20))+
  labs(x="Distance from Vetch (m)",y="Percentage of ASVs\nunique to Clover",fill="ASV type")
Upset_summary_plot

#### Making taxonomic composition plots of shared taxa ####
df_TV_low_upset[df_TV_low_upset== 0] <- "FALSE"
df_TV_low_upset[df_TV_low_upset== 1] <- "TRUE"
df_TV_med_upset[df_TV_med_upset==0] <- "FALSE"
df_TV_med_upset[df_TV_med_upset==1] <- "TRUE"
df_TV_high_upset[df_TV_high_upset==0] <- "FALSE"
df_TV_high_upset[df_TV_high_upset==1] <- "TRUE"

# low invasion sites
TV.low.rel.taxa <- tax_table(TV.low.rel.phy) %>% data.frame() %>% rownames_to_column()
df_TV_low_upset2 <- df_TV_low_upset %>%
  select(rowname, '0-1', '5-10', '50', Vetch) # removing duplicate taxonomy data
colnames(df_TV_low_upset2)[colnames(df_TV_low_upset2) == '0-1'] <- 'bin01' # renaming columns for correct subsetting (phyloseq doesn't like number strings as characters)
colnames(df_TV_low_upset2)[colnames(df_TV_low_upset2) == '5-10'] <- 'bin510'
colnames(df_TV_low_upset2)[colnames(df_TV_low_upset2) == '50'] <- 'bin50'

TV.low.rel.taxa2 <- merge(TV.low.rel.taxa, df_TV_low_upset2, by="rowname")
TV.low.rel.taxa2 <- column_to_rownames(TV.low.rel.taxa2, var = "rowname")
tax_table(TV.low.rel.phy) <- as.matrix(TV.low.rel.taxa2)
TV.low.shared.rel.phy <- subset_taxa(TV.low.rel.phy,bin01=="TRUE"&bin510=="TRUE"&
                                       bin50=="TRUE"&Vetch=="TRUE")
TV.low.shared.rel.phy <- prune_samples(sample_sums(TV.low.shared.rel.phy)>0, TV.low.shared.rel.phy)

TV.low.shared.merged.phy <- TV.low.shared.rel.phy %>%
  aggregate_top_taxa2(., 20, "genus")

brewerPlus <- distinct_palette()
scales::show_col(brewerPlus)

low.shared <- plot_bar(TV.low.shared.merged.phy, x="Sample", fill = "genus") + 
  geom_bar(aes(fill=genus), stat="identity", position="stack", color="gray33") +
  facet_wrap(Site~Plant_species, scales = "free_x", drop = TRUE)+ 
  theme_light()+
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=10), 
        axis.text.x = element_blank())+
       # axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1))+  
  guides(fill=guide_legend(ncol = 3)) +
  labs(title="Low invasion sites",y="Relative abundance (%)", fill="Genus")+
  scale_fill_manual(values = brewerPlus)
low.shared

# medium invasion sites
TV.med.rel.taxa <- tax_table(TV.med.rel.phy) %>% data.frame() %>% rownames_to_column()
df_TV_med_upset2 <- df_TV_med_upset %>%
  select(rowname, '0-1', '5-10', '50', Vetch) # removing duplicate taxonomy data
colnames(df_TV_med_upset2)[colnames(df_TV_med_upset2) == '0-1'] <- 'bin01' # renaming columns for correct subsetting (phyloseq doesn't like number strings as characters)
colnames(df_TV_med_upset2)[colnames(df_TV_med_upset2) == '5-10'] <- 'bin510'
colnames(df_TV_med_upset2)[colnames(df_TV_med_upset2) == '50'] <- 'bin50'

TV.med.rel.taxa2 <- merge(TV.med.rel.taxa, df_TV_med_upset2, by="rowname")
TV.med.rel.taxa2 <- column_to_rownames(TV.med.rel.taxa2, var = "rowname")
tax_table(TV.med.rel.phy) <- as.matrix(TV.med.rel.taxa2)
TV.med.shared.rel.phy <- subset_taxa(TV.med.rel.phy,bin01=="TRUE"&bin510=="TRUE"&
                                       bin50=="TRUE"&Vetch=="TRUE")
TV.med.shared.rel.phy <- prune_samples(sample_sums(TV.med.shared.rel.phy)>0, TV.med.shared.rel.phy)

TV.med.shared.merged.phy <- TV.med.shared.rel.phy %>%
  aggregate_top_taxa2(., 20, "genus")

brewerPlus <- distinct_palette()
scales::show_col(brewerPlus)

med.shared <- plot_bar(TV.med.shared.merged.phy, x="Sample", fill = "genus") + 
  geom_bar(aes(fill=genus), stat="identity", position="stack", color="gray33") +
  facet_wrap(Site~Plant_species, scales = "free_x", drop = TRUE)+ 
  theme_light()+
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=10), 
        axis.text.x = element_blank())+
        #axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1))+  
  guides(fill=guide_legend(ncol = 3)) +
  labs(title="Medium invasion sites",y="Relative abundance (%)", fill="Genus")+
  scale_fill_manual(values = brewerPlus)
med.shared

# high invasion sites
TV.high.rel.taxa <- tax_table(TV.high.rel.phy) %>% data.frame() %>% rownames_to_column()
df_TV_high_upset2 <- df_TV_high_upset %>%
  select(rowname, '0-1', '5-10', '50', Vetch) # removing duplicate taxonomy data
colnames(df_TV_high_upset2)[colnames(df_TV_high_upset2) == '0-1'] <- 'bin01' # renaming columns for correct subsetting (phyloseq doesn't like number strings as characters)
colnames(df_TV_high_upset2)[colnames(df_TV_high_upset2) == '5-10'] <- 'bin510'
colnames(df_TV_high_upset2)[colnames(df_TV_high_upset2) == '50'] <- 'bin50'

TV.high.rel.taxa2 <- merge(TV.high.rel.taxa, df_TV_high_upset2, by="rowname")
TV.high.rel.taxa2 <- column_to_rownames(TV.high.rel.taxa2, var = "rowname")
tax_table(TV.high.rel.phy) <- as.matrix(TV.high.rel.taxa2)
TV.high.shared.rel.phy <- subset_taxa(TV.high.rel.phy,bin01=="TRUE"&bin510=="TRUE"&
                                       bin50=="TRUE"&Vetch=="TRUE")
TV.high.shared.rel.phy <- prune_samples(sample_sums(TV.high.shared.rel.phy)>0, TV.high.shared.rel.phy)

TV.high.shared.merged.phy <- TV.high.shared.rel.phy %>%
  aggregate_top_taxa2(., 20, "genus")

brewerPlus <- distinct_palette()
scales::show_col(brewerPlus)

high.shared <- plot_bar(TV.high.shared.merged.phy, x="Sample", fill = "genus") + 
  geom_bar(aes(fill=genus), stat="identity", position="stack", color="gray33") +
  facet_wrap(Site~Plant_species, scales = "free_x", drop = TRUE)+ 
  theme_light()+
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=10), 
        axis.text.x = element_blank())+
        #axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1))+  
  guides(fill=guide_legend(ncol = 3)) +
  labs(title='High invasion sites',y="Relative abundance (%)", fill="Genus")+
  scale_fill_manual(values = brewerPlus)
high.shared

taxa.overlaps <- ggarrange(low.shared, med.shared, high.shared, nrow = 1)
taxa.overlaps

all.shared.merged.phy <- merge_phyloseq(TV.low.shared.merged.phy, TV.med.shared.merged.phy,
                                        TV.high.shared.merged.phy)
sample_data(all.shared.merged.phy)$Site_invasion_level <- factor(sample_data(all.shared.merged.phy)$Site_invasion_level,
                                                                 levels = c("Low","Medium","High"))

all.shared <- plot_bar(all.shared.merged.phy, x="Sample", fill = "genus") + 
  geom_bar(aes(fill=genus), stat="identity", position="stack", color="gray33") +
  facet_wrap(Site_invasion_level~Plant_species, scales = "free_x", 
             drop = TRUE, ncol=2)+ 
  theme_light()+
  theme(legend.position = "right", 
        legend.key.size = unit(0.5, 'cm'),legend.text = element_text(size=10), 
        axis.text.x = element_blank())+
  #axis.text.x = element_text(size=15, angle = 90, vjust = 0.5, hjust=1))+  
  guides(fill=guide_legend(ncol = 1)) +
  labs(y="Relative abundance (%)", fill="Genus")+
  scale_fill_manual(values = brewerPlus)
all.shared

#### Differential abundance analysis of shared taxa ####
# low invasion
output.low <- ancombc2(data = TV.low.shared.rel.phy, assay_name = "counts", tax_level = "genus",
                   fix_formula = "Plant_species", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,
                   group = "Plant_species", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(-1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE)),
                                        node = list(2, 2, 1),
                                        solver = "ECOS",
                                        B = 100)) # could add field to the fix_formula argument

# primary analysis (differential abundance in seeds)
tab_zero_low = output.low$zero_ind
print(head(tab_zero_low))
res_prim_low = output.low$res 

df_low_vetch = res_prim_low %>% # rel. abund. in vetch compared to Clover
  dplyr::select(taxon, ends_with("Vetch")) 
df_fig_low_vetch = df_low_vetch %>%
  dplyr::filter(diff_Plant_speciesVetch == 1) %>% 
  dplyr::arrange(desc(lfc_Plant_speciesVetch)) %>%
  dplyr::mutate(direct = ifelse(lfc_Plant_speciesVetch > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_Plant_speciesVetch == 1, "aquamarine3", "black"))

df_fig_low_vetch$taxon = factor(df_fig_low_vetch$taxon, levels = df_fig_low_vetch$taxon)
df_fig_low_vetch$direct = factor(df_fig_low_vetch$direct, 
                            levels = c("Positive LFC", "Negative LFC"))
# filter out anything that failed the pseudo-count sensitivity test
df_fig_low_vetch2 <- df_fig_low_vetch %>%
  filter(color == "aquamarine3") 
#df_fig_seed2$genus <- c("Shigella", "Kosakonia", "Acidovorax", "Bradyrhizobium")

fig_low_vetch = df_fig_low_vetch2 %>%
  ggplot(aes(x = taxon, y = lfc_Plant_speciesVetch)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Plant_speciesVetch - se_Plant_speciesVetch, 
                    ymax = lfc_Plant_speciesVetch + se_Plant_speciesVetch), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Genus", y = "Log fold change", 
       title = "Differential abundance in Vetch, low invasion") + 
  theme_light() + 
  theme(plot.title = element_text(size=20,hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15,angle = 60, hjust = 1))
fig_low_vetch

# medium invasion
output.med <- ancombc2(data = TV.med.shared.rel.phy, assay_name = "counts", tax_level = "genus",
                   fix_formula = "Plant_species", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,
                   group = "Plant_species", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(-1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE)),
                                        node = list(2, 2, 1),
                                        solver = "ECOS",
                                        B = 100)) # could add field to the fix_formula argument

# primary analysis (differential abundance in seeds)
tab_zero_med = output.med$zero_ind
print(head(tab_zero_med))
res_prim_med = output.med$res 

df_med_Vetch = res_prim_med %>%
  dplyr::select(taxon, ends_with("Vetch")) 
df_fig_med_Vetch = df_med_Vetch %>%
  dplyr::filter(diff_Plant_speciesVetch == 1) %>% 
  dplyr::arrange(desc(lfc_Plant_speciesVetch)) %>%
  dplyr::mutate(direct = ifelse(lfc_Plant_speciesVetch > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_Plant_speciesVetch == 1, "aquamarine3", "black"))

df_fig_med_Vetch$taxon = factor(df_fig_med_Vetch$taxon, levels = df_fig_med_Vetch$taxon)
df_fig_med_Vetch$direct = factor(df_fig_med_Vetch$direct, 
                            levels = c("Positive LFC", "Negative LFC"))
# filter out anything that failed the pseudo-count sensitivity test
df_fig_med_Vetch2 <- df_fig_med_Vetch %>%
  filter(color == "aquamarine3") 

fig_med_Vetch = df_fig_med_Vetch2 %>%
  ggplot(aes(x = taxon, y = lfc_Plant_speciesVetch)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Plant_speciesVetch - se_Plant_speciesVetch, 
                    ymax = lfc_Plant_speciesVetch + se_Plant_speciesVetch), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Genus", y = "Log fold change", 
       title = "Differential abundance in seeds") + 
  theme_light() + 
  theme(plot.title = element_text(size=20,hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15,angle = 60, hjust = 1))
fig_med_Vetch # no differentially abundant genera

# high invasion
output.high <- ancombc2(data = TV.high.shared.rel.phy, assay_name = "counts", tax_level = "genus",
                   fix_formula = "Plant_species", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,
                   group = "Plant_species", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(-1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE)),
                                        node = list(2, 2, 1),
                                        solver = "ECOS",
                                        B = 100)) # could add field to the fix_formula argument

# primary analysis (differential abundance in seeds)
tab_zero_high = output.high$zero_ind
print(head(tab_zero_high))
res_prim_high = output.high$res # differences for stigmas are the (intercept), since it's the reference group

df_high_Vetch = res_prim_high %>%
  dplyr::select(taxon, ends_with("Vetch")) 
df_fig_high_Vetch = df_high_Vetch %>%
  dplyr::filter(diff_Plant_speciesVetch == 1) %>% 
  dplyr::arrange(desc(lfc_Plant_speciesVetch)) %>%
  dplyr::mutate(direct = ifelse(lfc_Plant_speciesVetch > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_Plant_speciesVetch == 1, "aquamarine3", "black"))

df_fig_high_Vetch$taxon = factor(df_fig_high_Vetch$taxon, levels = df_fig_high_Vetch$taxon)
df_fig_high_Vetch$direct = factor(df_fig_high_Vetch$direct, 
                            levels = c("Positive LFC", "Negative LFC"))
# filter out anything that failed the pseudo-count sensitivity test
df_fig_high_Vetch2 <- df_fig_high_Vetch %>%
  filter(color == "aquamarine3") 

fig_high_Vetch = df_fig_high_Vetch2 %>%
  ggplot(aes(x = taxon, y = lfc_Plant_speciesVetch)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_Plant_speciesVetch - se_Plant_speciesVetch, 
                    ymax = lfc_Plant_speciesVetch + se_Plant_speciesVetch), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Genus", y = "Log fold change", 
       title = "Differential abundance in seeds") + 
  theme_light() + 
  theme(plot.title = element_text(size=20,hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15,angle = 60, hjust = 1))
fig_high_Vetch # no differentially abundant genera
