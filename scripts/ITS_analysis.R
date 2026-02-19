# ITS analysis

### code by RAN and GEB

### Date of most recent update to code: 2-19-2026

#### loading relevant packages ####
library(patchwork)
library(tidyverse)
library(phyloseq)
library(vegan)
library(plotly)
library(viridis)
library(ANCOMBC)
library(microbiome)
library(microbiomeutilities)
library(microViz)
library(FSA)
library(Rmisc)
library(ggpubr)
library(UpSetR)
library(ComplexUpset)
library(FUNGuildR)
library(reltools)
library(Biostrings)
library(RColorBrewer)

#### loading the data, subsetting for different questions ####
ITS.rel.phy <- readRDS("output/analysis/clean.ITS.rel.phy.rds")
ITS.abund.rel.phy <- readRDS("output/analysis/clean.ITS.abund.rel.phy.rds")

jaccard <- read.csv("data/vivi_vs_trfu_jaccard_values.csv")
sorensen <- read.csv("data/vivi_vs_trfu_sorensen_values.csv")
dist <- read.csv("data/poll.dist.T.csv")

ITS.rel.tax <- tax_table(ITS.rel.phy) %>% data.frame
ITS.rel.tax %>% write.csv("output/analysis/ITS.rel.ASVtaxonomy.csv")
refseq(ITS.rel.phy) %>% writeXStringSet(., file="output/analysis/ITS.rel.seqs.fasta", width=1000)
ITS.rel.tax2 <- read.csv("output/analysis/ITS.rel.ASV.taxonomy.corrected.csv")
ITS.rel.tax3 <- funguild_assign(ITS.rel.tax2)
ITS.rel.tax3$confidenceRanking <- factor(ITS.rel.tax3$confidenceRanking) # levels: Highly Probable, Probable, Possible
ITS.rel.tax3$trophicMode <- factor(ITS.rel.tax3$trophicMode)
ITS.rel.tax3$guild <- factor(ITS.rel.tax3$guild)
ITS.rel.tax3$growthForm <- factor(ITS.rel.tax3$growthForm)
ITS.rel.tax3 <- column_to_rownames(ITS.rel.tax3, var = "X")
tax_table(ITS.rel.phy) <- as.matrix(ITS.rel.tax3)

# subsetting by plant species
L.ITS.rel.phy <- subset_samples(ITS.rel.phy, Plant_species == "Goldfield")
L.ITS.rel.phy <- prune_taxa(taxa_sums(L.ITS.rel.phy)>0, L.ITS.rel.phy)
L.ITS.abund.rel.phy <- subset_samples(ITS.abund.rel.phy, Plant_species == "Goldfield")
L.ITS.abund.rel.phy <- prune_taxa(taxa_sums(L.ITS.abund.rel.phy)>0, L.ITS.abund.rel.phy)

T.ITS.rel.phy <- subset_samples(ITS.rel.phy, Plant_species == "Clover")
T.ITS.rel.phy <- prune_taxa(taxa_sums(T.ITS.rel.phy)>0, T.ITS.rel.phy)
T.ITS.abund.rel.phy <- subset_samples(ITS.abund.rel.phy, Plant_species == "Clover")
T.ITS.abund.rel.phy <- prune_taxa(taxa_sums(T.ITS.abund.rel.phy)>0, T.ITS.abund.rel.phy)

V.ITS.rel.phy <- subset_samples(ITS.rel.phy, Plant_species == "Vetch")
V.ITS.rel.phy <- prune_taxa(taxa_sums(V.ITS.rel.phy)>0, V.ITS.rel.phy)
V.ITS.abund.rel.phy <- subset_samples(ITS.abund.rel.phy, Plant_species == "Vetch")
V.ITS.abund.rel.phy <- prune_taxa(taxa_sums(V.ITS.abund.rel.phy)>0, V.ITS.abund.rel.phy)


#### summarizing the top 20 most abundant genera ####
brewerPlus <- distinct_palette(pal = "brewerPlus")

ITS.abund.merge <- ITS.abund.rel.phy %>%
  aggregate_top_taxa2(., 20, "genus") %>%
  plot_bar(., x="Sample", fill = "genus") +
  geom_bar(aes(fill=genus), stat="identity", position="stack", color="gray33") +
  facet_wrap(~Plant_species, scales = "free_x", nrow = 1) +
  theme(legend.position = "right", legend.key.size = unit(.5, 'cm'),
        legend.text = element_text(size=10), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  guides(fill=guide_legend(ncol = 1)) +
  ylab("Relative abundance (%)") +
  labs(title = "All samples")+
  scale_fill_manual(values = brewerPlus)
ITS.abund.merge

ITS.abund.merge.site <- ITS.abund.rel.phy %>%
  aggregate_top_taxa2(., 20, "genus") %>%
  plot_bar(., x="Sample", fill = "genus") +
  geom_bar(aes(fill=genus), stat="identity", position="stack", color="gray33") +
  facet_wrap(~Plant_species*Site, scales = "free_x", nrow = 3) +
  theme(legend.position = "right", legend.key.size = unit(.5, 'cm'),
        legend.text = element_text(size=10), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+  
  guides(fill=guide_legend(ncol = 1)) +
  ylab("Relative abundance (%)") +
  labs(title = "All samples")+
  scale_fill_manual(values = brewerPlus)
ITS.abund.merge.site

###without Unknown genera
ITS.abund.merge.site <- ITS.abund.rel.phy %>%
  aggregate_top_taxa2(20, "genus") %>%
  subset_taxa(!genus %in% c("Unknown")) %>%
  plot_bar(x = "Sample", fill = "genus") +
  geom_bar(aes(fill = genus), stat = "identity", position = "stack", color = "gray33") +
  facet_wrap(~Plant_species * Site, scales = "free_x", nrow = 3) +
  theme(
    legend.position = "right",
    legend.key.size = unit(.5, 'cm'),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  ylab("Relative abundance (%)") +
  labs(title = "All samples") +
  scale_fill_manual(values = brewerPlus)

ITS.abund.merge.site


#### Differential abundance of all taxa by plant species ####
output <- ancombc2(data = ITS.rel.phy, assay_name = "counts", tax_level = "genus",
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

# primary analysis 
tab_zero = output$zero_ind
print(head(tab_zero))
res_prim = output$res # Clover is the reference group

df_plants = res_prim %>%
  dplyr::select(taxon, contains("Plant_Species")) 
df_fig_plants1 = df_plants %>%
  dplyr::filter(diff_Plant_speciesGoldfield == T | 
                  diff_Plant_speciesVetch == T) %>%
  dplyr::mutate(lfc1 = ifelse(diff_Plant_speciesGoldfield == T, 
                              round(lfc_Plant_speciesGoldfield, 2), 0),
                lfc2 = ifelse(diff_Plant_speciesVetch == T, 
                              round(lfc_Plant_speciesVetch, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc2, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_plants2 = df_plants %>%
  dplyr::filter(diff_Plant_speciesGoldfield == T | 
                  diff_Plant_speciesVetch == T) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_Plant_speciesGoldfield == T & 
                                diff_Plant_speciesGoldfield == T, 
                              "aquamarine3", "black"),
                lfc2 = ifelse(passed_ss_Plant_speciesVetch == T &
                                diff_Plant_speciesVetch == T, 
                              "aquamarine3", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc2, 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_plants = df_fig_plants1 %>%
  dplyr::left_join(df_fig_plants2, by = c("taxon", "group"))

df_fig_plants$group = recode(df_fig_plants$group, 
                          `lfc1` = "Goldfield - Clover",
                          `lfc2` = "Vetch - Clover")
df_fig_plants$group = factor(df_fig_plants$group, 
                          levels = c("Vetch - Clover",
                                     "Goldfield - Clover"))

lo = floor(min(df_fig_plants$value))
up = ceiling(max(df_fig_plants$value))
mid = (lo + up)/2

df_fig_plants3 <- df_fig_plants %>%
  filter(color == "aquamarine3")

fig_plants <- df_fig_plants3 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to Clover samples") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_plants

### interpret results
global_df <- output$res_global
diff <- global_df %>% filter(passed_ss == TRUE) %>% filter(diff_abn == TRUE)
diff_clover_goldfield <- res_prim %>% filter(passed_ss_Plant_speciesGoldfield == TRUE) %>% filter(diff_Plant_speciesGoldfield == TRUE)
diff_clover_vetch <- res_prim %>% filter(passed_ss_Plant_speciesVetch == TRUE) %>% filter(diff_Plant_speciesVetch == TRUE)

#### Fungal ASV Richness #######
asv.div <- otu_table(ITS.rel.phy)
meta.div <- sample_data(ITS.rel.phy)
asv.rich <- vegan::specnumber(asv.div)
asv.rich.df <- as.data.frame(asv.rich)
asv.rich.df <- rownames_to_column(asv.rich.df, "sample")

asv.rich.df$Plant_species <- meta.div$Plant_species
asv.rich.df$Site <- meta.div$Site
asv.rich.df$Site_invasion_level <- meta.div$Site_invasion_level
asv.rich.df$Soil_type <- meta.div$Soil_type
asv.rich.df$dist_from_vetch_m <- meta.div$dist_from_vetch_m
asv.rich.df$NODFc <- meta.div$NODFc
asv.rich.df$Mean_Seed_Mass <- meta.div$Mean_Seed_Mass
asv.rich.df$Ratio_250_Log <- meta.div$Ratio_250_Log
asv.rich.df$Ind_Nestedness_Contribution <- meta.div$Ind_Nestedness_Contribution
asv.rich.df$Pollinator_Richness <- meta.div$Pollinator_Richness
asv.rich.df$Pollinator_Abundance <- meta.div$Pollinator_Abundance
asv.rich.df$Honeybee_Abundance <- meta.div$Honeybee_Abundance
asv.rich.df$Mean_Floral_Abundance <- meta.div$Mean_Floral_Abundance
asv.rich.df$Site_Plant_Richness <- meta.div$Site_Plant_Richness
asv.rich.df$Bombus_Abundance <- meta.div$Bombus_Abundance
asv.rich.df$Site_invasion_level <- factor(asv.rich.df$Site_invasion_level, 
                                          levels = c("Low", "Medium", "High"))
asv.rich.df$Plant_species <- factor(asv.rich.df$Plant_species,
                                    levels = c("Clover","Vetch","Goldfield"))
asv.rich.df$Ind_Nestedness_Contribution <- meta.div$Ind_Nestedness_Contribution

asv.rich.df$Site_lat <- meta.div$Site_lat
asv.rich.df$Site_long <- meta.div$Site_long

print(jaccard)
## add overlap in pollinators between vetch and clover
asv.rich.df <- asv.rich.df %>%
  mutate(Jaccard = case_when(
    Site == "Aikawa" ~ 0.28125000,
    Site == "Bertha" ~ 0.5,
    Site == "Goatgrass" ~ 0.15789474,
    Site == "Long" ~ 0.37500000,
    Site == "Quarry" ~ 0.3,
    Site == "Rock" ~ 0.53846154,
    TRUE ~ NA_real_  
  ))

print(head(sorensen))

asv.rich.df <- asv.rich.df %>%
  mutate(Sorensen = case_when(
    Site == "Aikawa" ~ 0.4390244,
    Site == "Bertha" ~ 0.6666667,
    Site == "Goatgrass" ~  0.2727273,
    Site == "Long" ~ 0.5454545,
    Site == "Quarry" ~ 0.4615385,
    Site == "Rock" ~ 0.7000000,
    TRUE ~ NA_real_  
  ))

asv.rich.df <- asv.rich.df %>%
  mutate(H2 = case_when(
    Site == "Aikawa" ~ 0.484937819117919,
    Site == "Bertha" ~ 0.45056284419381,
    Site == "Goatgrass" ~  0.498196535169117,
    Site == "Long" ~ 0.560694357715821,
    Site == "Quarry" ~ 0.719289128137327,
    Site == "Rock" ~ 0.468795896167944,
    TRUE ~ NA_real_  
  ))

asv.rich.df <- asv.rich.df %>%
  mutate(H_Shannon = case_when(
    Site == "Aikawa" ~ 2.96116274597892,
    Site == "Bertha" ~ 3.94049915379335,
    Site == "Goatgrass" ~  4.04506869844542,
    Site == "Long" ~ 3.07240408221974,
    Site == "Quarry" ~ 2.33866703525998,
    Site == "Rock" ~ 2.31306960029476,
    TRUE ~ NA_real_  
  ))

asv.rich.df <- asv.rich.df %>%
  mutate(d = case_when(
    Site == "Aikawa" & Plant_species == "Goldfield" ~ 0.644254208863882,
    Site == "Bertha" & Plant_species == "Goldfield" ~ 0.566120823676938,
    Site == "Goatgrass" & Plant_species == "Goldfield"  ~  0.561909612656833,
    Site == "Long" & Plant_species == "Goldfield"  ~ 0.816070060473425,
    Site == "Quarry" & Plant_species == "Goldfield"  ~ 0.62006693368836,
    Site == "Rock" & Plant_species == "Goldfield"  ~ NA,
    Site == "Aikawa" & Plant_species == "Vetch" ~ 0.253629434177368,
    Site == "Bertha" & Plant_species == "Vetch" ~ 0.57018425622141,
    Site == "Goatgrass" & Plant_species == "Vetch"  ~  0.789602352243406,
    Site == "Long" & Plant_species == "Vetch"  ~ 0.540154359006613,
    Site == "Quarry" & Plant_species == "Vetch"  ~ 0.343399065072067,
    Site == "Rock" & Plant_species == "Vetch"  ~ 0.151165000404385,
    Site == "Aikawa" & Plant_species == "Clover" ~ 0.423887959533828,
    Site == "Bertha" & Plant_species == "Clover" ~ 0.470262888886576,
    Site == "Goatgrass" & Plant_species == "Clover"  ~  0.326163337231031,
    Site == "Long" & Plant_species == "Clover"  ~ 0.205938156139874,
    Site == "Quarry" & Plant_species == "Clover"  ~ 0.791817000205232,
    Site == "Rock" & Plant_species == "Clover"  ~ 0.389155979130123,
    TRUE ~ NA_real_  
  ))
asv.rich.T <- asv.rich.df %>%
  filter(Plant_species == "Clover") 

asv.rich.V <- asv.rich.df %>%
  filter(Plant_species == "Vetch") 
asv.rich.L <- asv.rich.df %>%
  filter(Plant_species == "Goldfield") 



ggplot(asv.rich.df, aes(x=Plant_species, y=asv.rich)) + 
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  ylab("Richness") +
  facet_wrap(~Site_invasion_level) + theme_classic() +
  # scale_fill_viridis(discrete = T, option = "C") +
  theme(axis.text = element_text(size = rel(1.5)), axis.title = element_text(size = rel(1.75)), 
        legend.text = element_text(size = rel(2)))

asv.rich.T$dist_from_vetch_m <- as.factor(asv.rich.T$dist_from_vetch_m)
asv.rich.T$Site_Plant_Richness <- as.factor(asv.rich.T$Site_Plant_Richness)
asv.rich.T$Pollinator_Richness <- as.factor(asv.rich.T$Pollinator_Richness)
asv.rich.T$NODFc <- as.factor(asv.rich.T$NODFc)


kruskal.test(asv.rich ~ Plant_species, data = asv.rich.df)


#### Fungal Shannon Alpha diversity ####
asv.div <- otu_table(ITS.rel.phy)
meta.div <- sample_data(ITS.rel.phy)
asv.shannon <- vegan::diversity(asv.div, index = "shannon", MARGIN = 1, base = exp(1))
asv.shannon.df <- as.data.frame(asv.shannon)
asv.shannon.df <- rownames_to_column(asv.shannon.df, "sample")

# adding metadata variables to the diversity dataframe
asv.shannon.df$Plant_species <- meta.div$Plant_species
asv.shannon.df$Site <- meta.div$Site
asv.shannon.df$Site_invasion_level <- meta.div$Site_invasion_level
asv.shannon.df$Soil_type <- meta.div$Soil_type
asv.shannon.df$dist_from_vetch_m <- meta.div$dist_from_vetch_m
asv.shannon.df$NODFc <- meta.div$NODFc
asv.shannon.df$Mean_Seed_Mass <- meta.div$Mean_Seed_Mass
asv.shannon.df$Ratio_250_Log <- meta.div$Ratio_250_Log
asv.shannon.df$Ind_Nestedness_Contribution <- meta.div$Ind_Nestedness_Contribution
asv.shannon.df$Pollinator_Richness <- meta.div$Pollinator_Richness
asv.shannon.df$Pollinator_Abundance <- meta.div$Pollinator_Abundance
asv.shannon.df$Honeybee_Abundance <- meta.div$Honeybee_Abundance
asv.shannon.df$Mean_Floral_Abundance <- meta.div$Mean_Floral_Abundance
asv.shannon.df$Site_Plant_Richness <- meta.div$Site_Plant_Richness
asv.shannon.df$Bombus_Abundance <- meta.div$Bombus_Abundance
asv.shannon.df$Plant_species <- factor(asv.shannon.df$Plant_species,
                                       levels = c("Clover","Vetch","Goldfield"))
asv.shannon.df$Site_invasion_level <- factor(asv.shannon.df$Site_invasion_level,
                                             levels = c("Low","Medium","High"))

asv.shannon.df$Site_long <- meta.div$Site_long
asv.shannon.df$Site_lat <- meta.div$Site_lat

asv.shannon.df <- asv.shannon.df %>%
  mutate(Jaccard = case_when(
    Site == "Aikawa" ~ 0.28125000,
    Site == "Bertha" ~ 0.5,
    Site == "Goatgrass" ~ 0.15789474,
    Site == "Long" ~ 0.37500000,
    Site == "Quarry" ~ 0.3,
    Site == "Rock" ~ 0.53846154,
    TRUE ~ NA_real_  
  ))

print(head(sorensen))

asv.shannon.df <- asv.shannon.df %>%
  mutate(Sorensen = case_when(
    Site == "Aikawa" ~ 0.4390244,
    Site == "Bertha" ~ 0.6666667,
    Site == "Goatgrass" ~  0.2727273,
    Site == "Long" ~ 0.5454545,
    Site == "Quarry" ~ 0.4615385,
    Site == "Rock" ~ 0.7000000,
    TRUE ~ NA_real_  
  ))


print(head(asv.shannon.df))

asv.shannon.df <- asv.shannon.df %>%
  mutate(Distance = case_when(
    dist_from_vetch_m %in% c(0, 1) ~ "0-1",
    dist_from_vetch_m %in% c(5, 10) ~ "5-10",
    dist_from_vetch_m == 50 ~ "50",
    TRUE ~ NA_character_  # Preserve NA values
  ))

asv.shannon.df <- asv.shannon.df %>%
  mutate(H2 = case_when(
    Site == "Aikawa" ~ 0.484937819117919,
    Site == "Bertha" ~ 0.45056284419381,
    Site == "Goatgrass" ~  0.498196535169117,
    Site == "Long" ~ 0.560694357715821,
    Site == "Quarry" ~ 0.719289128137327,
    Site == "Rock" ~ 0.468795896167944,
    TRUE ~ NA_real_  
  ))

asv.shannon.df <- asv.shannon.df %>%
  mutate(H_Shannon = case_when(
    Site == "Aikawa" ~ 2.96116274597892,
    Site == "Bertha" ~ 3.94049915379335,
    Site == "Goatgrass" ~  4.04506869844542,
    Site == "Long" ~ 3.07240408221974,
    Site == "Quarry" ~ 2.33866703525998,
    Site == "Rock" ~ 2.31306960029476,
    TRUE ~ NA_real_  
  ))

asv.shannon.df <- asv.shannon.df %>%
  mutate(d = case_when(
    Site == "Aikawa" & Plant_species == "Goldfield" ~ 0.644254208863882,
    Site == "Bertha" & Plant_species == "Goldfield" ~ 0.566120823676938,
    Site == "Goatgrass" & Plant_species == "Goldfield"  ~  0.561909612656833,
    Site == "Long" & Plant_species == "Goldfield"  ~ 0.816070060473425,
    Site == "Quarry" & Plant_species == "Goldfield"  ~ 0.62006693368836,
    Site == "Rock" & Plant_species == "Goldfield"  ~ NA,
    Site == "Aikawa" & Plant_species == "Vetch" ~ 0.253629434177368,
    Site == "Bertha" & Plant_species == "Vetch" ~ 0.57018425622141,
    Site == "Goatgrass" & Plant_species == "Vetch"  ~  0.789602352243406,
    Site == "Long" & Plant_species == "Vetch"  ~ 0.540154359006613,
    Site == "Quarry" & Plant_species == "Vetch"  ~ 0.343399065072067,
    Site == "Rock" & Plant_species == "Vetch"  ~ 0.151165000404385,
    Site == "Aikawa" & Plant_species == "Clover" ~ 0.423887959533828,
    Site == "Bertha" & Plant_species == "Clover" ~ 0.470262888886576,
    Site == "Goatgrass" & Plant_species == "Clover"  ~  0.326163337231031,
    Site == "Long" & Plant_species == "Clover"  ~ 0.205938156139874,
    Site == "Quarry" & Plant_species == "Clover"  ~ 0.791817000205232,
    Site == "Rock" & Plant_species == "Clover"  ~ 0.389155979130123,
    TRUE ~ NA_real_  
  ))

asv.shannon.T <- asv.shannon.df %>%
  filter(Plant_species == "Clover") 
asv.shannon.V <- asv.shannon.df %>%
  filter(Plant_species == "Vetch") 
asv.shannon.L <- asv.shannon.df %>%
  filter(Plant_species == "Goldfield") 


shannon.aov <- aov(asv.shannon~Plant_species, data = asv.shannon.df)
summary(shannon.aov)
TukeyHSD(shannon.aov)

#### Comparing community composition with ordination and perMANOVA ####
# making the NMDS plot - all samples
set.seed(12)
fun.bray <- ITS.rel.phy %>% phyloseq::distance("bray") %>% sqrt()
fun.nmds <- metaMDS(fun.bray, trymax=200, parallel=10, k=2) 
fun.nmds.dat <- scores(fun.nmds, display = "site") %>% data.frame %>% 
  rownames_to_column("sampID") %>%
  full_join(sample_data(ITS.rel.phy) %>% data.frame %>% rownames_to_column("sampID"))


fun.nmds.dat$Plant_species <- factor(fun.nmds.dat$Plant_species,
                                      levels = c("Vetch", "Clover", "Goldfield"))  # Specify the order 

nmds.plot <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(shape=Site, color=Plant_species),size=3)+
  stat_ellipse(aes(color=Plant_species))+
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(fill="Plant Species", color = "Plant Species", shape="Site")
nmds.plot

fun.nmds.dat$Site_invasion_level <- factor(fun.nmds.dat$Site_invasion_level,
                                           levels=c("Low","Medium","High"))

nmds.plot.2 <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(shape=Site_invasion_level, color=Plant_species),size=3)+
  stat_ellipse(aes(color=Plant_species))+
  facet_grid(~Site_invasion_level)+
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(fill="Plant Species", color = "Plant Species", shape="Site Invasion Level")
nmds.plot.2 # Gillian: I facetted the plot by invasion level to make the comparisons more clear

fun.nmds.dat$dist_from_vetch_m <- as.factor(fun.nmds.dat$dist_from_vetch_m)
nmds.plot.3 <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(color=dist_from_vetch_m, shape=Site),size=3)+
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(fill="Site", shape= "Site", color="Distance from vetch") + facet_wrap(~Site)
nmds.plot.3

fun.nmds.dat <- fun.nmds.dat %>%
  mutate(dist_from_vetch_m = case_when(
    Plant_species == "Vetch" ~ "0",  
    Plant_species == "Goldfield" ~ "50",  
    TRUE ~ dist_from_vetch_m
  ))

fun.nmds.dat <- fun.nmds.dat %>%
  mutate(Distance = case_when(
    dist_from_vetch_m == "0" ~ "0-1",  
    dist_from_vetch_m == "1" ~ "0-1", 
    dist_from_vetch_m == "50" ~ "50",  
    TRUE ~ "5-10"
  ))

nmds.plot <- ggplot(fun.nmds.dat, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(color = Plant_species, shape = Distance), size = 3) + 
  scale_color_viridis(discrete = TRUE, option = "C") +  # Color by species
  ggthemes::theme_few() +
  labs(color = "Plant Species", shape = "Distance from Vetch") +
  facet_wrap(~Site)  
nmds.plot


T.bray <- T.ITS.rel.phy %>% phyloseq::distance("bray") %>% sqrt()
T.nmds <- metaMDS(T.bray, trymax=200, parallel=10, k=2) 
T.nmds.dat <- scores(T.nmds, display = "site") %>% data.frame %>% 
  rownames_to_column("sampID") %>%
  full_join(sample_data(T.ITS.rel.phy) %>% data.frame %>% rownames_to_column("sampID"))

T.nmds.dat$dist_from_vetch_m <- as.factor(T.nmds.dat$dist_from_vetch_m)
T.nmds.dat$Site_invasion_level <- factor(T.nmds.dat$Site_invasion_level,
                                         levels = c("Low","Medium","High"))

T.nmds.dat <- T.nmds.dat %>%
  mutate(Distance = case_when(
    dist_from_vetch_m == "0" ~ "0-1",  
    dist_from_vetch_m == "1" ~ "0-1", 
    dist_from_vetch_m == "50" ~ "50",  
    TRUE ~ "5-10"
  ))

distance.nmds.plot <- ggplot(T.nmds.dat, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(color = Distance), size = 3) + 
  scale_color_viridis(discrete = TRUE, option = "D") +  # Color by species
  stat_ellipse(aes(color = Distance))+
  ggthemes::theme_few() +
  labs(color = "Distance\nfrom Vetch (m)") +
  facet_wrap(~Site_invasion_level)+
  theme(axis.text = element_text(size=15),
        axis.title=element_text(size=18),
        legend.text = element_text(size=15),
        legend.title = element_text(size=18),
        strip.text = element_text(size=20))
distance.nmds.plot

nmds.plot <- ggplot(fun.nmds.dat, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(color = Plant_species, shape = Site), size = 3) + 
  scale_color_viridis(discrete = TRUE, option = "C") +  # Color by species
  ggthemes::theme_few() +
  labs(color = "Plant Species", shape = "Site") +
  facet_wrap(~Distance)  
nmds.plot


#version with site
fun.nmds.dat$dist_from_vetch_m <- as.factor(fun.nmds.dat$dist_from_vetch_m)
nmds.plot.3 <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(shape=dist_from_vetch_m, color=Site),size=3)+
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  labs(fill="Site", color = "Site", shape="Distance from Vetch")
nmds.plot.3

fun.nmds.dat$Honeybee_Abundance <- as.factor(fun.nmds.dat$Honeybee_Abundance)
nmds.plot.4 <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2))+ 
  geom_point(aes(shape=Site, color=Honeybee_Abundance),size=3)+
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  ggthemes::theme_few() +
  facet_wrap(~Plant_species) +
  labs(fill="Honeybee_Abundance", color = "Honeybee_Abundance", shape="Site")
nmds.plot.4

# perMANOVA - all samples
fun.dist <- ITS.rel.phy %>% phyloseq::distance("bray")  
perm <- adonis2(fun.dist~sample_data(ITS.rel.phy)$Plant_species) 

   

fun.dist <- T.ITS.rel.phy %>% phyloseq::distance("bray")  
perm <- adonis2(fun.dist~sample_data(T.ITS.rel.phy)$Site_invasion_level) 

fun.dist <- T.ITS.rel.phy %>% phyloseq::distance("bray")  
perm <- adonis2(fun.dist~sample_data(T.ITS.rel.phy)$Ratio_250_Log) 


# Extract sample data and remove NA values
sample_data_filtered <- sample_data(T.ITS.rel.phy) %>%
  data.frame() %>%
  filter(!is.na(dist_from_vetch_m))

# Prune the phyloseq object to keep only samples with non-NA Distance
T.ITS.rel.phy.filtered <- prune_samples(
  rownames(sample_data(T.ITS.rel.phy)) %in% rownames(sample_data_filtered),
  T.ITS.rel.phy
)

# Recalculate the distance matrix on the filtered phyloseq object
fun.dist <- phyloseq::distance(T.ITS.rel.phy.filtered, "bray")



T.ITS.rel.meta <- sample_data(T.ITS.rel.phy) %>% data.frame
T.ITS.rel.meta2 <- T.ITS.rel.meta %>%
  mutate(Distance = case_when(
    dist_from_vetch_m == "0" ~ "0-1",  
    dist_from_vetch_m == "1" ~ "0-1", 
    dist_from_vetch_m == "50" ~ "50",  
    TRUE ~ "5-10"
  ))
T.ITS.rel.meta2$Distance <- factor(T.ITS.rel.meta2$Distance,
                                   levels = c("0-1",'5-10','50'))

sample_data(T.ITS.rel.phy) <- T.ITS.rel.meta2
fun.dist <- T.ITS.rel.phy %>% phyloseq::distance("bray")  
perm <- adonis2(fun.dist~sample_data(T.ITS.rel.phy)$Distance*
                  sample_data(T.ITS.rel.phy)$Site_invasion_level)
perm


T.low.ITS.rel.phy <- subset_samples(T.ITS.rel.phy, Site_invasion_level == "Low")
T.low.ITS.rel.phy<-prune_taxa(taxa_sums(T.low.ITS.rel.phy)>0, T.low.ITS.rel.phy)
T.mid.ITS.rel.phy <- subset_samples(T.ITS.rel.phy, Site_invasion_level == "Medium")
T.mid.ITS.rel.phy<-prune_taxa(taxa_sums(T.mid.ITS.rel.phy)>0, T.mid.ITS.rel.phy)
T.high.ITS.rel.phy <- subset_samples(T.ITS.rel.phy, Site_invasion_level == "High")
T.high.ITS.rel.phy<-prune_taxa(taxa_sums(T.high.ITS.rel.phy)>0, T.high.ITS.rel.phy)

low.dist <- T.low.ITS.rel.phy %>% phyloseq::distance("bray")  
low.perm <- adonis2(low.dist~sample_data(T.low.ITS.rel.phy)$Distance)
low.perm


mid.dist <- T.mid.ITS.rel.phy %>% phyloseq::distance("bray")  
mid.perm <- adonis2(mid.dist~sample_data(T.mid.ITS.rel.phy)$Distance)
mid.perm
 

high.dist <- T.high.ITS.rel.phy %>% phyloseq::distance("bray")  
high.perm <- adonis2(high.dist~sample_data(T.high.ITS.rel.phy)$Distance)
high.perm


#### Vector analysis to check correlations with site variables and ASVs ####
# analyzing metadata
fun.var <- as.data.frame(scores(fun.nmds, display = 'sites'))
fun.meta <- as.matrix(sample_data(ITS.rel.phy)) %>% as.data.frame()
fun.var <- fun.meta %>% 
  select(Site, Plant_species, Site_invasion_level, NODFc, Pollinator_Abundance, Pollinator_Richness, Bombus_Abundance, Honeybee_Abundance, Mean_Floral_Abundance, Ind_Nestedness_Contribution, Site_Plant_Richness) %>% 
  cbind(fun.var)
fun.fit <- envfit(fun.nmds, fun.var, na.rm = TRUE)
fun.filt <- data.frame(r = fun.fit$vectors$r,pvals = fun.fit$vectors$pvals) %>%
  rownames_to_column(var = "var") 

fun.var.scrs <- as.data.frame(scores(fun.fit, display = "vectors"))
fun.var.scrs <- cbind(fun.var.scrs, variables = rownames(fun.var.scrs)) 
#fun.var.scrs %>% filter(variables %in% fun.filt$var)

vector.plot <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2, color=Plant_species)) +
  geom_point() +
  coord_fixed() +
  scale_color_viridis(discrete = T, option = "C") +
  geom_segment(data = fun.var.scrs, linewidth = 2,
               aes(x = 0, xend = NMDS1*.8, y = 0, yend = NMDS2*.8), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  geom_text(data = fun.var.scrs, aes(x = NMDS1, y = NMDS2, label = variables), inherit.aes = FALSE, size = 5, hjust = 1, vjust = 0)
vector.plot

# analyzing if certain ASVs drive community composition (top 25 ASVs from above)
top.asv <- names(sort(taxa_sums(ITS.rel.phy), TRUE)[1:25])
top.asv.phy <- prune_taxa(top.asv, ITS.rel.phy)
top.asv.meta <- otu_table(top.asv.phy) %>% data.frame
fun.asv.var <- as.data.frame(scores(fun.nmds, display = 'sites'))
fun.asv.var <- top.asv.meta %>% cbind(fun.asv.var) 
fun.asv.fit <- envfit(fun.nmds, fun.asv.var, na.rm = TRUE)
fun.asv.filt <- data.frame(r = fun.asv.fit$vectors$r,pvals = fun.asv.fit$vectors$pvals) %>%
  rownames_to_column(var = "var") %>% 
  filter(r >= 0.1, pvals <= 0.002) %>%
  filter(var!="NMDS1") %>%
  filter(var!="NMDS2")

fun.asv.scrs <- as.data.frame(scores(fun.asv.fit, display = "vectors"))
fun.asv.scrs <- cbind(fun.asv.scrs, variables = rownames(fun.asv.scrs)) 
fun.asv.scrs <- fun.asv.scrs %>%
  filter(variables %in% fun.asv.filt$var)
fun.asv.scrs$ASV.taxa <- c("Vishniacozyma_victoriae", "Filobasidium_stepposum", "Papiliotrema_sp",
                           "Cladosporium_cladosporioides","Cladosporium_sp1", "Cladosporium_sp2",
                           "Podosphaera_sp", "Alternaria_eureka",
                            "Alternaria_sp", "Unk_Nectriaceae", "Sporobolomyces_sp",
                           "Unk")

asv.vector.plot <- ggplot(fun.nmds.dat,aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(color=Plant_species),size=4) +
  coord_fixed() +
  geom_segment(data = fun.asv.scrs, size = 1,
               aes(x = 0, xend = NMDS1*.8, y = 0, yend = NMDS2*.8), inherit.aes = FALSE,
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  geom_text(data = fun.asv.scrs, aes(x = NMDS1, y = NMDS2, label = ASV.taxa), inherit.aes = FALSE, size = 2, hjust = 0, vjust = 0) +
  theme(text = element_text(size=15)) +
  ggthemes::theme_few()+
  scale_color_viridis(discrete = T, option = "C")
asv.vector.plot

###### Beta dispersion ########
table(asv.shannon.df$Plant_species)

plant_species <- factor(c(rep(1,33), rep(2,21), rep(3,98)), 
                        labels = c("Vetch", "Goldfield", "Clover"))
fung.dist <- vegdist(otu_table(ITS.rel.phy))
fung.disper <- betadisper(fung.dist, plant_species)
permutest(fung.disper) 

TukeyHSD(fung.disper)

plot(fung.disper) 

asv.shannon.T %>% filter(Site == "Rock") #to get replicates per site
table(asv.shannon.T$Site)

sites <- factor(c(rep(1,17), rep(2,11), rep(3,21), rep(4,11), rep(5,25), rep(6,13)),
                labels = c("Aikawa", "Quarry", "Bertha", "Goatgrass", "Long", "Rock"))
t.fung.dist <- vegdist(otu_table(T.ITS.rel.phy))
site.disper <- betadisper(t.fung.dist, sites)
permutest(site.disper)
TukeyHSD(site.disper)
plot(site.disper)


table(asv.shannon.T$Site_invasion_level)

sites <- factor(c(rep(1,32), rep(2,28), rep(3,38)),
                labels = c("Low", "Medium", "High"))
t.fung.dist <- vegdist(otu_table(T.ITS.rel.phy))
invas.disper <- betadisper(t.fung.dist, sites)
permutest(invas.disper)
TukeyHSD(invas.disper)
plot(invas.disper)

############# Diversity Begets Diversity analysis ##############
## run the shannon and richness parts first 

### datasets: asv.rich.L, asv.rich.T, asv.rich.V,asv.shannon.L, asv.shannon.T, asv.shannon.V
### asv.rich.df T, L, an V together -- Plant_species, distinguishes them, asv.shannon.df shows them all together 
## response variabels asv.rich for asv.rich datasets and asv.shannon for asv.shannon datasets
## predictor variables: Site_Plant_Richness, Mean_Floral_Abundance, NODFc, Pollinator Richness, Pollinator_Abundance, Honeybee_Abundance, Bombus_Abundance, Ind_Nestedness_Contribution, H2', d'


library(purrr)
library(broom)

responses <- list(
  rich    = list(var = "asv.rich",    prefix = "asv.rich"),
  shannon = list(var = "asv.shannon", prefix = "asv.shannon")
)

predictors <- c(
  "Site_Plant_Richness",
  "Mean_Floral_Abundance",
  "NODFc",
  "Pollinator_Richness",
  "Pollinator_Abundance",
  "Honeybee_Abundance",
  "Bombus_Abundance",
  "Ind_Nestedness_Contribution",
  "H2",
  "H_Shannon",
  "d"
)

datasets <- list(
  L = list(rich = asv.rich.L,    shannon = asv.shannon.L),
  T = list(rich = asv.rich.T,    shannon = asv.shannon.T),
  V = list(rich = asv.rich.V,    shannon = asv.shannon.V)
)

## aggregate site means for fungal diversity across samples, then fit linear models 
fit_site_lm <- function(df, response, predictor) {
  
  site_means <- df %>%
    dplyr::group_by(Site, .data[[predictor]]) %>%
    dplyr::summarise(
      mean_response = mean(.data[[response]], na.rm = TRUE),
      .groups = "drop"
    )
  

  form <- reformulate(predictor, response = "mean_response")
  
  mod <- lm(form, data = site_means)
  
  list(
    data   = site_means,
    model  = mod,
    tidy   = broom::tidy(mod),
    glance = broom::glance(mod)
  )
}


results <- tidyr::expand_grid(
  Species   = names(datasets),
  Response  = names(responses),
  Predictor = predictors
) %>%
  dplyr::mutate(
    fit = purrr::pmap(
      list(Species, Response, Predictor),
      function(Species, Response, Predictor) {
        
        df <- datasets[[Species]][[Response]]
        
        fit_site_lm(
          df        = df,
          response  = responses[[Response]]$var,
          predictor = Predictor
        )
      }
    )
  )



check_assumptions <- function(mod) {
  par(mfrow = c(1, 2))
  plot(fitted(mod), resid(mod),
       xlab = "Fitted", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, lty = 2)
  
  qqnorm(resid(mod))
  qqline(resid(mod))
}

results$fit[[31]]$model %>% check_assumptions()


model_table <- results %>%
  mutate(
    estimate = map_dbl(fit, ~ .x$tidy$estimate[2]),
    std.error = map_dbl(fit, ~ .x$tidy$std.error[2]),
    tval = map_dbl(fit, ~ .x$tidy$statistic[2]),
    pval  = map_dbl(fit, ~ .x$tidy$p.value[2]),
    r2    = map_dbl(fit, ~ .x$glance$r.squared)
  ) %>%
  select(Species, Response, Predictor, estimate, std.error, tval, pval, r2)


library(flextable)
library(officer)

ft <- flextable(model_table) %>%
  colformat_double(j = c("estimate","std.error","tval","pval","r2"), digits = 3) %>%
  autofit() %>%
  set_header_labels(
    Species    = "Plant Species",
    Response   = "Response",
    Predictor  = "Predictor",
    estimate      = "Effect Size (Slope)",
    std.error  = "Std. Error",
    tval    = "t-value",
    pval       = "p-value",
    r2         = "R²"
  )

read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = "Fungal_Diversity_Model_Table.docx")



results <- results %>%
  mutate(
    slope  = map_dbl(fit, ~ {
      if (!is.null(.x)) .x$tidy$estimate[2] else NA_real_
    }),
    pval   = map_dbl(fit, ~ {
      if (!is.null(.x)) .x$tidy$p.value[2] else NA_real_
    }),
    signif = ifelse(pval < 0.05, "sig", "ns")
  )


site_means_all <- results %>%
  mutate(data = map(fit, "data")) %>%
  select(Species, Response, Predictor, signif, data) %>%
  unnest(data)


species_colors <- c(
  L = "#DAA520", 
  V = "#800080", 
  T = "#FF69B4"  
)

plot_model <- function(resp, pred) {
  
  df_plot <- site_means_all %>%
    filter(Response == resp, Predictor == pred)
  

  if(nrow(df_plot) == 0) return(NULL)
  
  ggplot(df_plot, aes(x = .data[[pred]], y = mean_response,
                      color = Species, linetype = signif)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_color_manual(values = species_colors) +
    scale_linetype_manual(values = c(sig = "solid", ns = "dashed")) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = rel(1.4)),
      axis.title = element_text(size = rel(1.6)),
      legend.text = element_text(size = rel(1.3))
    ) +
    labs(
      x = pred,
      y = ifelse(resp == "rich", "Mean ASV richness", "Mean ASV Shannon"),
      color = "Plant species",
      linetype = "Significance"
    )
}

plots <- expand_grid(
  Response  = names(responses),
  Predictor = predictors
) %>%
  mutate(plot = map2(Response, Predictor, ~ plot_model(.x, .y)))



walk2(plots$plot, paste0(plots$Response, "_", plots$Predictor, ".png"), 
      ~ if(!is.null(.x)) ggsave(filename = .y, plot = .x, width = 6, height = 5, dpi = 300))


ggplot(
  site_means_all %>%
    filter(Response == "rich",
           Predictor == "Site_Plant_Richness"),
  aes(
    x = Site_Plant_Richness,
    y = mean_response,
    color = Species
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  theme_classic() +
  theme(
    axis.text  = element_text(size = rel(1.4)),
    axis.title = element_text(size = rel(1.6)),
    legend.text = element_text(size = rel(1.3))
  ) +
  labs(
    x = "Plant species richness (site-level)",
    y = "Mean ASV richness",
    color = "Plant species"
  )





