# R code for the publication "More fungi than legs: the first
# fungal microbiome for a fungus-eating millipede (Colobognatha)"

# Code authored by Dr.Brian Lovett and Angie Macias
# Written with R v. 4.0.3; R Studio v. 1.4.1103

#### Set working directory ####
setwd("/Users/angiemacias/Desktop/brachyillumina081221/Rwork")

#### Load libraries ####
library(tidyverse)       # v. 1.3.1
library(vegan)           # v. 2.5-7
library(ggfortify)       # v. 0.4.14
library(rstatix)         # v. 0.7.0
library(pairwiseAdonis)  # v. 0.4
library(ape)             # v. 5.6-1
library(BSDA)            # v. 1.2.1
library(tidytext)        # v. 0.3.2
library(janitor)         # v. 2.2.0
library(effectsize)      # v. 0.8.9
library(ggrepel)         # v. 0.9.4
library(grid)            # v. 4.0.3
library(gridExtra)       # v. 2.3

#### Import data ####

# General imports
bst2 <- read.delim("brachy_illumina_v3_wTax_corr.txt", header=T) # BASTA results
bm2 <- read.delim("brachy_miseq_v3.biom.txt", comment.char = "#") # amptk biom file
tx <- read.delim("brachy_miseq_v3.otu_table.taxonomy.txt", comment.char = "#", header=F) # just OTU number and taxonomy
metadata <- read.delim("brachy_miseq_v3.mapping_file_trimmed.txt", header=T) # CF Metadata
# Functional analysis imports
FungalTraits <- read.csv("FungalTraits1.2v16.csv", header=T) # FungalTraits db
CultureBased <- read.delim("CultureBased_OTUtable.txt", header=T, check.names=F) # CB OTU table
CultureFree <- read.delim("CultureFree_OTUtable.txt", header=T, check.names=F) # CF OTU table (dup for convenience)
CultureBasedwFunGuild <- read.delim("CultureBased_FunGuild.txt", header=T, check.names=F) %>%
  filter(Guild != "-") # CB FunGuild results
CultureFreewFunGuild <- read.delim("CultureFree_FunGuild.txt", header=T, check.names=F) %>%
  filter(Guild != "-") # CF FunGuild results
CultureBasedMeta <- read_tsv("CultureBased_Metadata.txt") # CB metadata
CultureFreeMeta <- read_tsv("CultureFree_Metadata.txt") # CF metadata

#### Pre-analysis data tidying ####

# This block modifies the original amptk taxonomy file to a useable format.
# Note, this is NOT the OTU table, just OTUs and taxonomy.
# Bug fixes: notax error with OTU2233, replace empty internal taxonomy with
# incertae sedis while leaving trailing empty taxonomy blank.
colnames(tx)[1] <- "OTU.ID"
colnames(tx)[62] <- "tx"
tx2 <- tx %>%
  select(colnames(tx)[c(1, length(colnames(tx)))]) %>%
  separate(tx, into=c("toss", "taxonomy"), sep=";") %>%
  select(-toss) %>%
  mutate(taxonomy=ifelse(OTU.ID=="OTU2233", "k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Sordariales,f:Chaetomiaceae", taxonomy)) %>% #AMPTK'S FAULT NOT OURS
  separate(taxonomy, into=as.character(c(1:7)), sep=",") %>%
  pivot_longer(cols=-OTU.ID, values_to="Words", names_to="Level") %>%
  drop_na(Words) %>%
  separate(Words, into=c("Level", "Words"), sep=":") %>%
  mutate(Words=ifelse(Level=="No hit", "No hit", Words), Level=ifelse(Words=="No hit", "k", Level)) %>%
  pivot_wider(id_cols="OTU.ID", names_from="Level", values_from="Words") %>%
  mutate(`s`=ifelse(is.na(`s`), "drop", `s`),
         `g`=ifelse(is.na(`g`) & `s`=="drop", "drop", `g`),
         `f`=ifelse(is.na(`f`) & `g`=="drop", "drop", `f`),
         `o`=ifelse(is.na(`o`) & `f`=="drop", "drop", `o`),
         `c`=ifelse(is.na(`c`) & `o`=="drop", "drop", `c`),
         `p`=ifelse(is.na(`p`) & `c`=="drop", "drop", `p`),
         `k`=ifelse(is.na(`k`) & `p`=="drop", "drop", `k`)) %>%
  replace(is.na(.), "incertae sedis") %>%
  pivot_longer(cols=-OTU.ID, values_to="Words", names_to="Level") %>%
  filter(Words!="drop") %>%
  mutate(Label=paste(Level, Words, sep="__")) %>%
  mutate_at("Label", str_replace, "k__No hit", "No hit") %>%
  select(OTU.ID, Label) %>%
  group_by(OTU.ID) %>%
  summarize(taxonomy=toString(Label), depth=ifelse(Label=="No hit", 0, length(Label))) %>%
  distinct()

# This block modifies the BASTA output file to a usable format.
# Bug fix: replaces "unknown" (how BASTA handles empty internal taxonomy)
# with incertae sedis
bst3 <- bst2 %>%
  separate(taxonomy, into=c("k", "p", "c", "o", "f", "g", "s"), sep=";") %>%
  pivot_longer(cols=-OTU.ID, values_to="Words", names_to="Level") %>%
  filter(Words!="") %>%
  mutate_at("Words", str_replace_all, "unknown", "incertae sedis") %>%
  drop_na(Words) %>%
  mutate(Label=paste(Level, Words, sep="__")) %>%
  select(OTU.ID, Label) %>%
  mutate_at("Label", str_replace, "k__No Hit", "No hit") %>%
  group_by(OTU.ID) %>%
  summarize(taxonomy_bst=toString(Label), depth_bst=ifelse(Label=="No hit", 0, length(Label))) %>%
  distinct()

# This block combines the amptk taxonomy file and the BASTA output
# file, then prunes out duplicates (made by adding the BASTA data),
# keeping the row with the greatest taxonomic depth.
comb.tax <- left_join(tx2, bst3)
comb.tax2 <- comb.tax %>%
  mutate(depth_bst=ifelse(is.na(depth_bst), 0, depth_bst)) %>%
  mutate(better=depth_bst>depth) %>%
  mutate(taxonomy_orig=taxonomy, taxonomy=ifelse(better==T, taxonomy_bst, taxonomy)) %>%
  select(OTU.ID, taxonomy)%>%
  mutate_at("taxonomy", str_replace_all, ",", ";")

# This block modifies the BIOM file from amptk, pasting in the new
# taxonomy set constructed above. It also modifies some column names.
bm3 <- left_join(bm2[1:(length(colnames(bm2))-1)], comb.tax2)
colnames(bm3)[1]="OTU ID"
colnames(bm3)[grep("X", colnames(bm3))]=str_remove(colnames(bm3)[grep("X", colnames(bm3))], "X")

# Now we should determine what samples and/or OTUs need to be dropped.

# This file is just the OTU name & read counts, no tax.
OTUs <- bm3[1:(length(colnames(bm3))-1)]

# tx3 and 4 are the OTU name and taxonomy in different formats.
tx3 <- bm3[c(1,length(colnames(bm3)))] %>%
  separate(taxonomy, into=c("k", "p", "c", "o", "f", "g", "s"), sep=";") %>%
  pivot_longer(cols=-`OTU ID`, values_to="Words", names_to="Level") %>%
  drop_na(Words) %>%
  separate(Words, into=c("toss", "Taxons"), sep="__") %>%
  mutate(Taxons=ifelse(toss=="No hit", "No hit", Taxons)) %>%
  select(-toss)
tx4 <- tx3 %>%
  rowwise() %>%
  mutate(Label=paste(Level, Taxons, sep="_")) %>%
  group_by(`OTU ID`) %>%
  summarize(Taxonomy=toString(Label))

# This block reformats the OTU table, for examining OTUs found in
# the negative controls.
OTUs2 <- OTUs %>%
  pivot_longer(cols=-`OTU ID`, names_to="Sample", values_to="Reads")
Negs <- OTUs2 %>%
  filter(Sample=="Neg1" | Sample=="Neg2" | Sample=="Neg3") %>%
  filter(Reads>0) %>%
  pivot_wider(names_from="Sample", values_from="Reads")

# 14 OTUs were found in 1+ negative controls.

# Calculating stats to help us make decisions about which OTUs to keep
Negs2 <- OTUs2 %>%
  filter(`OTU ID` %in% unique(Negs$`OTU ID`)) %>%
  pivot_wider(names_from="Sample", values_from="Reads") %>%
  rowwise() %>%
  mutate(Sum=sum(across(2:61)), NegSum=sum(Neg1, Neg2, Neg3), Prop=NegSum/Sum)

# We will drop the 10 OTUs where 10% or more  of all reads were from
# the negative controls: 4, 15, 546, 693, 1585, 2135, 2218, 2557, 2738, 2838
Negs3tax <- left_join(Negs2, tx3) %>%
  filter(Prop>0.1)

# The 4 remaining OTUs were present in >40% of real samples, and were present
# in only 1 (n=3) or 2 (n=1) of the negative controls, in relatively low incidence.
# Keep these OTUs: 1, 5, 17, 70
Negs3 <- OTUs2 %>%
  filter(`OTU ID` %in% unique(Negs$`OTU ID`)) %>%
  group_by(`OTU ID`) %>%
  filter(Sample!="Neg1", Sample!="Neg2", Sample!="Neg3", Sample!="BioMock1", Sample!="BioMock2") %>%
  mutate(Sample_prop=(length(Sample)-length(Sample[Reads==0]))/length(Sample))

# Dropping the 10 OTUs and dropping the negative control samples
OTUse <- OTUs2 %>%
  filter(!`OTU ID` %in% unique(Negs3tax$`OTU ID`)) %>%
  filter(Sample!="Neg1", Sample!="Neg2", Sample!="Neg3")
OTUse2 <- OTUse %>%
  pivot_wider(names_from="Sample", values_from="Reads") %>%
  left_join(tx4)
# Above are the final taxonomy tables, with all necessary corrections.
# OTUse is the "long" form, OTUse2 is more compact and human friendly.

# Reformatting the final tax table to be more analysis-friendly
OTUs3 <- OTUse %>%
  pivot_wider(names_from="Sample", values_from="Reads") %>%
  left_join(tx3) %>%
  dplyr::rename(Taxonomy="Taxons", OTU=`OTU ID`) %>%
  mutate(Level=as.factor(Level)) %>%
  mutate(Level=recode_factor(Level, k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")) %>%
  group_by(OTU) %>%
  pivot_wider(names_from="Level", values_from="Taxonomy") %>%
  pivot_longer(cols=2:58, names_to="Sample", values_to="Reads")

### Data analysis begins here ###

#### Misc summary stat calculations ####

# Calculate how many reads we have for a given taxonomic level
OTUs.stat <- OTUs3 %>%
  filter(Kingdom=="Fungi") %>%   # Drop non-fungal OTUs
  group_by(Order) %>%            # Change to other tax levels for other info.
  summarize(n=sum(Reads))

# Calculate avg & median # of orders per millipede
samples2drop <- c("BioMock1", "BioMock2", "19BA14M", "19CO07F", "19CO08M") # Outliers
OTUseX <- OTUs3 %>%                       # Take whole dataset  
  filter(Kingdom == "Fungi") %>%          # Drop non-fungi
  filter(!(Sample %in% samples2drop)) %>% # Drop controls & outliers
  group_by(Order) %>%
  mutate(n = sum(Reads)) %>%
  mutate(Order = replace_na(Order, "Other")) %>%
  mutate(Order = replace(Order, Order == "incertae sedis", "Other")) %>%
  group_by(Order, Sample) %>%
  summarize(Sum = sum(Reads)) %>%
  mutate(Sum = replace(Sum, Sum > 0, 1)) %>% # Change values to 1
  group_by(Sample) %>%
  mutate(OrderCount = sum(Sum)) %>%
  select(Sample,OrderCount) %>% # Take columns Sample and GenusCount
  distinct() # Drop duplicate rows
mean(OTUseX$OrderCount) # Mean calculation
median(OTUseX$OrderCount) # Median calculation

# Check in on the biomock compositions.
biomocks <- c("BioMock1", "BioMock2")
biomocks_OTU <- OTUs3 %>%
  filter(Sample %in% biomocks) %>%
  filter(Reads != 0)
mock1 <- biomocks_OTU %>% filter(Sample == "BioMock1")
mock2 <- biomocks_OTU %>% filter(Sample == "BioMock2")
biomock_barplot <- ggplot(biomocks_OTU, aes(fill=OTU, y=Reads, x=Sample)) +
  geom_bar(position="stack", stat="identity")
biomock_stackbarplot <- ggplot(biomocks_OTU, aes(fill=OTU, y=Reads, x=Sample)) +
  geom_bar(position="fill", stat="identity")
biomock_summ <- biomocks_OTU %>%
  group_by(Sample,OTU,Family,Genus,Species) %>%
  summarize(n=sum(Reads))
biomock1_summ <- mock1 %>%
  group_by(OTU,Phylum) %>%
  summarize(mock1reads=sum(Reads), read_pc=(mock1reads/105095)*100) # sum of biomock1 reads
biomock2_summ <- mock2 %>%
  group_by(OTU,Phylum) %>%
  summarize(mock2reads=sum(Reads), read_pc=(mock2reads/154752)*100) # sum of biomock2 reads
biomock_comb <- biomock1_summ %>% cbind(biomock2_summ$mock2reads, biomock2_summ$read_pc)
colnames(biomock_comb) <- c("OTU","Phylum","M1_reads","M1_readpc","M2_reads","M2_readpc")
  
#### Alpha diversity ####

# We will use Shannon's Diversity Index.

# Get data into cooperative format.
OTUs.mat2 <- OTUs3 %>%
  filter(Kingdom == "Fungi") %>% # remove non-fungi
  group_by(OTU, Sample) %>%
  summarize(Value = sum(Reads)) %>%
  filter(Sample != "BioMock1" & Sample != "BioMock2") %>% # remove biomocks
  pivot_wider(id_cols = "Sample", names_from = "OTU", values_from = "Value") %>%
  column_to_rownames("Sample")

# Calc SDI at the sample level.
a_sha_sample2 <- read.csv("SDI_by_sample.csv", header = T) %>%
  filter(Sample != "19CO08M")  # Outlier
a_sha_sample <- diversity(OTUs.mat2, index="shannon")

# Prepare file similar to OTUs3, but w/ metadata on the end.
metadata2 <- metadata %>%
  filter(SampleID != "BioMock1") %>%
  filter(SampleID != "BioMock2")
collist <- metadata2$SampleID
OTUs.mat3 <- cbind(collist, OTUs.mat2)
names(OTUs.mat3)[1] <- "SampleID"
OTUswMeta <- left_join(OTUs.mat3, metadata)

# Calculate SDI for sites (at OTU level).
OTUsbySite3 <- OTUswMeta[,1:2865] %>% 
  filter(SampleID != "19CO07F") %>%  # Outlier
  filter(SampleID != "19BA14M") %>%  # Outlier
  filter(SampleID != "19CO08M") %>%  # Outlier
  subset(select=-c(SampleID)) %>%
  pivot_longer(cols = 1:2863, names_to = "OTU", values_to = "Reads") %>%
  group_by(Site, OTU) %>%
  summarize(Value = sum(Reads)) %>%
  pivot_wider(id_cols = "Site", names_from = "OTU", values_from = "Value")
collist <- OTUsbySite3$Site # Next 3 lines are my messy way of renaming tibble rows.
OTUsbySite3$Site = NULL
rownames(OTUsbySite3) <- collist
a_sha_site <- diversity(OTUsbySite3, index="shannon")

# SDI at the colony level
OTUsbyColony3 <- cbind(OTUswMeta[,2866], OTUswMeta[1:2864]) %>%
  filter(SampleID != "19CO07F") %>%  # Outlier
  filter(SampleID != "19BA14M") %>%  # Outlier
  filter(SampleID != "19CO08M") %>%  # Outlier
  subset(select=-c(SampleID))
names(OTUsbyColony3)[1] <- "Colony"
OTUsbyColony4 <- OTUsbyColony3 %>%
  pivot_longer(cols=2:2864, names_to="OTU", values_to="Reads") %>%
  group_by(Colony, OTU) %>%
  summarize(Value=sum(Reads)) %>%
  pivot_wider(id_cols="Colony", names_from="OTU", values_from="Value")
collist <- OTUsbyColony4$Colony
OTUsbyColony4$Colony = NULL
rownames(OTUsbyColony4) <- collist
a_sha_col <- diversity(OTUsbyColony4, index="shannon")

# SDI by sex
OTUsbySex <- cbind(OTUswMeta[,2867], OTUswMeta[1:2864]) %>%
  filter(SampleID != "19CO07F") %>%  # Outlier
  filter(SampleID != "19BA14M") %>%  # Outlier
  filter(SampleID != "19CO08M") %>%  # Outlier
  subset(select=-c(SampleID))
names(OTUsbySex)[1] <- "Sex"
OTUsbySex2 <- OTUsbySex %>%
  pivot_longer(cols=2:2864, names_to="OTU", values_to="Reads") %>%
  group_by(Sex, OTU) %>%
  summarize(Value=sum(Reads)) %>%
  pivot_wider(id_cols="Sex", names_from="OTU", values_from="Value")
collist <- OTUsbySex2$Sex
OTUsbySex2$Sex = NULL
rownames(OTUsbySex2) <- collist
a_sha_sex <- diversity(OTUsbySex2, index="shannon")

# To statistically compare Shannon indices, we need multiple samples (=individuals),
# so we can only do this for sites & sex. Colony is possible, in theory, but only
# for a handful of special comparisons, like col 6 and 7 (both large colonies
# from 1 site), or 11 and 1 (both n=4, from diff sites).

# Check for normality and equality of variance before analysis.

# Sortable table with alpha diversity and metadata.
holder <- cbind(OTUswMeta[,2865:2867], OTUswMeta[1:2864]) %>%
  filter(SampleID != "19CO07F") %>%  # Outlier
  filter(SampleID != "19BA14M") %>%  # Outlier
  filter(SampleID != "19CO08M")      # Outlier
names(holder)[4] <- "Sample"
a_sha_sample5 <- left_join(a_sha_sample2, holder[,1:4], by="Sample")

# We will use a Welch's ANOVA (rstatix) & Games-Howell posthoc test (rstatix)
# to test if sites are significantly different in terms of SDI, and a standard
# t-test (rstatix) to test if sexes are significantly different in terms of SDI.

# Sites:
alphasite <- welch_anova_test(a_sha_sample5, SDI ~ Site)
alphasite2 <- oneway.test(SDI ~ Site, data = a_sha_sample5, var.equal = FALSE)
omega_squared(alphasite2)
games_howell_test(a_sha_sample5, SDI ~ Site)
# Sites are significantly different in terms of SDI.
# Specifically, CH and CO are different, but BA-CO are very close as well.

# Sexes:
t_test(a_sha_sample5, SDI ~ Sex, paired = F, var.equal = T)
# Sexes are not significantly different in terms of SDI.

# For colonies there are two sets of comparisons: 
# 1) Compare all of col 6 and all of col 7: 2 big colonies from the same site.
# 2) Compare all cols with n=4 (drop individuals from BA megacolonies
#    until n=4). Make all comparisons if sig diff.

# Setting up tables for each comparison.

# Comparison 1 (All indivs in cols 6 and 7):
a_sha_sample_comp1 <- a_sha_sample5 %>%
  filter(Site != "CH") %>%
  filter(Site != "CO")

# Comparison 2 (n=4 cols, dropping randomly from BA megacolonies):
sample(c("1", "2", "4", "5", "7", "8", "9", "10"), 4, replace=FALSE) # Col 6
sample(c("11", "12", "13", "15", "16", "17", "18", "19", "20"), 4, replace=FALSE) # Col 7
# Keep 2, 4, 5, 8 for col6, and 11, 15, 16, 20 for col7.
keepers <- c("19BA02M","19BA04M","19BA05M","19BA08F","19BA11F","19BA15M","19BA16F","19BA20F")
n4colonies <- c("11", "12", "1", "5")
nonBAindivs <- a_sha_sample5 %>%
  filter(Site != "BA") %>%
  subset(Colony %in% n4colonies)
a_sha_sample_comp2 <- a_sha_sample5 %>%
  filter(Site != "CH") %>%
  filter(Site != "CO") %>%
  subset(Sample %in% keepers)
a_sha_sample_comp2 <- rbind(a_sha_sample_comp2, nonBAindivs)

# Remember to check ANOVA assumptions.

# We will use a Wilcoxon test (rstatix) to test if 2 large colonies from the
# same site are significantly different in terms of SDI.
wilcox_test(SDI ~ Colony, data = a_sha_sample_comp1)
# Colonies 6 and 7 are not significantly different.

# Check normality and equality of variance for Comparison 2 dataset.

# We will use a one-way ANOVA (rstatix) to test if colonies
# of the same size are significantly different.
# Note the design is nested, even though I'm not interested in site here.
oneway_anova <- anova_test(SDI ~ Colony + Colony %in% Site, data = a_sha_sample_comp2, effect.size = 'ges')
get_anova_table(oneway_anova)
## Colonies of the same size are not significantly different.

# Is millipede weight correlated with alpha diversity?
# Let's start with a scatter plot, then try Pearson correlation.
metadata3 <- metadata2 %>%
  filter(SampleID != "19CO07F") %>%  # Outliers
  filter(SampleID != "19BA14M") %>%  # Outliers
  filter(SampleID != "19CO08M")      # Outliers
shapiro.test(a_sha_sample2$SDI)
shapiro.test(metadata3$DryWeight)
meta <- metadata3
colnames(meta)[1]="Sample"
weight <- merge(as.data.frame(a_sha_sample2), meta, by = "Sample")
plot(weight$SDI, weight$DryWeight)
cor.test(weight$SDI, weight$DryWeight, method="pearson")
# Weight and SDI are not correlated.

#Is DNA concentration correlated with weight?
shapiro.test(metadata3$DNAconc)
plot(weight$DryWeight, weight$DNAconc)
# Weight and DNA concentration are correlated, unsurprisingly.

# How about the same question on subset of males? Females?
fweight <- filter(weight, Sex == "F")
mweight <- filter(weight, Sex == "M")
plot(fweight$SDI, fweight$DryWeight)
plot(mweight$SDI, mweight$DryWeight)
cor.test(fweight$SDI, fweight$DryWeight, method="pearson")
cor.test(mweight$SDI, mweight$DryWeight, method="pearson")
# SDI is not correlated with weight for males or females.

#### Beta diversity (stats) ####

# We will use Bray-Curtis (abundance-based) & Bsim (presence-absence based). 
# Bsim is Simpson dissimilarity (= turnover component of Sørensen dissimilarity).

# Drop outliers & make a presence-absence version of the dataset.
OTUs.mat4 <- OTUs.mat2[-c(12, 44, 45),]
OTUs.mat4.presabs <- OTUs.mat4                  # Duplicate the dataset
OTUs.mat4.presabs[OTUs.mat4.presabs > 0] <- 1   # Convert to presence/absence
OTUs.mat4.presabs <- OTUs.mat4.presabs %>%      # Drop OTUs with 0 occurrences
  select_if(colSums(.) != 0)  # Note that in presabs, row 1 = BA, 2 = CH, 3 = CO
meta <- metadata %>%
  filter(!(SampleID %in% samples2drop))
OTUs.mat4.presabs2 <- cbind(OTUs.mat4.presabs, meta)
rownames(OTUs.mat4.presabs2) <- NULL
OTUs.mat4.presabs3 <- OTUs.mat4.presabs2 %>%
  pivot_longer(cols=1:2848, names_to="OTU", values_to="Reads") %>%
  group_by(Site, OTU) %>%
  summarize(Value = sum(Reads)) %>%
  pivot_wider(id_cols = "Site", names_from = "OTU", values_from = "Value") 
collist <- OTUs.mat4.presabs3$Site
OTUs.mat4.presabs3$Site <- NULL
rownames(OTUs.mat4.presabs3) <- collist # Format of OTUs.mat4.presabs3 now == OTUswMeta.
OTUs.mat4.presabs4 <- OTUs.mat4.presabs3 # Converting back to pres/abs.
OTUs.mat4.presabs4[OTUs.mat4.presabs4 > 0] <- 1

# Individual distances:
bray_dist_indiv <- vegdist(OTUs.mat4, "bray")
bsim_dist_indiv <- betadiver(OTUs.mat4.presabs, method = "w")

# What are the beta diversity distances between our sites?
# Tidy table (betacomps) of distances between sites:
beta_braycurtis_site <- vegdist(OTUsbySite3, "bray")
beta_bsim_site <- betadiver(OTUs.mat4.presabs4, method = "w")
comps <- c("BA-CH", "BA-CO", "CO-CH")
betacomps <-rbind(comps, beta_braycurtis_site, beta_bsim_site)

# Are the sites significantly different, based on Bray & Bsim?
# Which statistical test, ANOVA/Tukey HSD or PERMANOVA/Kruskal-Wallis?
# Microbial data is rarely normally distributed, so we should use the latter.
# (Hugerth & Anderson, 2017)

# We will use PERMANOVA (adonis2, vegan) and "pairwise multilevel 
# comparison using adonis" (pkg pairwiseAdonis)

# Bray-Curtis:
OTUdist_bray <- vegdist(OTUs.mat4, method = "bray")
betadispSite_bray <- betadisper(d = OTUdist_bray, group = metadata3$Site, bias.adjust=TRUE)
scores(betadispSite_bray, display="centroids") 
adonis2(OTUdist_bray ~ metadata3$Site, permutations = 9999)
pairwise.adonis(x = OTUdist_bray, factors = metadata3$Site, p.adjust.m = "bonferroni")

# Bsim:
OTUdist_bsim <- betadiver(OTUs.mat4.presabs, method="w")
betadispSite_bsim <- betadisper(d = OTUdist_bsim, group = meta$Site, bias.adjust = TRUE)
scores(betadispSite_bsim, display="centroids") 
adonis2(OTUdist_bsim ~ meta$Site, permutations = 9999)
pairwise.adonis(x = OTUdist_bsim, factors = meta$Site, p.adjust.m = "bonferroni")

# What are the distances between males and females, using Bray-Curtis & bsim?
# Are they statistically different?

# Bray-Curtis:
betadispSex_bray <- betadisper(d = OTUdist_bray, group = metadata3$Sex, bias.adjust=TRUE)
scores(betadispSex_bray, display="centroids")
adonis2(OTUdist_bray ~ metadata3$Sex, permutations = 9999)

# Bsim:
betadispSex_bsim <- betadisper(d = OTUdist_bsim, group = meta$Sex, bias.adjust=TRUE)
scores(betadispSex_bsim, display="centroids")
adonis2(OTUdist_bsim~ meta$Sex, permutations = 9999)

# What are the beta diversity distances between all colonies? 
# Note, no statistical testing here, just collecting data. Can't
# do stats on all cols because col size varies dramatically.

# Need to modify OTUsbyColony4 to pres/abs first.
OTUsbyColony4.presabs <- OTUsbyColony4[-2,]             # Dup dataset, drop bad colony
OTUsbyColony4.presabs[OTUsbyColony4.presabs > 0] <- 1   # Convert to pres/abs
OTUsbyColony4.presabs <- OTUsbyColony4.presabs %>%      # Drop OTUs w/ 0 occurrences
  select_if(colSums(.) != 0)
rownames(OTUsbyColony4.presabs) <- c(1,3,4,5,6,7,8,9,10,11,12,14,15)

# Distances between sites, by metric:
beta_bray_col <- vegdist(OTUsbyColony4, "bray")
beta_bsim_col <- betadiver(OTUsbyColony4.presabs, method = "w")

# What are the beta diversity distances between the below comparisons?
# 1) Compare all of col 6 and all of col 7: 2 big colonies from the same site.
# 2) Compare all colonies with n=4 (drop individuals from BA megacolonies
#    until n=4). Make all comparisons if sig diff.

# Comparison 1:

# Bray-Curtis:
meta2 <- meta[1:17,]
dist_bray_comp1 <- vegdist(OTUs.mat4[1:17,], "bray")
betadispCol_bray_comp1 <- betadisper(d = dist_bray_comp1, group = meta2$Colony, bias.adjust=TRUE)
adonis2(dist_bray_comp1 ~ meta2$Colony, permutations = 9999)

# Bsim:
dist_bsim_comp1 <- betadiver(OTUs.mat4.presabs[1:17,], method = "w")
betadispCol_bsim_comp1 <- betadisper(d = dist_bsim_comp1, group = meta2$Colony, bias.adjust=TRUE)
adonis2(dist_bsim_comp1 ~ meta2$Colony, permutations = 9999)

# Comparison 2:

# Bray-Curtis:
OTUs.mat4_comp2 <- cbind(meta, OTUs.mat4) %>% # Add metadata to make filtering easy
  subset(SampleID %in% keepers) %>% # Pull out randomly selected BA indivs
  rbind(cbind(meta, OTUs.mat4) %>% filter(Site != "BA") %>% subset(Colony %in% n4colonies)) # Pull out indivs from other sites
OTUs.mat4_comp2 <- OTUs.mat4_comp2[,8:2870] # Trim off metadata
meta_comp2 <- meta %>%
  subset(SampleID %in% keepers) %>%
  rbind(meta %>% filter(Site != "BA") %>% subset(Colony %in% n4colonies))
dist_bray_comp2 <- vegdist(OTUs.mat4_comp2, "bray")
betadispCol_bray_comp2 <- betadisper(d = dist_bray_comp2, group = meta_comp2$Colony, bias.adjust=TRUE)
adonis2(dist_bray_comp2 ~ meta_comp2$Colony, permutations = 9999)

# Bsim:
OTUs.mat4.presabs_comp2 <- cbind(meta, OTUs.mat4.presabs) %>%
  subset(SampleID %in% keepers) %>%
  rbind(cbind(meta, OTUs.mat4.presabs) %>% filter(Site != "BA") %>% subset(Colony %in% n4colonies))
OTUs.mat4.presabs_comp2 <- OTUs.mat4.presabs_comp2[,8:2855]
dist_bsim_comp2 <- betadiver(OTUs.mat4.presabs_comp2, method = "w")
betadispCol_bsim_comp2 <- betadisper(d = dist_bsim_comp2, group = meta_comp2$Colony, bias.adjust=TRUE)
adonis2(dist_bsim_comp2 ~ meta_comp2$Colony, permutations = 9999)

#### Beta diversity (ordinations) ####

# Bray-Curtis based NMDS

# Base analysis
Bray_NMDS <- metaMDS(comm = OTUs.mat4, distance = "bray", k = 2) # If you rerun this line, update best_soln_stress!
Shepard <- stressplot(Bray_NMDS)
best_soln_stress <- 0.177 # updated 3-1-22

# Validate Bray NMDS using simplified permutation-based method from 
# Dexter et al. 2018.
# Console only holds 1000 lines, so do 500 iterations 2x, and paste stress
# results together manually into a csv. "Text to Columns" and filtering in
# Excel makes csv prep fairly easy.
Bray_NMDSp1 <- metaMDS(comm = OTUs.mat4, distance = "bray", k = 2, try = 500,
                       trymax = 501, engine = "monoMDS", autotransform = T)
Bray_NMDSp2 <- metaMDS(comm = OTUs.mat4, distance = "bray", k = 2, try = 500,
                       trymax = 501, engine = "monoMDS", autotransform = T)
Bray_NMDS_stress <- read.csv("NMDS_stress.csv", header = T)
sigmaX <- sd(Bray_NMDS_stress$Stress)
# one sample z test via BSDA package
z.test(Bray_NMDS_stress$Stress, alternative = 'two.sided', mu = best_soln_stress, sigma.x = sigmaX, conf.level = .95)
# NMDS is reliable!

# Generating Bray NMDS figure
scores <- as.data.frame(scores(Bray_NMDS))
scores$Site <- metadata3$Site # Set up possible comparisons, tho we will likely
scores$Sex <- metadata3$Sex   #     only use site.
scores$Colony <- as.factor(metadata3$Colony)
lab <- paste0("stress = ", best_soln_stress)
bray_nmds_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS2, shape = Site, color = Site),size = 4) +
  scale_color_manual(values=c('#009E73','#0072B2', '#D55E00'))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  annotate(geom="text", x=0.8, y=-0.7, label = lab, color = "black", size = 6) +
  stat_ellipse(aes(x = scores$NMDS1, y = scores$NMDS2, lty = scores$Site), show.legend = T, level = 0.95) +
  scale_shape_manual(values = c(15, 16, 17), name = "Site", breaks = c("BA", "CH", "CO"), labels = c("BA", "CH", "CO")) 
bray_nmds_plot

# Bsim based NMDS

# Base analysis
#bsim_dist <- betadiver(OTUs.mat4.presabs, method = "w")
bsim_NMDS <- metaMDS(comm = OTUs.mat4.presabs, distfun = betadiver, distance = "w", k = 2)
bsim_Shep <- stressplot(bsim_NMDS)
best_soln_stress2 <- 0.164 # Updated 6-21-22

# Validate Bsim NMDS
Bsim_NMDSp1 <- metaMDS(comm = OTUs.mat4.presabs, distfun = betadiver, distance = "w", k = 2, try = 500,
                       trymax = 501, engine = "monoMDS", autotransform = T)
Bsim_NMDSp2 <- metaMDS(comm = OTUs.mat4.presabs, distfun = betadiver, distance = "w", k = 2, try = 500,
                       trymax = 501, engine = "monoMDS", autotransform = T)
Bsim_NMDS_stress <- read.csv("NMDS_stress2.csv", header = T)
sigmaX <- sd(Bsim_NMDS_stress$Stress)
z.test(Bsim_NMDS_stress$Stress, alternative = 'two.sided', mu = best_soln_stress2, sigma.x = sigmaX, conf.level = .95)
# NMDS is reliable!

# Generating Bsim NMDS figure
scores <- as.data.frame(scores(bsim_NMDS))
scores$Site <- metadata3$Site
scores$Sex <- metadata3$Sex
scores$Colony <- as.factor(metadata3$Colony)
lab <- paste0("stress = ", best_soln_stress2)
bsim_nmds_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS2, shape = Site, color = Site),size = 4) +
  scale_color_manual(values=c('#009E73','#0072B2', '#D55E00'))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  annotate(geom="text", x=1, y=-0.8, label = lab, color = "black", size = 6) +
  stat_ellipse(aes(x = scores$NMDS1, y = scores$NMDS2, lty = scores$Site), show.legend = T, level = 0.95) +
  scale_shape_manual(values = c(15, 16, 17), name = "Site", breaks = c("BA", "CH", "CO"), labels = c("BA", "CH", "CO"))
bsim_nmds_plot

# Bray-Curtis based PCoA

bray_pcoa <- pcoa(OTUdist_bray)
axes <- bray_pcoa$vectors[,1:2]
meta_axes <- cbind(metadata3, axes)
#head(brachy_pcoa) # Run to see the % variation explained by ea axis
bray_pcoa_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(data = meta_axes, aes(x = Axis.1, y = Axis.2, shape = Site, color = Site),size = 4) +
  xlab("PCo1 (14.5%)") +
  ylab("PCo2 (10.2%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  stat_ellipse(aes(x = meta_axes$Axis.1, y = meta_axes$Axis.2, lty = meta_axes$Site), show.legend = T, level = 0.95) + # Turn off for no ellipses, modify for sex or col
  scale_shape_manual(values = c(15, 16, 17), name = "Site", breaks = c("BA", "CH", "CO"), labels = c("BA", "CH", "CO")) +  # Modify this if sex or col is desired
  scale_color_manual(values=c('#009E73','#0072B2', '#D55E00'))
# Below line not used, but saved as a template for using ggsave.
#ggsave("Bdiv_Bray_PCoA.jpeg", plot = bray_pcoa_plot, device = "jpeg", width = 4, height = 3, units = "in", dpi = 600)

# Bsim based PCoA
bsim_pcoa <- pcoa(OTUdist_bsim)
axes <- bsim_pcoa$vectors[,1:2]
meta_axes <- cbind(metadata3, axes)
#head(bray_pcoa) # Run to see the % variation explained by ea axis
bsim_pcoa_plot <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(data = meta_axes, aes(x = Axis.1, y = Axis.2, shape = Site, color = Site),size = 4) +
  xlab("PCo1 (14.9%)") +
  ylab("PCo2 (10.5%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  stat_ellipse(aes(x = meta_axes$Axis.1, y = meta_axes$Axis.2, lty = meta_axes$Site), show.legend = T, level = 0.95) + # Turn off for no ellipses, modify for sex or col
  scale_shape_manual(values = c(15, 16, 17), name = "Site", breaks = c("BA", "CH", "CO"), labels = c("BA", "CH", "CO")) +  # Modify this if sex or col is desired
  scale_color_manual(values=c('#009E73','#0072B2', '#D55E00'))
bsim_pcoa_plot

#### 100% stacked bar plot ####

# Setting up data for the plot
bar.dat <- OTUs3 %>%
  filter(Kingdom == "Fungi") %>%
  group_by(Class) %>%
  mutate(n = sum(Reads)) %>%
  mutate(Class = replace_na(Class, "Other")) %>%
  mutate(Class = replace(Class, Class == "incertae sedis", "Other")) %>%
  mutate(Class = replace(Class, Class == "No Hit", "Other")) %>%
  mutate(Class = replace(Class, n<1000, "Low")) %>%
  group_by(Class, Sample) %>% 
  summarize(Sum = sum(Reads))
colnames(bar.dat) <- c("Class", "SampleID", "Sum")
bar.meta <- left_join(bar.dat, metadata3) %>%
  filter(!(SampleID %in% samples2drop))
bar.dat.stat <- bar.dat %>%
  group_by(Class) %>% # Change "Class" to diff lvl if desired
  summarize(Sum = sum(Sum))
# 18 bins, 16 classes, 2 others:
# Low = OTU w/ <1000 reads, Other = not IDed to class

# Organizing the order of the classes in the plot.
seq <- bar.dat.stat$Class[c(9,13,2,10,16,4,11,14,3,7,18,12,17,6,1,8,15,5)]
bar.meta$Class=as.factor(bar.meta$Class)
bar.meta$Class=factor(bar.meta$Class, levels=seq)

# Collapsing individuals into colonies
collapse.bar = bar.meta %>%
  group_by(Class, Site, Colony) %>%
  summarize(Sum = sum(Sum))
collapse.bar$Colony <- as.factor(collapse.bar$Colony)

# revising colony names to include n
collapse.bar2 <- collapse.bar
collapse.bar2$Colony <- as.character(collapse.bar2$Colony)
collapse.bar2 <- collapse.bar2 %>%
  mutate(Colony = replace(Colony, Colony == "6", "06 (8)")) %>%
  mutate(Colony = replace(Colony, Colony == "7", "07 (9)")) %>%
  mutate(Colony = replace(Colony, Colony == "8", "08 (2)")) %>%
  mutate(Colony = replace(Colony, Colony == "9", "09 (2)")) %>%
  mutate(Colony = replace(Colony, Colony == "10", "10 (3)")) %>%
  mutate(Colony = replace(Colony, Colony == "11", "11 (4)")) %>%
  mutate(Colony = replace(Colony, Colony == "12", "12 (4)")) %>%
  mutate(Colony = replace(Colony, Colony == "14", "14 (2)")) %>%
  mutate(Colony = replace(Colony, Colony == "15", "15 (3)")) %>%
  mutate(Colony = replace(Colony, Colony == "1", "01 (4)")) %>%
  mutate(Colony = replace(Colony, Colony == "2", "02 (1)")) %>%
  mutate(Colony = replace(Colony, Colony == "3", "03 (3)")) %>%
  mutate(Colony = replace(Colony, Colony == "4", "04 (3)")) %>%
  mutate(Colony = replace(Colony, Colony == "5", "05 (4)"))
collapse.bar2$Site <- factor(collapse.bar2$Site, levels = c("CO", "BA", "CH"))

# Assembling the plot
pal <- c("#FFFF00","#004D43","#7A4900","#A30059","#809693","#997D87","#8FB0FF","#63FFAC","#008941","#4FC601","#1CE6FF","#B79762","#FF34FF","#006FA6","#FF4A46","#3B5DFF","#0000A6","#5A0007")
brachy_stackbar <- ggplot(collapse.bar2, aes(x = Colony, y = Sum, fill = Class)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = pal) +
  ylab("Proportion of reads") +
  xlab("Colony (n individuals)") +
  facet_grid(~Site, scales = "free_x", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "white", colour = "black"))
brachy_stackbar

#### Species accumulation curve ####

spec_accum_CI <- specaccum(OTUs.mat4, method="random", permutations=100)
spec_accum_CI_plot <- plot(spec_accum_CI$sites, spec_accum_CI$richness,
                           xlab="Individuals sampled", ylab="Species richness")

#### FungalTraits & FunGuild ####

# Data modifications for FungalTraits analysis to follow
# Eric Morrison's method:
# https://github.com/ewmorr/FungalTraits_Polme

# CB taxonomy & processing:
CultureBased_Taxonomy <- CultureBased[,c(1,300)] %>%
  separate(Taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="; ")
CultureBased_Taxonomy$Kingdom <- gsub("k__", "", CultureBased_Taxonomy$Kingdom)
CultureBased_Taxonomy$Phylum <- gsub("p__", "", CultureBased_Taxonomy$Phylum)
CultureBased_Taxonomy$Class <- gsub("c__", "", CultureBased_Taxonomy$Class)
CultureBased_Taxonomy$Order <- gsub("o__", "", CultureBased_Taxonomy$Order)
CultureBased_Taxonomy$Family <- gsub("f__", "", CultureBased_Taxonomy$Family)
CultureBased_Taxonomy$Genus <- gsub("g__", "", CultureBased_Taxonomy$Genus)
CultureBased_Taxonomy$Species <- gsub("s__", "", CultureBased_Taxonomy$Species)

# CF taxonomy & processing
CultureFree_Taxonomy <- CultureFree[,c(1,57)] %>%
  separate(Taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=", ")
CultureFree_Taxonomy$Kingdom <- gsub("k_", "", CultureFree_Taxonomy$Kingdom)
CultureFree_Taxonomy$Phylum <- gsub("p_", "", CultureFree_Taxonomy$Phylum)
CultureFree_Taxonomy$Class <- gsub("c_", "", CultureFree_Taxonomy$Class)
CultureFree_Taxonomy$Order <- gsub("o_", "", CultureFree_Taxonomy$Order)
CultureFree_Taxonomy$Family <- gsub("f_", "", CultureFree_Taxonomy$Family)
CultureFree_Taxonomy$Genus <- gsub("g_", "", CultureFree_Taxonomy$Genus)
CultureFree_Taxonomy$Species <- gsub("s_", "", CultureFree_Taxonomy$Species)
CultureFree_Taxonomy[CultureFree_Taxonomy == "incertae sedis"] <- NA

# Creating a dataframe of FungalTraits columns of interest.
lifestyles = data.frame(
  Genus.join = FungalTraits$GENUS,
  primary = FungalTraits$primary_lifestyle,
  secondary = FungalTraits$Secondary_lifestyle,
  lifestyle_c = FungalTraits$Comment_on_lifestyle_template,
  endo_capac = FungalTraits$Endophytic_interaction_capability_template,
  plpa_capac = FungalTraits$Plant_pathogenic_capacity_template,
  decay_subst = FungalTraits$Decay_substrate_template,
  decay_type = FungalTraits$Decay_type_template,
  animal_intrxn = FungalTraits$Animal_biotrophic_capacity_template,
  hosts = FungalTraits$Specific_hosts,
  stringsAsFactors = F)

# Running CB thru FT:
# New column to join by (without losing taxonomy columns)
CultureBased_Taxonomy$Genus.join <- CultureBased_Taxonomy$Genus
CultureBasedwFungalTraits <- left_join(CultureBased_Taxonomy, lifestyles, by="Genus.join")
# Replace blanks w/ NA for consistency
CultureBasedwFungalTraits[CultureBasedwFungalTraits == ""] <- NA
# Getting OTUs and final table
CultureBased_OTUs <- CultureBased[,1:299]
CultureBasedwFungalTraits2 <- left_join(CultureBased_OTUs, CultureBasedwFungalTraits, by="Isolate")

# Running CF thru FT:
CultureFree_Taxonomy$Genus.join <- CultureFree_Taxonomy$Genus
CultureFreewFungalTraits <- left_join(CultureFree_Taxonomy, lifestyles, by="Genus.join")
CultureFreewFungalTraits[CultureFreewFungalTraits == ""] <- NA
CultureFree_OTUs <- CultureFree[,1:56]
CultureFreewFungalTraits2 <- left_join(CultureFree_OTUs, CultureFreewFungalTraits, by="OTU ID")

# FunGuild analysis:
# Lets isolate the Guild column from FunGuild results.

# CB:
CultureBasedwFunGuild_Guilds <- data.frame(
  Isolate = CultureBasedwFunGuild$Isolate,
  Guild = CultureBasedwFunGuild$Guild,
  stringsAsFactors = F)
CultureBasedwFunGuild_Guilds_Meta <- merge(CultureBasedwFunGuild_Guilds, CultureBasedMeta) # Add metadata
# CF:
CultureFreewFunGuild_Guilds <- data.frame(
  OTU = CultureFreewFunGuild$`OTU ID`,
  Guild = CultureFreewFunGuild$Guild,
  stringsAsFactors = F)
# CF metadata needs more modification, & guilds need revision.
CultureFreewFunGuild_Guilds_Meta <- merge(CultureFreewFunGuild_Guilds, CultureFreeMeta)
CFwFG_CorrGuilds <- CultureFreewFunGuild_Guilds %>%
  separate(Guild, into=as.character(c(1:15)), sep="-") %>%
  pivot_longer(cols=as.character(1:15)) %>%
  filter(!is.na(value))
colnames(CFwFG_CorrGuilds)[3] = "OldGuild"

CFSeqCounts <- pivot_longer(CultureFree[,1:56], !`OTU ID`, names_to="Millipede", values_to="SeqCount")
CFwFG_Meta <- CultureFreewFunGuild_Guilds_Meta %>%
  separate(Guild, into = as.character(c(1:15)), sep = "-") %>%
  pivot_longer(cols = as.character(1:15)) %>%
  filter(!is.na(value))
CFwFG <- merge(CFwFG_Meta, CFSeqCounts, by.x="OTU", by.y="OTU ID") %>% # Slow, wait!
  filter(Millipede == SampleID) %>%
  subset(select = -c(SampleID, ExtGroup, DryWeight, DNAconc, name))
CFwFG2 <- CFwFG[,c(6,1,5,7,2,3,4)] 
colnames(CFwFG2)[3] = "OldGuild"  
Shifts <- data.frame( # Look here to see which old guilds were binned where.
  OldGuild = unique(CFwFG_CorrGuilds$OldGuild),
  NewGuild = c("Plant symbiont","Saprotroph","Epiphyte","Animal associate",
               "Fungal parasite","Plant pathogen","Plant saprotroph",
               "Plant saprotroph","Plant symbiont","Lichen parasite",
               "Saprotroph","Plant parasite","Plant parasite","Plant saprotroph",
               "Saprotroph","DROP","Animal associate","Lichenized",
               "Plant symbiont","Animal associate","Plant symbiont",
               "Animal associate","Plant saprotroph","Plant parasite",
               "Animal associate","Plant symbiont","Plant saprotroph",
               "Plant symbiont","Plant symbiont","DROP"),
              stringsAsFactors = F)
CFwFG3 <- merge(CFwFG2, Shifts) 

# Use "separate" to break up the guild-string into units.
CultureBasedwFunGuild_Guilds2 <- CultureBasedwFunGuild_Guilds %>%
  separate(Guild, into = as.character(c(1:15)), sep = "-") %>% # overshoot
  pivot_longer(cols = as.character(1:15)) %>%
  filter(!is.na(value)) %>%
  group_by(value) %>%
  summarize(Count = length(value), Prop = Count/595) %>% # N isolates returned = 595
  mutate(Method = "CB")
CultureFreewFunGuild_Guilds2 <- CultureFreewFunGuild_Guilds %>%
  separate(Guild, into = as.character(c(1:15)), sep = "-") %>%
  pivot_longer(cols = as.character(1:15)) %>%
  filter(!is.na(value)) %>%
  group_by(value) %>%
  summarize(Count = length(value), Prop = Count/1009) %>% # N OTUs returned = 1009
  mutate(Method = "CF")
All_Guilds <- rbind(CultureBasedwFunGuild_Guilds2, CultureFreewFunGuild_Guilds2)
# Need a smaller number of guilds. Merging some & fixing bugs:
All_Guilds2 <- All_Guilds
All_Guilds2$RevisedGuild <- c("Plant parasite","Animal associate","Animal associate",
                              "Plant parasite","Plant symbiont","Saprotroph",
                              "Plant symbiont","Plant symbiont","DROP","Epiphyte",
                              "Plant symbiont","Fungal parasite","Animal associate",
                              "Animal associate","Plant saprotroph","Lichen parasite",
                              "Lichenized","Plant saprotroph","Plant pathogen",
                              "Plant saprotroph","Saprotroph","Saprotroph",
                              "Plant saprotroph","Plant parasite","Animal associate",
                              "Animal associate","Plant symbiont","Plant parasite",
                              "Plant symbiont","Saprotroph","Plant symbiont",
                              "Plant symbiont","DROP","Epiphyte","Plant symbiont",
                              "Fungal parasite","Animal associate","Animal associate",
                              "Plant saprotroph","Lichen parasite","Lichenized",
                              "Plant saprotroph","Animal associate","DROP","Plant symbiont",
                              "Plant parasite","Plant pathogen","Plant saprotroph",
                              "Plant symbiont","Saprotroph","Plant saprotroph",
                              "Saprotroph","Plant saprotroph")
All_Guilds2 <- All_Guilds2 %>%
  filter(RevisedGuild!="DROP")
All_Guilds2[31,2] <- 153
# Merge by RevisedGuild
All_Guilds3 <- All_Guilds2 %>%
  group_by(RevisedGuild, Method) %>%
  summarize(NumOTUs = sum(Count), Prop2 = sum(Prop), Method=Method) %>%
  distinct()

# Remember that these guild assignments can overlap--so, a
# visualization like a 100% stacked bar doesn't make sense!
# We'll instead use a scatterplot, with % of millipedes and % of
# isolates/OTUs as the axes, and dot size from % of genera.

# At one point we were using unassigneds, but later dropped.
# I filtered from the same table later though, so the below blocks are kept.

# Culture-Free:

# Go back to raw data, as we filtered unassigneds out in the orig import step.
CF_Unassgn <- read.delim("CultureFree_FunGuild.txt", header=T, check.names=F) %>%
  filter(Guild == "-") %>%
  subset(select = -c(Taxonomy, Taxon, `Taxon Level`, `Trophic Mode`,`Growth Morphology`, Trait, `Confidence Ranking`, Notes, `Citation/Source`))
CF_Unassgn <- CF_Unassgn[,1:56] %>%
  pivot_longer(!`OTU ID`, names_to="Millipede", values_to="SeqCount")
colnames(CF_Unassgn)[2] = "SampleID"
CF_Unassgn2 <- merge(CF_Unassgn, CultureFreeMeta) %>%
  subset(select = -c(ExtGroup, DryWeight, DNAconc))
CF_Unassgn2$Guild <- "Unassigned"
colnames(CFwFG3)[2] = "SampleID"  # cleaning up cols
colnames(CFwFG3)[3] = "OTU ID"    # cleaning up cols
colnames(CFwFG3)[8] = "Guild"     # cleaning up cols
CFwFG3 <- CFwFG3[-1]
Full_CFwFG <- rbind(CFwFG3, CF_Unassgn2) %>%
  filter(Guild != "DROP")

# Culture-Based:

# Re-import data without filtering. Drop unneeded cols.
CBwFG <- read.delim("CultureBased_FunGuild.txt", header=T, check.names=F) %>%
  subset(select = -c(Taxonomy, Taxon, `Taxon Level`, `Trophic Mode`, `Growth Morphology`, Trait, `Confidence Ranking`, Notes, `Citation/Source`)) #Drop unneeded cols
# Add in metadata, drop unneeded cols again
CBwFG2 <- merge(CBwFG, CultureBasedMeta, by.x = "Isolate", by.y= "Isolate") %>%
  subset(select = -c(Sample, Qced, Coords))
CBwFG2["Guild"][CBwFG2["Guild"] == "-"] <- "Unassigned"
# CBwFG2 contains all the data, even unassigned guilds.
# Metadata, individual, isolate #, and counts.
# Next we need to split up the Guild-strings & go long.
CBwFG3 <- CBwFG2 %>%
  separate(Guild, into = as.character(c(1:15)), sep = "-") %>%
  pivot_longer(cols = as.character(1:15)) %>%
  filter(!is.na(value))
CBwFG4 <- pivot_longer(CBwFG3[,1:299], !Isolate, names_to="Individual", values_to="SeqCount") %>%
  filter(SeqCount != 0) %>%
  distinct() %>%
  separate(Individual, into=as.character(c(1:2)), sep = "X")
CBwFG4 <- CBwFG4[-2]
colnames(CBwFG4)[2] = "Individual"
# Next, the counts from CBwFG4 should be merged w CBwFG3_Meta
CBwFG3_Meta <-   CBwFG3 %>%
  subset(select = -c(2:299)) %>%
  subset(select = -c(name))
CBwFG_Please <- merge(CBwFG4, CBwFG3_Meta, by=c("Isolate" ,"Individual"))
colnames(CBwFG_Please)[8] = "OldGuild"
# Next, merge Shifts with CBwFG_Please
Shifts[31,] <- c("Unassigned", "Unassigned")
CBwFG_Please2 <- merge(CBwFG_Please, Shifts) %>%
  filter(NewGuild != "DROP")

# Starting datasets for the scatterplots:
# CBwFG_Please2 : CB dataset, counts, metadata, etc
# Full_CFwFG : CF dataset, counts, metadata, etc

# Gathering OTU/iso# and taxonomy. Splice from OTU tables.
CB_Tax <- CultureBased_Taxonomy %>% # CB taxonomy
  subset(select = -c(2:8))
CF_Tax <- CultureFree_Taxonomy %>% # CF taxonomy
  subset(select = -c(2:8))
# Merge CBwFG_Please2 and Full_CFwFG with CB_Tax and CF_Tax, by OTU/iso#
CB_Tax2 <- merge(CB_Tax, CBwFG_Please2)
CF_Tax2 <- merge(CF_Tax, Full_CFwFG)

# Making overall summary tables. (incl. unresolved guild)
CB_Tax2_OverallSummary <- CB_Tax2 %>% # CB summary table
  filter(SeqCount != 0) %>%   # Drop rows with seq count = 0
  filter(NewGuild != "Unassigned") %>% # Drop all "Unassigned"
  group_by(NewGuild) %>%      # Group by the new guild for below math
  summarize(ISOpc = ((length(unique(Isolate))/595)*100),       # % of num isolates in dataset
            PEDEpc = (((length(unique(Individual)))/258)*100), # % of num millipedes in dataset
            GENUSpc = ((length(unique(Genus.join)))/168)*100,  # % of num genera (& fams in certain cases)in dataset. Includes "NA"
            `Genera (n)` = length(unique(Genus.join)))
CF_Tax2_OverallSummary <- CF_Tax2 %>% # CF summary table
  filter(SeqCount != 0) %>%  # Drop rows with seq count = 0
  filter(Guild != "Unassigned") %>% # Drop all "Unassigned"
  group_by(Guild) %>%        # Group by the new guild for below math
  summarize(SEQpc = (((sum(SeqCount))/2691541)*100), # % of num seqs in dataset
            PEDEpc = (((length(unique(SampleID)))/55) *100), # % of num millipedes in dataset
            GENUSpc = ((length(unique(Genus.join)))/482) *100, # % of num genera in dataset. Includes "NA"
            `Genera (n)` = length(unique(Genus.join)), # How many genera have this guild assignment
            OTUpc = ((length(unique(`OTU ID`)))/1008) *100) # % of num OTUs in dataset

# Generating plots:

# First, setting up a theme I can use repeatedly
GuildPlotTheme <- theme_bw() +                     # Nice set of overall edits for publication
  theme(text=element_text(size=12, color="black"), # Change font size & color overall
        axis.title.x=element_text(size=17),        # Change size of X-axis main label
        axis.text.x=element_text(size=14, color="black"),  # Change size of labels on X-axis
        axis.title.y=element_text(size=17),        # Change size of Y-axis main label
        axis.text.y=element_text(size=14, color="black"),  # Change size of labels on Y-axis
        plot.title=element_text(size=20, hjust=0.5),       # Change size & position of title text
        legend.title=element_text(size=12),        # Change size of legend title
        legend.text=element_text(size=12),         # Change size of legend label text
        panel.grid.minor = element_blank(),        # no minor grid lines
        panel.border=element_rect(color="black", size=1))  # Black panel border 

# CB scatterplot:
CB_Tax2_Scatterplot <- ggplot(CB_Tax2_OverallSummary, aes(x=PEDEpc, y=ISOpc)) +
  GuildPlotTheme +
  theme(legend.position = c(0.125, 0.8), legend.background = element_rect(color="grey")) +
  labs(x="Millipedes (%)", y="Isolates (%)") +
  geom_point(aes(fill = `Genera (n)`), size = 6, shape=21, color= "black") +
  scale_fill_gradient(low="white", high="black", space="Lab", guide="colourbar", aesthetics="fill") +
  geom_text_repel(aes(label=NewGuild, size=9), color="black", box.padding=0.7, min.segment.length=0.2, point.padding = 0.8) +
  guides(size = "none")
CB_Tax2_Scatterplot

# CF scatterplot
CF_Tax2_Scatterplot <- ggplot(CF_Tax2_OverallSummary, aes(x=PEDEpc, y=OTUpc)) +
  GuildPlotTheme +
  theme(legend.position = c(0.125, 0.8), legend.background = element_rect(color="grey")) +
  labs(x="Millipedes (%)", y="OTUs (%)") +
  geom_point(aes(fill = `Genera (n)`), size=6, shape=21, color="black") +
  scale_fill_gradient(low="white", high="black", space="Lab", guide="colourbar", aesthetics="fill") +
  geom_text_repel(aes(label=Guild, size=9), color="black", box.padding=0.7, min.segment.length=0.2, point.padding=0.8) +
  guides(size = "none")
CF_Tax2_Scatterplot

# Statistical comparisons w/ PERMANOVA & visualizations w/ NMDS

# Culture-free:

# Make an "OTU table", where instead of OTUs we have guilds.
CF_Tax3 <- CF_Tax2 %>%
  filter(SeqCount != 0) %>%   # Drop rows with seq count = 0
  filter(Guild != "Unassigned") %>% # Drop all "Unassigned"
  group_by(Guild, SampleID) %>%
  summarize(SeqCount = (sum(SeqCount))) %>%
  pivot_wider(names_from=SampleID, values_from=SeqCount)
CF_Tax3[is.na(CF_Tax3)] <- 0
collist <- CF_Tax3$Guild      # row name fix
CF_Tax3$Guild = NULL          # row name fix
rownames(CF_Tax3) <- collist  # row name fix
CF_Tax3T <- as.data.frame(t(CF_Tax3)) # flip table
braycurtis_CF_FG <- vegdist(CF_Tax3T, "bray") # distance matrix

# Generate NMDS
Bray_CF_FG_NMDS <- metaMDS(comm = braycurtis_CF_FG, distance = "bray", k = 2, trymax=999) #If you rerun this line, update the best_soln_stress!
Shepard <- stressplot(Bray_CF_FG_NMDS)
best_soln_stress3 <- 0.067 # updated 2-18-24

# Validate NMDS
BrayFG_NMDSp1 <- metaMDS(comm = OTUs.mat4, distance = "bray", k = 2, try = 500,
                       trymax = 501, engine = "monoMDS", autotransform = T)
BrayFG_NMDSp2 <- metaMDS(comm = OTUs.mat4, distance = "bray", k = 2, try = 500,
                       trymax = 501, engine = "monoMDS", autotransform = T)
BrayFG_NMDS_stress <- read.csv("BrayFG_NMDS_stress.csv", header = T)
sigmaX <- sd(BrayFG_NMDS_stress$Stress)
z.test(BrayFG_NMDS_stress$Stress, alternative = 'two.sided', mu = best_soln_stress, sigma.x = sigmaX, conf.level = .95)
# NMDS is reliable!

# Generate plot
scores <- as.data.frame(scores(Bray_CF_FG_NMDS))
scores$Site <- CultureFreeMeta$Site
scores$Sex <- CultureFreeMeta$Sex
scores$Colony <- as.factor(CultureFreeMeta$Colony)
lab <- paste0("stress = ", best_soln_stress3)
CF_FG_nmds_site <- ggplot() +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS2, shape = Site, color = Site),size = 4) + # change "colour" to show site, sex, or colony
  scale_color_manual(values=c('#009E73','#0072B2', '#D55E00'))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  annotate(geom="text", x=-0.5, y=0.13, label=lab, color="black", size=6) +
  stat_ellipse(aes(x = scores$NMDS1, y = scores$NMDS2, lty = scores$Site), show.legend = T, level = 0.95) +
  scale_shape_manual(values = c(15, 16, 17), name = "Site", breaks = c("BA", "CH", "CO"), labels = c("BA", "CH", "CO"))
CF_FG_nmds_site

# PERMANOVA by site:
betadispSite_bray <- betadisper(d = braycurtis_CF_FG, group = CultureFreeMeta$Site, bias.adjust=TRUE)
adonis2(braycurtis_CF_FG ~ CultureFreeMeta$Site, permutations = 9999)
pairwise.adonis(x = braycurtis_CF_FG, factors = CultureFreeMeta$Site, p.adjust.m = "bonferroni")

# PERMANOVA by sex:
betadispSex_bray <- betadisper(d = braycurtis_CF_FG, group = CultureFreeMeta$Sex, bias.adjust=TRUE)
adonis2(braycurtis_CF_FG ~ CultureFreeMeta$Sex, permutations = 9999)

# To compare colonies, need to slice up the data and redo vegdist.

# Start w cols 6-7 (Comparison 1 from above): drop rows where Site != BA.
CF_Tax3T_BigCols <- CF_Tax3T[1:18,]
# Cols where n=4
keepers <- c("19CO01F","19CO02F","19CO03F","19CO04F","19CO17F","19CO18F",
             "19CO19F","19CO20M","19CH10M","19CH11M","19CH13F","19CH14F",
             "19BA02M","19BA04M","19BA05M","19BA08F","19BA11F","19BA15M",
             "19BA16F","19BA20F","19CH18F","19CH30F","19CH31F","19CH32F")
CF_Tax3T_n4Cols <- CF_Tax3T %>% filter(row.names(CF_Tax3T) %in% keepers)

# Stats on big cols:
CF_Tax3T_BigCols_dist <- vegdist(CF_Tax3T_BigCols, "bray")
CFMeta_BigCols <- CultureFreeMeta %>% filter(Site == "BA")
betadispBigCols_bray <- betadisper(d = CF_Tax3T_BigCols_dist, group=CFMeta_BigCols$Colony, bias.adjust=T)
adonis2(CF_Tax3T_BigCols_dist ~ CFMeta_BigCols$Colony, permutations = 9999)

# Stats on n=4 cols (Comparison 2 from above):
CF_Tax3T_n4Cols_dist <- vegdist(CF_Tax3T_n4Cols, "bray")
CFMeta_n4Cols <- CultureFreeMeta %>% filter(SampleID %in% keepers)
betadispn4Cols_bray <- betadisper(d = CF_Tax3T_n4Cols_dist, group=CFMeta_n4Cols$Colony, bias.adjust=T)
adonis2(CF_Tax3T_n4Cols_dist ~ CFMeta_n4Cols$Colony + CFMeta_n4Cols$Colony %in% CFMeta_n4Cols$Site, permutations = 9999)
