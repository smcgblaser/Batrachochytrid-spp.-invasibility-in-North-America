#Analysis code for *INSERT MS TITLE HERE*
#Sarah McGrath-Blaser
#Summer 2022

#First, load required packages
Packages <- c("tidyverse", "ggplot2", "reshape2", "dada2", "phyloseq")
lapply(Packages, library, character.only = TRUE)

#Set working directory to project directory
setwd("C:/Users/sarah/OneDrive - University of Florida/Soil_Experiment/pathogen_growth_soils_exp")

#Useful tutorials and links
#Tutorial https://usda-ars-gbru.github.io/Microbiome-workshop/tutorials/amplicon/

################## Read in data ################################
# These files are from running the PacBio_data_analysis.R script on University of Florida's computing core, HiPerGator, and transferred locally
seqtab <- readRDS("Inputs/seqtab_final.rds")
taxtab <- readRDS("Inputs/tax_final.rds")
metadata <- read.table("Inputs/Sample_metadata_updated.csv", header=TRUE, sep = ",", stringsAsFactors=FALSE)

#Inspect distribution of sequence lengths
hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")

table(nchar(getSequences(seqtab)))

samples.out <- rownames(seqtab)
head(samples.out)

#Assigning rownames in the dataframe to SampleID https://github.com/joey711/phyloseq/issues/1020.
rownames(metadata) <- metadata$SampleID

################## Construct phyloseq object #############################
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
               sample_data(metadata),
               tax_table(taxtab))
ps

####################### Look for contamination ##############################
#Tutorial for using decontam package
#https://bioconductor.org/packages/release/bioc/vignettes/decontam/inst/doc/decontam_intro.html

#Install and load package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("decontam")
library(decontam)

ps
sample_data(ps)

head(sample_data(ps))

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_control)) + geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

head(which(contamdf.prev$contaminant))

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", 
                                 threshold=0.5,detailed=TRUE)
contamdf.prev05
table(contamdf.prev05$contaminant)


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_control == "True_sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Controls)") + ylab("Prevalence (True Samples)")

ps

ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)
ps.noncontam
sum(sample_sums(ps.noncontam))

#Number of taxa in OTU (ASV) tables before and after decontamination
36741-36648
#Removed 93 taxa that were found in DNA extraction water control, 
#PCR water control (both of which library size was 0, so I don't think there really 
#was any contamination there), Broth and Autoclaved (positive control) samples

#What were these taxa?
ps.contam <- prune_taxa(contamdf.prev05$contaminant, ps)
ps.contam
tax_talbe_ps.contam <- tax_table(ps.contam)
#Write to a csv file to more closely inspect taxa removed as 'contamination'
write.csv(tax_talbe_ps.contam,file="Outputs/tax_table_ps.contam.csv")

top5gen <- sort(tapply(taxa_sums(ps.contam), tax_table(ps.contam)[, "Genus"],
                       sum), decreasing = TRUE)[1:5]
top5gen
###################### Check Mock Community control ########################
#check the mock community. Is it what we expect to see? 
#Tutorial
#https://benjjneb.github.io/dada2/tutorial.html 
#Go to section 'Evaluate accuracy'

#Make object for mock community sequence data
mock <- seqtab["demultiplex.bc1024--bc1105.hifi_reads.fastq",]

#Drop ASVs absent in the Mock
mock <- sort(mock[mock>0], decreasing=TRUE)
cat("DADA2 inferred", length(mock), "sample sequences present in the Mock community.\n")
#Says DADA2 inferred 38 sequences from the Mock community

#refseqs <- "ZymoBIOMICS.STD.refseq.v2/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/Bacillus_subtilis_16S_170923.fasta|"
mock.ref <- getSequences("Inputs/ZymoBIOMICS.STD.refseq.v2/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/Concat_file")
mock.ref
match.ref <- sum(sapply(names(mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
#Says Of those, 15 were exact matches to the expected reference sequences
#ZymoBIOMICS Microbial Community Standard Lot# 190633
#files.zymoresearch.com%2Fprotocols%2F_d6300_zymobiomics_microbial_community_standard.pdf

#Make stacked bar plot of mock community
#Subset phyloseq object (all taxa) for mock community samples
ps_mock <- subset_samples(ps, Site=='Mock_community')
tax_table(ps_mock)
mock_comm <- ps_mock  %>% #Create new object
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #Condense by species      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.00) %>%  #Filter out low abundance counts 
  arrange(Genus) #Order species

ggplot(mock_comm, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

#Another stacked barplot of mock community data from decontaminated phyloseq object
ps_mock <- subset_samples(ps.noncontam, Site=='Mock_community')
tax_table(ps_mock)
mock_comm <- ps_mock  %>% #Create new object
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #Condense by species      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  filter(Abundance > 0.00) %>%  #Filter out low abundance counts 
  arrange(Genus) #Order species

ggplot(mock_comm, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(direction = -1) +
  scale_x_discrete(labels = NULL) +
  geom_hline(yintercept = 0, size = 1, colour = "#333333") +
  theme(legend.position = "right", 
        legend.justification = "left",
        legend.text=element_text(size=8)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))
#Removed Pseudomonas! Interesting.

#I think the phyloseq object is decontaminated and mock community looks good so can proceed with looking at core samples

######################### Summary of core samples ################
hist(sample_sums(ps.noncontam), main="Histogram: Read Counts", xlab="Total Reads", las=1, breaks=12)
#Prune ASVs that are not present in any samples
SB <- prune_taxa(taxa_sums(ps.noncontam) > 0, ps.noncontam)
sum(sample_sums(SB))
#954,440

#Subset out controls and autoclaved samples
SB <- subset_samples(SB, Soil_treatment %in% c("Non_inoc","Non_autoclaved") & Sample_or_control == "True_sample")
hist(sample_sums(SB), main="Histogram: Read Counts", xlab="Total Reads", las=1, breaks=12)
## total number of sequences per sample, sorted in order
sort(sample_sums(SB))
19084/338
# = 56.46154
19084/864
# = 22.08796
19084/338
# > 10-fold increase in difference between highest and lowest sample sequences
#Do I need to rarefy/normalize data before looking at diversity estimates? Here is some info. Not normalizing for now.
#https://www.researchgate.net/post/Should_I_rarefy_the_microbiome_count_data_before_calculating_alpha_diversity
#Info on SRS normalization: https://www.researchgate.net/publication/343392051_Improved_normalization_of_species_count_data_in_ecology_by_scaling_with_ranked_subsampling_SRS_Application_to_microbial_communities
#In /PacBio_metabarcoding saved an R file that has the script information for running SRS

top15ph <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Phylum"],sum), TRUE)[1:15]
top15g <- sort(tapply(taxa_sums(SB), tax_table(SB)[, "Genus"], sum), TRUE)[1:15]
top15ph
top15g

sum(taxa_sums(SB))
#875,171 total sequence count
#Proteobacteria = 356,644 sequences
(356644/875171)*100
#Proteobacteria make up 40.75% of all sequences
#Acidobacteria
(144554/875171)*100
#Acidobacteria make up 16.51% of all sequences
#Planctomycetota
(102681/875171)*100
#Plactomycetota make up 11.73% of all sequences
#Actinobacteriota
(88408/875171)*100
#Actinobacteriota make up 10.10% of all sequences
#Bacteroidota
(70914/875171)*100
#Bacteroidota make up 8.10% of all sequences

## 
get_taxa_unique(SB, "Kingdom")

SB_Archae = subset_taxa(SB, Kingdom == "Archaea")
#Look at what Archaea are present
get_taxa_unique(SB_Archae, "Phylum")
#Euryarchaeota

SB_Eukaryota = subset_taxa(SB, Kingdom == "Eukaryota")
get_taxa_unique(SB_Eukaryota, "Phylum")
#NA?

## 46 bacterial and archael phyla (one NA)
get_taxa_unique(SB, "Phylum")
tax_table(SB)[,"Phylum"] %>% unique %>% na.exclude %>% length
#45

#### 36,648 ASVs
SB
#This says there are 871 taxa that have a sum of 0.
sum(taxa_sums(SB) == 0)

########################## Alpha Diversity Analyses ##################################
#https://microbiome.github.io/tutorials/PlotDiversity.html
library(microbiome)
library(ggpubr)
#Prune ASVs that are not present in any samples
SB <- prune_taxa(taxa_sums(ps.noncontam) > 0, ps.noncontam)

#Subset out controls and autoclaved samples
SB <- subset_samples(SB, Soil_treatment %in% c("Non_inoc", "Non_autoclaved") & Sample_or_control == "True_sample" & Day %in% c("0","1","7","14"))

#BD_w_d0 <- subset_samples(BD_w_d0, Soil_treatment %in% c("Non_inoc","Non_autoclaved") & Sample_or_control == "True_sample" & Day %in% c("0","1","7","14"))

plot_richness(SB, x="Day", measures = c("Observed","Shannon","Simpson"))
#Get alpha diversity data for SB phyloseq object, compare across all indices
SB_alpha <- microbiome::alpha(SB, index = "all")
head(SB_alpha)
str(SB_alpha)
#Get metadata from phyloseq object
SB_meta <- meta(SB)
head(SB_meta)
SB_meta$Day
#Add diversity table to metadata
SB_meta$Observed <- SB_alpha$observed
SB_meta$Shannon <- SB_alpha$diversity_shannon
SB_meta$Simpson <- SB_alpha$diversity_gini_simpson

# create a list of pairwise comaprisons
# get the variables
SB_meta$Day <- as.factor(SB_meta$Day)
day <- levels(SB_meta$Day)
day_pairs <- combn(seq_along(day), 2, simplify = FALSE, FUN = function(i)day[i])
print(day_pairs)

SB_meta$Site <- as.factor(SB_meta$Site)
site <- levels(SB_meta$Site)
site_pairs <- combn(seq_along(site), 2, simplify = FALSE, FUN = function(i)site[i])

hist(SB_meta$Shannon)
qqnorm(SB_meta$Shannon)
shapiro.test(SB_meta$Shannon)

hist(SB_meta$Simpson)
qqnorm(SB_meta$Simpson)
shapiro.test(SB_meta$Simpson)
#Not normally distributed, will need to do use another statistical test besides ANOVA

library(dplyr)
group_by(SB_meta, Day) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(Shannon, na.rm = TRUE),
    sd = sd(Shannon, na.rm = TRUE),
    median = median(Shannon, na.rm = TRUE),
    IQR = IQR(Shannon, na.rm = TRUE)
  )

kruskal.test(Observed ~ Day, data = SB_meta)
kruskal.test(Shannon ~ Day, data = SB_meta)
kruskal.test(Simpson ~ Day, data = SB_meta)
#No significant difference by day

#Now testing by site
kruskal.test(Observed ~ Site, data = SB_meta)
kruskal.test(Simpson ~ Site, data = SB_meta)
#Both show significant differences, but where? Need to do Dunn's test.

#install.packages("FSA")
library(FSA)
dunnTest(Observed ~ Site, data = SB_meta, method = "bh")
dunnTest(Simpson ~ Site, data = SB_meta, method = "bh")


SB_meta_Bd <- subset(SB_meta, Pathogen_type =="Bd")
hist(SB_meta_Bd$Shannon)
qqnorm(SB_meta_Bd$Shannon)
shapiro.test(SB_meta_Bd$Shannon)

hist(SB_meta_Bd$Simpson)
qqnorm(SB_meta_Bd$Simpson)
shapiro.test(SB_meta_Bd$Simpson)

kruskal.test(Observed ~ Day, data = SB_meta_Bd)
kruskal.test(Simpson ~ Day, data = SB_meta_Bd)
kruskal.test(Observed ~ Site, data = SB_meta_Bd)
kruskal.test(Simpson ~ Site, data = SB_meta_Bd)
#All p > 0.05

SB_meta_Bsal <- subset(SB_meta, Pathogen_type =="Bsal")
hist(SB_meta_Bsal$Shannon)
qqnorm(SB_meta_Bsal$Shannon)
shapiro.test(SB_meta_Bsal$Shannon)

hist(SB_meta_Bsal$Simpson)
qqnorm(SB_meta_Bsal$Simpson)
shapiro.test(SB_meta_Bsal$Simpson)

kruskal.test(Observed ~ Day, data = SB_meta_Bsal)
kruskal.test(Simpson ~ Day, data = SB_meta_Bsal)
kruskal.test(Observed ~ Site, data = SB_meta_Bsal)
kruskal.test(Simpson ~ Site, data = SB_meta_Bsal)
#Bsal Simpson by site is significant. Run Dunn's test.
dunnTest(Simpson ~ Site, data = SB_meta_Bsal, method = "bh")

####################### Summary stats Bd and Bsal #####################################
SB_Bd <- subset_samples(SB, Pathogen_type == "Bd" & Day == "14")
SB_Bd@sam_data
top15ph_Bd <- sort(tapply(taxa_sums(SB_Bd), tax_table(SB_Bd)[, "Phylum"],sum), TRUE)[1:15]
top15ph_Bd
sum(taxa_sums(SB_Bd))
#161233 total sequence count
#Proteobacteria = 63489 sequences
(63489/161233)*100
#39.37%
#Acidobacteriota = 27984
(27984/161233)*100
#17.35%
#Planctomycetota = 18559
(18559/161233)*100
#11.5%
#Actinobacteriota = 16045
(16045/161233)*100
#9.95%
#Bacteroidota = 12724
(12724/161233)*100
#7.89%
SB_Bsal <- subset_samples(SB, Pathogen_type == "Bsal" & Day == "14")
SB_Bsal@sam_data
top15ph_Bsal <- sort(tapply(taxa_sums(SB_Bsal), tax_table(SB_Bsal)[, "Phylum"],sum), TRUE)[1:15]
top15ph_Bsal
sum(taxa_sums(SB_Bsal))
#167853 total sequence count
#Proteobacteria = 64752 sequences
(64752/167853)*100
#38.57%
#Acidobacteriota = 28164
(28164/167853)*100
#16.77%
#Planctomycetota = 19460
(19460/167853)*100
#11.59%
#Actinobacteriota = 17789
(17789/167853)*100
#10.59%
#Bacteroidota = 13840 
(13840/167853)*100
#8.24%



################################ Alpha diversity non-inoculated soils ###########################
#Subset to just the non-inoculated samples
NonInoc <- subset_samples(ps.noncontam, Soil_treatment == "Non_inoc")

## total number of sequences per sample, sorted in order
sort(sample_sums(NonInoc))
13499/10779

top15ph_NonInoc <- sort(tapply(taxa_sums(NonInoc), tax_table(NonInoc)[, "Phylum"],sum), TRUE)[1:15]
top15ph_NonInoc

sum(taxa_sums(NonInoc))
#59,586 total sequence count
#Proteobacteria = 23660 sequences
(23660/59586)*100
#Proteobacteria is 39.7% of Non Inoc seqs
#Acidobacteriota = 9592 sequences
(9592/59586)*100
#Acidobacteriota is 16% of Non Inoc seqs
#Planctomycetota = 6743 sequences
(6743/59586)*100
#Planctomycetota is 11.3% of Non Inoc seqs
#Bacteroidota = 6561 sequences
(6561/59586)*100
#Bacteroidota is 11% of Non Inoc seqs
#Actinobacteriota = 4939 sequences
(4939/59586)*100
#Actinobacteriota is 8.2% of Non Inoc seqs

########################################### RDA Analysis + PCoA ###############################
SB@sam_data
SB@sam_data$Elevation[grepl("Fontana_Dam", SB@sam_data$Site)] <- "656"
SB@sam_data$Elevation[grepl("Newfound_Gap", SB@sam_data$Site)] <- "1577"
SB@sam_data$Elevation[grepl("Clingmans_Dome", SB@sam_data$Site)] <- "2027"
SB@sam_data$Elevation[grepl("Mount_Cammerer", SB@sam_data$Site)] <- "1278"
SB@sam_data$Elevation[grepl("Davenport_Gap", SB@sam_data$Site)] <- "744"
#Do I want non-inoc in this calculation or not??
SB@sam_data$Soil_name[grepl("Fontana_Dam", SB@sam_data$Site)] <- "Mesic_soft_metasandstone_soco_stecoah_soils"
SB@sam_data$Soil_name[grepl("Newfound_Gap", SB@sam_data$Site)] <- "Frigid_Anakeesta_Slate_Luftee_Anakeesta_soils"
SB@sam_data$Soil_name[grepl("Clingmans_Dome", SB@sam_data$Site)] <- "Frigid_hard_sandstone_breakneck_pullback_soils"
SB@sam_data$Soil_name[grepl("Mount_Cammerer", SB@sam_data$Site)] <- "Mesic_soft_metasandstone_soco_stecoah_soils"
SB@sam_data$Soil_name[grepl("Davenport_Gap", SB@sam_data$Site)] <- "Mesic_soft_metasandstone_soco_stecoah_soils"

asv_tab <- otu_table(SB)
if (taxa_are_rows(asv_tab)) {
  asv_tab <- t(asv_tab)
}
asv_mat <- as(asv_tab, "matrix")
asv_df <- as.data.frame(asv_mat)
SB_sam_data <- as.data.frame(sample_data(SB))
SB_df_RDA <- cbind(SB_sam_data, asv_df)
dim(SB_df_RDA)
rownames(SB_df_RDA)<-NULL
SB_df_RDA <- subset (SB_df_RDA, select = -c(6:9))


fix(SB_df_RDA)
#row.names(SB_df_RDA) = SB_df_RDA$Sample_Name
species = SB_df_RDA[,8:36]
environment = SB_df_RDA[,3:7]
fix(species)
fix(environment)
species001= (species + 0.001)
#SB_df_RDA$Day[SB_df_RDA$Day == 0] <- 0.5
fix(species001)

library(vegan)
dbRDA = rda(species001 ~ Day+Pathogen_type+Site+Elevation+Soil_name, environment, dist="bray", na.action = na.exclude)
anova(dbRDA, by="terms", permu=200)
#Site is significant
dbRDA = capscale(species001 ~ Day+Pathogen_type+Elevation+Soil_name, environment, dist="bray", na.action = na.exclude)
anova(dbRDA, by="terms", permu=200)
#Elevation is significant
dbRDA = capscale(species001 ~ Day+Pathogen_type+Soil_name, environment, dist="bray", na.action = na.exclude)
anova(dbRDA, by="terms", permu=200)
#Soil name is significant
dbRDA = capscale(species001 ~ Day+Pathogen_type, environment, dist="bray", na.action = na.exclude)
anova(dbRDA, by="terms", permu=200)
#Neither Day nor Pathogen type are significant
plot(dbRDA) # use base plot, might be done with ggplot2
#devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
autoplot(dbRDA, arrows = T, legend = "none")

#################################################################
otu_ps <- as(otu_table(ps), "matrix")
#if(taxa_are_rows(ps)){otu_ps <- t(otu_ps)}
# Coerce to data.frame
OTUdf <- as.data.frame(otu_ps)
OTUdf <- t(OTUdf)
Index <- (1:36741)
Index
OTUdf <- cbind(Index,OTUdf)
tax_table_ps <- tax_table(ps)
taxa <- read.csv("taxonomy.csv",header = T)
taxa <- as.data.frame(taxa)
taxa
OTUdf <- cbind(taxa, OTUdf)
OTUdf <- OTUdf[,c(2:107,1)]
rownames(OTUdf)<-NULL
sample_info<-metadata$Sample_Name
sample_info <- as.list(sample_info)
colnames(OTUdf)[2:106]<-sample_info

write_tsv(OTUdf, file = "OTUdf.tsv")

otus_ps <- otu_table(ps)
write.csv(otus_ps,file = 'c:/Users/sarah/OneDrive - University of Florida/Soil_Experiment/PacBio_metabarcoding/otu_table_ps_1.csv')

write.csv(tax_table_ps,file = 'tax_table_ps.csv')

############################ Beta Diversity Analyses ##################################
library(phyloseq)
library(vegan)
library(RVAideMemoire)

#Going to use Bray-Curtis for beta diversity estimates since we are interested in the presence/absence of taxa as well as the abundance of taxa between samples

#Make distance matrix
bray_SB <- distance(SB, "bray")
#Make dataframe for vegan to use for analysis
df_SB <- as(sample_data(SB), "data.frame")

#Set permutations
perm <- how(nperm = 999)
#Set blocks for stratification
setBlocks(perm) <- with(df_SB, Site)
#Run permutaion test
adonis2(bray_SB ~ Day, data = df_SB, permutations = perm)
#No significant differences between days
adonis2(bray_SB ~ Day/Site, data = df_SB, permutations = perm)

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",force = T)
library(pairwiseAdonis)
pairwise.adonis2(bray_SB ~ Day/Site, data = df_SB, strata = 'Site')
#Differences between 1&7 (p = 0.045) and 14&7 (p = 0.031)

#Bd only
SB_Bd <- subset_samples(SB, Pathogen_type == "Bd")
bray_SB_Bd <- distance(SB_Bd, "bray")
df_SB_Bd <- as(sample_data(SB_Bd), "data.frame")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df_SB_Bd, Site)
adonis2(bray_SB_Bd ~ Day, data = df_SB_Bd, permutations = perm)
pairwise.adonis2(bray_SB_Bd ~ Day/Site, data = df_SB_Bd, strata = 'Site')
#No significant differences
pairwise.adonis2(bray_SB_Bd ~ Site, data = df_SB_Bd)
#Except when comparing by site, all sites are significantly different from one another

#Bsal only
SB_Bsal <- subset_samples(SB, Pathogen_type == "Bsal")
bray_SB_Bsal <- distance(SB_Bsal, "bray")
df_SB_Bsal <- as(sample_data(SB_Bsal), "data.frame")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df_SB_Bsal, Site)
adonis2(bray_SB_Bsal ~ Day, data = df_SB_Bsal, permutations = perm)
pairwise.adonis2(bray_SB_Bsal ~ Day/Site, data = df_SB_Bsal, strata = 'Site')
#Significant difference between days 1&14 (p = 0.004)
pairwise.adonis2(bray_SB_Bsal ~ Site, data = df_SB_Bsal)

#SB <- subset_samples(SB, Soil_treatment %in% c("Non_autoclaved") & Sample_or_control == "True_sample" & Day %in% c("1","14"))
SB <- subset_samples(SB, Soil_treatment %in% c("Non_autoclaved") & Sample_or_control == "True_sample" & Day == "14")
SB@sam_data
meta <- as(sample_data(SB), "data.frame")
meta$PathDay <- paste(meta$Pathogen_type, meta$Day)
#Add pathogen type and day to make new variable to test variance
bray_SB <- distance(SB, "bray")
mod<-betadisper(bray_SB, meta$PathDay)
permutest(mod)
plot(mod)

SB <- subset_samples(SB, Soil_treatment %in% c("Non_autoclaved") & Sample_or_control == "True_sample" & Day == "14")
meta <- as(sample_data(SB), "data.frame")
bray_SB <- distance(SB, "bray")
mod<-betadisper(bray_SB, meta$Pathogen_type)
permutest(mod)
plot(mod)

SB <- subset_samples(SB, Soil_treatment %in% c("Non_autoclaved") & Sample_or_control == "True_sample" & Day == "1")
meta <- as(sample_data(SB), "data.frame")
bray_SB <- distance(SB, "bray")
mod<-betadisper(bray_SB, meta$Pathogen_type)
permutest(mod)
plot(mod)

SB <- subset_samples(SB, Soil_treatment %in% c("Non_autoclaved") & Sample_or_control == "True_sample")
SB_Bd <- subset_samples(SB, Pathogen_type == "Bd")
bray_SB_Bd <- distance(SB_Bd, "bray")
Bd_meta <- as(sample_data(SB_Bd), "data.frame")
Bd_meta
mod<-betadisper(bray_SB_Bd, Bd_meta$Day)
permutest(mod)
first_Bd <- plot(mod)
first_Bd

SB_Bd <- subset_samples(SB_Bd, Day == "1")
bray_SB_Bd <- distance(SB_Bd, "bray")
Bd_meta <- as(sample_data(SB_Bd), "data.frame")
mod<-betadisper(bray_SB_Bd, Bd_meta$Site)
permutest(mod)
first_Bd <- plot(mod)
first_Bd

SB_Bd <- subset_samples(SB, Pathogen_type == "Bd")
SB_Bd <- subset_samples(SB_Bd, Day == "14")
bray_SB_Bd <- distance(SB_Bd, "bray")
Bd_meta <- as(sample_data(SB_Bd), "data.frame")
mod<-betadisper(bray_SB_Bd, Bd_meta$Site)
permutest(mod)
second_Bd <- plot(mod)
second_Bd

#install.packages("patchwork")
library(patchwork)
all_Bd <- (first_Bd + second_Bd)

Bd_meta <- as(sample_data(SB_Bd), "data.frame")
mod<-betadisper(bray_SB_Bd, Bd_meta$Site)
permutest(mod)
plot(mod)

SB_Bsal <- subset_samples(SB, Pathogen_type = "Bsal")
SB_Bsal <- subset_samples(SB_Bsal, Day == "14")
bray_SB_Bsal <- distance(SB_Bsal, "bray")
Bsal_meta <- as(sample_data(SB_Bsal), "data.frame")
mod<-betadisper(bray_SB_Bsal, Bsal_meta$Site)
permutest(mod)
plot(mod)