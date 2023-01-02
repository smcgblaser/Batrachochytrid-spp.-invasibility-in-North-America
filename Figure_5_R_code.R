#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("microbiomeMarker")
library(microbiomeMarker)

set.seed(4444)

#Run Analyses.R script until you get phyloseq object ps.noncontam
#Get phyloseq object ready for LDA analysis
#SB stands for Soil Bacteria

SB_for_lda <- subset_samples(ps.noncontam, Soil_treatment %in% c("Non_inoc","Non_autoclaved") & Sample_or_control == "True_sample" & Day %in% c("1","14"))
SB_for_lda@sam_data #Check to make sure all samples are non-autocalved, no controls, and just days 1 and 14 for Bd and Bsal
SB_for_lda <- prune_taxa(taxa_sums(SB_for_lda) > 0, SB_for_lda) #Get rid of any ASVs that contain 0's for every sample
SB_for_lda

#Used this tutorial for changing the taxa names to something other than the hierarchical taxonomy so we can see exactly what and where those seqs are coming from
#https://david-barnett.github.io/microViz/reference/tax_name.html#:~:text=Simple%20way%20to%20set%20unique%20taxa_names%20for%20phyloseq,%28pick%20an%20appropriate%20prefix%20or%20set%20your%20own%29
library(microViz)

#Save original taxa names
old_taxa_names <- taxa_names(SB_for_lda)


SB_for_lda <- tax_name(SB_for_lda, prefix = "ASV", pad_number = FALSE, rank = "Family")
head(taxa_names(SB_for_lda))

#Store new names with old names in dataframe for reference
names_df <- tibble::tibble(old = old_taxa_names, new = taxa_names(SB_for_lda))
taxonomy_table <- SB_for_lda@tax_table
write.csv(taxonomy_table, "taxonomy_table.csv")

mm_lefse <- run_lefse(
  SB_for_lda,
  taxa_rank = "none",
  wilcoxon_cutoff = 0.05,
  group = "Day",
  kw_cutoff = 0.05,
  norm = "CPM", #This is the same normalization that the Huttenhower Galaxy portal uses so I chose it
  multigrp_strat = FALSE,
  lda_cutoff = 2
)

mm_table <- mm_lefse@marker_table
write.csv(mm_table, "mm_table.csv")
getwd()

features <- mm_table$feature
features

mm <- tax_select(SB_for_lda, tax_list = features)
mm@tax_table
mm@otu_table
mm_lefse@sam_data

library(ggplot2)
plot_bar(mm, "Order", fill="Pathogen_type", facet_grid=Day~Site)
plot_heatmap(mm_lefse, group = "Day", transform = "log10p", sample_label = T)

#Both Days
bubble_marker  <- mm %>% #Create new object
  tax_glom(taxrank = "Order", NArm=FALSE) %>% #Condense by Genus      
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  #unite("PathSite", Pathogen_type:Site, remove = FALSE) %>%
  group_by(Pathogen_type,Site,Day,Order,OTU) %>%
  summarize(Avg_abund = mean(Abundance,na.rm=T)) %>%
  ungroup() %>%
  write.csv("bubble_marker.csv")
  #mutate(Avg_abund = ifelse(Day == "14",
  #                               Avg_abund,
  #                               -1*Avg_abund))

bubble_marker1  <- mm %>% #Create new object
  tax_glom(taxrank = "Order", NArm=FALSE) %>% #Condense by Genus      
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  #unite("PathSite", Pathogen_type:Site, remove = FALSE) %>%
  group_by(Pathogen_type,Site,Day,OTU) %>%
  summarize(Avg_abund = mean(Abundance,na.rm=T)) %>%
  ungroup() %>%
  spread(key = Day, value = Avg_abund) %>%
  mutate(Diff = bubble_marker1$`14`- bubble_marker1$`1`)

#Change levels of site so facets are by elevation
bubble_marker <- as.data.frame(bubble_marker)
bubble_marker$Site <- as.factor(bubble_marker$Site)
levels(bubble_marker$Site)
bubble_marker$Site <- factor(bubble_marker$Site, levels = c("Fontana_Dam","Davenport_Gap","Mount_Cammerer","Newfound_Gap","Clingmans_Dome"))

bubble_marker1 <- as.data.frame(bubble_marker1)
bubble_marker1$Site <- as.factor(bubble_marker1$Site)
levels(bubble_marker1$Site)
bubble_marker1$Site <- factor(bubble_marker1$Site, levels = c("Fontana_Dam","Davenport_Gap","Mount_Cammerer","Newfound_Gap","Clingmans_Dome"))

#dev.off()
ggplot(bubble_marker[which(bubble_marker$Avg_abund>0),], aes(x = factor(Day), y = OTU, fill = Pathogen_type, size = Avg_abund)) + 
  geom_point(width = 0.2, height = 0, alpha=0.5, shape=21, color="black") +
  scale_fill_manual(values = c("steelblue","lightseagreen")) +
  scale_size(range = c(.1, 20), name="Abundance") +
  scale_x_discrete(breaks = c(1,14)) +
  facet_grid(~Site*Pathogen_type,scales = "free_x") +
  #scale_size(name = "Abundance") +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  guides(fill="none")

breaks_values <- pretty(bubble_marker$Avg_abund)
breaks_values <- c(-50, -25, 0, 25, 50)
breaks_values

ggplot(bubble_marker, aes(x = OTU, y = Avg_abund, fill = Day)) + 
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values),
                     limits = c(-50,50)) +
  geom_hline(yintercept = 0) +
  #scale_fill_discrete(values = c("steelblue","lightseagreen")) +
  #scale_size(range = c(.1, 20), name="Abundance") +
  #scale_x_discrete(breaks = c(1,14)) +
  facet_grid(~Site*Pathogen_type,scales = "free_x") +
  theme(panel.grid.major=element_line(linetype=1,color="lightgrey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank(),
        text = element_text(size = 12)) +
  guides(fill="none")

#install.packages("Cairo")
library(Cairo)

cairo_ps("Figure_5.eps", width = 15, height = 10)

ggplot(bubble_marker1, aes(x = OTU, y = Diff, fill = Pathogen_type)) + 
  geom_bar(stat = "identity",alpha = 0.5,color = "black") +
  coord_flip() +
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values),
                     limits = c(-50,50)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("cyan4","seagreen3")) +
  facet_grid(~Site*Pathogen_type,scales = "free_x") +
  theme(panel.grid.major=element_line(linetype=1,color="lightgrey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank(),
        text = element_text(size = 12)) +
  guides(fill="none")
dev.off()


Day1 <- subset_samples(mm, Day == "1")
taxD1 <- Day1@tax_table
write.csv(taxD1, "taxD1.csv")

bubble_marker  <- Day1 %>% #Create new object
  tax_glom(taxrank = "Order", NArm=FALSE) %>% #Condense by Genus      
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() #%>% #Reformat to long form instead of wide
  #filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
  #arrange(Genus) #Order Genus

bubble_marker <- as.data.frame(bubble_marker)
bubble_marker$Site <- as.factor(bubble_marker$Site)
levels(bubble_marker$Site)
bubble_marker$Site <- factor(bubble_marker$Site, levels = c("Fontana_Dam","Davenport_Gap","Mount_Cammerer","Newfound_Gap","Clingmans_Dome"))


ggplot(bubble_marker[which(bubble_marker$Abundance>0),], aes(x = Sample_Name, y = Order, fill = Pathogen_type, size = Abundance)) + 
  geom_point(alpha=0.5, shape=21, color="black") +
  #scale_size(range = c(.1, 24), name="Abundance") +
  facet_grid(Day~Site,scales = "free_x") +
  #scale_size(name = "Abundance") +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  guides(fill="none")

ggsave("Day1_lefse_results.jpg")

Day14 <- subset_samples(mm, Day == "14")
taxD14 <- Day14@tax_table
write.csv(taxD14, "taxD14.csv")

#oldData <- Day14 %>%
#  psmelt() %>%
#  unite("PathSite", Pathogen_type:Site, remove = FALSE)
  
bubble_marker_2  <- Day14 %>% #Create new object
  tax_glom(taxrank = "Order", NArm=FALSE) %>% #Condense by Genus      
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% #Transform counts
  psmelt() %>% #Reformat to long form instead of wide
  unite("PathSite", Pathogen_type:Site, remove = FALSE) %>%
  group_by(Pathogen_type,Site,Day,Order,OTU) %>%
  summarize(Avg_abund = mean(Abundance,na.rm=T)) %>%
  ungroup() #%>%
 
#newData <- right_join(oldData, bubble_marker_2, by = "PathSite")

#filter(Abundance > 0.02) %>%  #Filter out low abundance counts        
#arrange(Genus) #Order Genus

bubble_marker_2 <- as.data.frame(bubble_marker_2)
bubble_marker_2$Site <- as.factor(bubble_marker_2$Site)
levels(bubble_marker_2$Site)
bubble_marker_2$Site <- factor(bubble_marker_2$Site, levels = c("Fontana_Dam","Davenport_Gap","Mount_Cammerer","Newfound_Gap","Clingmans_Dome"))
#bubble_marker_2$Average <- paste(bubble_marker_2$Pathogen_type, bubble_marker_2$Site)


#dev.off()
ggplot(bubble_marker_2[which(bubble_marker_2$Avg_abund>0),], aes(x = Pathogen_type, y = OTU, fill = Pathogen_type, size = Avg_abund)) + 
  geom_point(alpha=0.5, shape=21, color="black") +
  scale_fill_manual(values = c("steelblue","lightseagreen")) +
  scale_size(range = c(.1, 24), name="Abundance") +
  facet_grid(Day~Site,scales = "free_x") +
  #scale_size(name = "Abundance") +
  theme(panel.grid.major=element_line(linetype=1,color="grey"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0),
        panel.background = element_blank()) +
  guides(fill="none")

ggsave("Day14_lefse_results_1.jpg")
