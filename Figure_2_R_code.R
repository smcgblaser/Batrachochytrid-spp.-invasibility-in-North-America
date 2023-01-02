#### Stacked barplots for non-inoculated samples ####

#Check sample metadata
NonInoc@sam_data
table(tax_table(NonInoc)[, "Phylum"], exclude = NULL)
NonInoc <- subset_taxa(NonInoc, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
library(microViz)
#NonInoc %>% tax_fix() %>% tax_agg(rank = "Phylum", add_unique = F)

#install.packages('Cairo')
library(Cairo)

#tiff("Figures/Non_inoculated_soils_barplot.tiff", units = "mm", width = 120, height = 120, res = 300)
#setEPS()                                             # Set postscript arguments
#postscript("Figures/Non_inoculated_soils_barplot.eps", width = 120, height = 120)                           # Start graphics device driver
NonInoc %>%
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    Site = factor(
      Site, levels = c("Fontana_Dam", "Davenport_Gap", "Mount_Cammerer","Newfound_Gap","Clingmans_Dome"), labels = c("FD","DG","MC","NG","CD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 15,
    bar_outline_colour = NA, facet_by = "Site"
  ) + labs(x = NULL, y = NULL) +
  theme_classic() +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) 
#dev.off()

############### Stacked barplots x pathogen x day for Figure 2 ########################
library(microViz)
#NonInoc %>% tax_fix() %>% tax_agg(rank = "Phylum", add_unique = F)
SB_Bd_taxa_plot <- subset_samples(SB, Pathogen_type == "Bd" & Day == "14")
#Check sample metadata
SB_Bd_taxa_plot@sam_data
table(tax_table(SB_Bd_taxa_plot)[, "Phylum"], exclude = NULL)
SB_Bd_taxa_plot <- subset_taxa(SB_Bd_taxa_plot, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#SB_Bd_taxa_plot <- merge_samples(SB_Bd_taxa_plot, "Site")
SB_Bd_taxa_plot@sam_data
#tiff("Figures/Bd_Day14_barplots.tiff", units = "mm", width = 120, height = 120, res = 300)
SB_Bd_taxa_plot %>%
  #ps_select(Sample_Name, Site) %>% 
  #phyloseq::merge_samples(group = "Site") %>%
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    Site = factor(
      Site, levels = c("Fontana_Dam", "Davenport_Gap", "Mount_Cammerer","Newfound_Gap","Clingmans_Dome"), labels = c("FD","DG","MC","NG","CD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 16,
    bar_outline_colour = NA, facet_by = "Site",
    palette = distinct_palette(16, pal = "kelly")
  ) + labs(x = NULL, y = NULL) +
  theme_classic() +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
  ) 
dev.off()

#NonInoc %>% tax_fix() %>% tax_agg(rank = "Phylum", add_unique = F)
SB_Bsal_taxa_plot <- subset_samples(SB, Pathogen_type == "Bsal" & Day == "14")
#Check sample metadata
SB_Bsal_taxa_plot@sam_data
table(tax_table(SB_Bsal_taxa_plot)[, "Phylum"], exclude = NULL)
SB_Bsal_taxa_plot <- subset_taxa(SB_Bsal_taxa_plot, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#SB_Bd_taxa_plot <- merge_samples(SB_Bd_taxa_plot, "Site")
SB_Bsal_taxa_plot@sam_data
#tiff("Figures/Bsal_Day14_barplots.tiff", units = "mm", width = 120, height = 120, res = 300)
SB_Bsal_taxa_plot %>%
  #ps_select(Sample_Name, Site) %>% 
  #phyloseq::merge_samples(group = "Site") %>%
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    Site = factor(
      Site, levels = c("Fontana_Dam", "Davenport_Gap", "Mount_Cammerer","Newfound_Gap","Clingmans_Dome"), labels = c("FD","DG","MC","NG","CD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 16,
    bar_outline_colour = NA, facet_by = "Site",
    palette = distinct_palette(16, pal = "kelly")
  ) + labs(x = NULL, y = NULL) +
  theme_classic() +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) 
dev.off()

library(microViz)
#NonInoc %>% tax_fix() %>% tax_agg(rank = "Phylum", add_unique = F)
SB_taxa_plot <- subset_samples(SB, Day == "14")
#Check sample metadata
SB_taxa_plot@sam_data
table(tax_table(SB_taxa_plot)[, "Phylum"], exclude = NULL)
SB_taxa_plot <- subset_taxa(SB_taxa_plot, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#SB_Bd_taxa_plot <- merge_samples(SB_Bd_taxa_plot, "Site")
SB_taxa_plot@sam_data
#tiff("Figures/Day14_barplots.tiff", units = "mm", width = 120, height = 120, res = 300)
SB_taxa_plot %>%
  #ps_select(Sample_Name, Site) %>% 
  #phyloseq::merge_samples(group = "Site") %>%
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    Site = factor(
      Site, levels = c("Fontana_Dam", "Davenport_Gap", "Mount_Cammerer","Newfound_Gap","Clingmans_Dome"), labels = c("FD","DG","MC","NG","CD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 15,
    bar_outline_colour = NA, palette = distinct_palette(16, pal = "kelly")
  ) +
  facet_wrap(Pathogen_type ~ Site, 
             scales = "free" # these options are critically important!
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
  ) 
dev.off()

#tiff("Figures/Day14_barplots_no_labels.tiff", units = "mm", width = 120, height = 120, res = 300)
SB_taxa_plot %>%
  #ps_select(Sample_Name, Site) %>% 
  #phyloseq::merge_samples(group = "Site") %>%
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    Site = factor(
      Site, levels = c("Fontana_Dam", "Davenport_Gap", "Mount_Cammerer","Newfound_Gap","Clingmans_Dome"), labels = c("FD","DG","MC","NG","CD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 15,
    bar_outline_colour = NA, palette = distinct_palette(16, pal = "kelly")
  ) +
  facet_wrap(Pathogen_type ~ Site, 
             scales = "free" # these options are critically important!
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        strip.text.x = element_blank()
  ) 
dev.off()

#### Boxplots ####
SB_meta_Bd <- subset(SB_meta, Pathogen_type =="Bd")

levels(SB_meta_Bd$Site)
levels(SB_meta_Bd$Site) <- c("CD","DG","FD","MC","NG")

#tiff("Figures/Bd_alpha_boxplots_Simpson.tiff", units = "mm", width = 200, height = 120, res = 300)
ggplot(SB_meta_Bd, aes(x = Day, y = Simpson)) +
  geom_boxplot() +
  facet_wrap(~Site, scales = "free") +
  scale_y_continuous (name = NULL, limits = c (0.984,1)) +
  xlab(NULL) +
  theme_classic() +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())
dev.off()

SB_meta_Bsal <- subset(SB_meta, Pathogen_type =="Bsal")

levels(SB_meta_Bsal$Site)
levels(SB_meta_Bsal$Site) <- c("CD","DG","FD","MC","NG")

#tiff("Figures/Bsal_alpha_boxplots_Simpson.tiff", units = "mm", width = 200, height = 120, res = 300)
ggplot(SB_meta_Bsal, aes(x = Day, y = Simpson)) +
  geom_boxplot() +
  facet_wrap(~Site, scales = "free") +
  scale_y_continuous (name = NULL, limits = c (0.984,1)) +
  xlab(NULL) +
  theme_classic() +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())
dev.off()

#### Alpha diversity plots for soils before pathogen introduction, Observed & Simpson ####

library(microbiome)
library(ggpubr)

#tiff("Figures/Non_inoculated_soils_alpha_metrics.tiff", units = "mm", width = 120, height = 120, res = 300)

plot_richness(NonInoc, x="Site", measures = c("Observed","Simpson")) +
  scale_x_discrete(
    "Site",
    limits = c("Fontana_Dam","Davenport_Gap","Newfound_Gap","Mount_Cammerer","Clingmans_Dome"),
    labels = c(
      "Fontana_Dam" = "FD",
      "Davenport_Gap" = "DG",
      "Newfound_Gap" = "NG",
      "Mount_Cammerer" = "MC",
      "Clingmans_Dome" = "CD"
    )) +
  theme_classic() +
  theme_bw() +
  theme(axis.title.x=element_blank()) +
  theme(text=element_text(size=20))

ggsave("Figures/Fig2_partB.eps",dpi = 300)



