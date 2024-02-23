## Visualization comparing raw and rarefied datasets

#Run Analyses.R script file to get the SB phyloseq object with the raw sequences (can highlight lines 1-181 and hit run)
SB

#Copy code from section titled 'rarefied dataset'
#Rename phyloseq object to SB_rarefied to keep objects separate
SB_rarefied <- rarefy_even_depth(SB, rngseed=1, sample.size=0.9*min(sample_sums(SB)), replace=F)
SB_rarefied

####################### Alpha diversity ########################
plot_richness(SB_rarefied, x="Site", measures = c("Observed","Simpson")) +
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

ggsave("Figures/Supp_fig_rarefied.eps",dpi = 300, width = 5, height = 4)


##Re-run code for SB phyloseq object to get non-rarefied dataset plot
plot_richness(SB, x="Site", measures = c("Observed","Simpson")) +
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

ggsave("Figures/Supp_fig_non-rarefied.eps",dpi = 300, width = 5, height = 4)


############################ Beta diversity #####################
#Bray-Curtis PCoA plots
#Non-rarefied
SB_bray <- ordinate(SB, "PCoA", "bray")
p1 = plot_ordination(SB, SB_bray, type="samples", color="Site")
p1 + theme_bw()

ggsave("Figures/Supp_fig_bray_non-rarefied.eps",dpi = 300)

#Rarefied
SB_rarefied_bray <- ordinate(SB_rarefied, "PCoA", "bray")
p2 = plot_ordination(SB_rarefied, SB_rarefied_bray, type="samples", color="Site")
p2 + theme_bw()

ggsave("Figures/Supp_fig_bray_rarefied.eps",dpi = 300)