#################################################################################

#Make another phyloseq object but include day 0 samples
#Object shorthand -> Beta Diversity with day 0
BD_w_d0 <- prune_taxa(taxa_sums(ps.noncontam) > 0, ps.noncontam)
#Subset out controls and autoclaved samples but keep day 0 samples
BD_w_d0 <- subset_samples(BD_w_d0, Soil_treatment %in% c("Non_inoc","Non_autoclaved") & Sample_or_control == "True_sample" & Day %in% c("0","1","7","14"))
sample_data(BD_w_d0)

#ASV table 
asv_tab <- otu_table(BD_w_d0)
if (taxa_are_rows(asv_tab)) {
  asv_tab <- t(asv_tab)
}

#Make it a matrix
asv_mat <- as(asv_tab, "matrix")
asv_df <- as.data.frame(asv_mat)
BD_sam_data <- as.data.frame(sample_data(BD_w_d0))
BD_w_d0_df <- cbind(BD_sam_data, asv_df)
dim(BD_w_d0_df)
rownames(BD_w_d0_df)<-NULL
#Change column names from sequences (too long) to ASV
names(BD_w_d0_df)[10:36657] <- paste("ASV", 1:36648, sep="")

#Need to duplicate Non_inoc rows and make separate rows for each pathogen type (Bd and Bsa) so that non-inoc/day 0 data is in both Bd and Bsal plots
#Make new df and duplicate each row of the first 5 rows (non-inoc samples)
#BD_w_d0_df_1 <- BD_w_d0_df[rep(seq_len(5), each = 2), ]
#Now change pathogen type from NA to Bd and Bsal respectively
#There has literally got to be a one-line way to do this, I just can't figure out what it is...so here we are
#BD_w_d0_df_1[1,4] = "Bd"
#BD_w_d0_df_1[2,4] = "Bsal"
#BD_w_d0_df_1[3,4] = "Bd"
#BD_w_d0_df_1[4,4] = "Bsal"
#BD_w_d0_df_1[5,4] = "Bd"
#BD_w_d0_df_1[6,4] = "Bsal"
#BD_w_d0_df_1[7,4] = "Bd"
#BD_w_d0_df_1[8,4] = "Bsal"
#BD_w_d0_df_1[9,4] = "Bd"
#BD_w_d0_df_1[10,4] = "Bsal"

#Remove first 5 rows of original df
#BD_w_d0_df <- BD_w_d0_df[-c(1:5), ]

#Use rbind function to add both dataframes back together
#BD_w_d0_df <- rbind(BD_w_d0_df_1,BD_w_d0_df)

#Now need to make Bray-Curtis dissimilarity matrix
#Get df into correct format
bcdf <- BD_w_d0_df %>%
  select(Sample_Name, starts_with("ASV")) 

bcdf %>%
  as_tibble() %>%
  pivot_longer(-Sample_Name) %>%
  group_by(Sample_Name) %>%
  summarize(n_seqs = sum(value)) %>%
  arrange(n_seqs) %>%
  print(n=15)

bcdf <- bcdf %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames("Sample_Name")

library(vegan)
# generate distance matrix - bray-curtis (default)
#Used sample and iteration numbers from the example on this page..
#https://rdrr.io/cran/vegan/man/avgdist.html#:~:text=Description%20The%20function%20computes%20the%20dissimilarity%20matrix%20of,that%20represents%20the%20average%20of%20multiple%20subsampling%20iterations.
soil_dist <- avgdist(bcdf, sample = 50, iterations = 10)

soil_dist %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(day_a = case_when(str_detect(sample, "non-inoc") ~ 0,
                           str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D7_") ~ 7,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(day_b = case_when(str_detect(name, "non-inoc") ~ 0,
                           str_detect(name, "_D1_") ~ 1,
                           str_detect(name, "_D7_") ~ 7,
                           str_detect(name, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_b = case_when(str_detect(name, "FD") ~ "FD",
                            str_detect(name, "NG") ~ "NG",
                            str_detect(name, "CD") ~ "CD",
                            str_detect(name, "MC") ~ "MC",
                            str_detect(name, "DG") ~ "DG")) %>%
  mutate(pathogen_a = case_when(str_detect(sample, "Bd") ~ "Bd", TRUE ~ "Bsal")) %>%
  mutate(pathogen_b = case_when(str_detect(name, "Bsal") ~ "Bsal", TRUE ~ "Bd")) %>%
  mutate(day = if_else(day_b > day_a, day_b, day_a)) %>%
  filter(site_a == site_b) %>%
  view()
ggplot(aes(x=day, y=value, color=site_a)) +
  geom_line() +
  facet_wrap(~pathogen_a, scales="free_x", nrow=1) +
  theme_bw() +
  theme_classic()

####################### Bd #####################################
#Now need to make Bray-Curtis dissimilarity matrix
#Get df into correct format
bcdf_Bd <- BD_w_d0_df %>%
  replace(is.na(.), "Bd") %>%
  filter(Pathogen_type == "Bd") %>%
  select(Sample_Name, starts_with("ASV"))


bcdf_Bd %>%
  as_tibble() %>%
  pivot_longer(-Sample_Name) %>%
  group_by(Sample_Name) %>%
  summarize(n_seqs = sum(value)) %>%
  arrange(n_seqs) %>%
  print(n=15)

bcdf_Bd <- bcdf_Bd %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames("Sample_Name")

library(vegan)
# generate distance matrix - bray-curtis (default)
#Used sample and iteration numbers from the example on this page..
#https://rdrr.io/cran/vegan/man/avgdist.html#:~:text=Description%20The%20function%20computes%20the%20dissimilarity%20matrix%20of,that%20represents%20the%20average%20of%20multiple%20subsampling%20iterations.
soil_dist_Bd <- avgdist(bcdf_Bd, sample = 50, iterations = 10)

library(RColorBrewer)
library(LaCroixColoR)
#install.packages("nord")
library(nord)
#devtools::install_github("ropenscilabs/ochRe")
#install.packages("ochRe")
library(ochRe)
soil_dist_Bd %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(day_a = case_when(str_detect(sample, "non-inoc") ~ 0,
                           str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(day_b = case_when(str_detect(name, "non-inoc") ~ 0,
                           str_detect(name, "_D1_") ~ 1,
                           str_detect(name, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_b = case_when(str_detect(name, "FD") ~ "FD",
                            str_detect(name, "NG") ~ "NG",
                            str_detect(name, "CD") ~ "CD",
                            str_detect(name, "MC") ~ "MC",
                            str_detect(name, "DG") ~ "DG")) %>%
  mutate(day = if_else(day_b > day_a, day_b, day_a)) %>%
  mutate(diff = abs(day_a-day_b),
         period = if_else(day < 10, "early", "late")) %>%
  filter(site_a == site_b & day_a != day_b) %>%
  filter(!row_number() %in% c(13,14,23,24,33,34,43,44,53,54,63,68,69,74,75)) %>% #Remove rows that were comparing a day 0 sample (non-inoc) to a day 14 sample. Didn't know how else to do this besides viewing the tibble (pipe to View()) and removing those specific rows. Filtering by both conditions removed all non-inoc samples in name column which I didn't want.
  mutate(site_a = factor(site_a, levels = c("FD","DG","MC","NG","CD"))) %>%
  ggplot(aes(x=day, y=value, color=site_a, label = sample)) +
  geom_jitter(width=1.9, height=0, shape = 21, colour = "black", size = 5, stroke = 1, aes(fill = site_a)) +
  scale_fill_manual(values = nord(12)) +
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1))) + #https://stackoverflow.com/questions/40600824/how-to-apply-geom-smooth-for-every-group
  scale_x_continuous(breaks = c(1,14)) +
  scale_y_continuous(limits = c (0.84,0.99)) +
  scale_color_manual(values = nord(12)) +
  #scale_color_manual(values = lacroix_palette("PommeBaya")) +
  labs(x="Days after pathogen introduction",
       y="Bray-Curtis distance to previous day",
       color = "Site") +
  #geom_text() +
  #geom_line() +
  #facet_wrap(~period, scales="free_x") +
  theme_bw() +
  theme_classic() +
  theme(text=element_text(size=20))

#ggsave("Figures/Bd_bray_curtis_x_day.tiff", width=5, height=5)
ggsave("Figures/Bd_bray_curtis_x_day.eps")

##For presentation
soil_dist_Bd %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(day_a = case_when(str_detect(sample, "non-inoc") ~ 0,
                           str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(day_b = case_when(str_detect(name, "non-inoc") ~ 0,
                           str_detect(name, "_D1_") ~ 1,
                           str_detect(name, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_b = case_when(str_detect(name, "FD") ~ "FD",
                            str_detect(name, "NG") ~ "NG",
                            str_detect(name, "CD") ~ "CD",
                            str_detect(name, "MC") ~ "MC",
                            str_detect(name, "DG") ~ "DG")) %>%
  mutate(day = if_else(day_b > day_a, day_b, day_a)) %>%
  mutate(diff = abs(day_a-day_b),
         period = if_else(day < 10, "early", "late")) %>%
  filter(site_a == site_b & day_a != day_b) %>%
  filter(!row_number() %in% c(13,14,23,24,33,34,43,44,53,54,63,68,69,74,75)) %>% #Remove rows that were comparing a day 0 sample (non-inoc) to a day 14 sample. Didn't know how else to do this besides viewing the tibble (pipe to View()) and removing those specific rows. Filtering by both conditions removed all non-inoc samples in name column which I didn't want.
  mutate(site_a = factor(site_a, levels = c("FD","DG","MC","NG","CD"))) %>%
  ggplot(aes(x=day, y=value, color=site_a, label = sample)) +
  geom_point(size = 3,alpha = 0.2) +
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1))) + #https://stackoverflow.com/questions/40600824/how-to-apply-geom-smooth-for-every-group
  scale_x_discrete(limits = c(1,14)) +
  scale_y_continuous(limits = c (0.84,0.99)) +
  scale_color_manual(values = nord(12)) +
  #scale_color_manual(values = lacroix_palette("PommeBaya")) +
  labs(x="Days after pathogen introduction",
       y="Bray-Curtis distance to previous day",
       color = "Site") +
  #geom_text() +
  #geom_line() +
  #facet_wrap(~period, scales="free_x") +
  theme_bw() +
  theme_classic()

ggsave("Figures/Bd_bray_curtis_x_day_presentation.tiff", width=5, height=5)


####################### Bsal #####################################
#Now need to make Bray-Curtis dissimilarity matrix
#Get df into correct format
bcdf_Bsal <- BD_w_d0_df %>%
  replace(is.na(.), "Bsal") %>%
  filter(Pathogen_type == "Bsal") %>%
  select(Sample_Name, starts_with("ASV"))

bcdf_Bsal %>%
  as_tibble() %>%
  pivot_longer(-Sample_Name) %>%
  group_by(Sample_Name) %>%
  summarize(n_seqs = sum(value)) %>%
  arrange(n_seqs) %>%
  print(n=15)

bcdf_Bsal <- bcdf_Bsal %>%
  `row.names<-`(., NULL) %>% 
  column_to_rownames("Sample_Name")

library(vegan)
# generate distance matrix - bray-curtis (default)
#Used sample and iteration numbers from the example on this page..
#https://rdrr.io/cran/vegan/man/avgdist.html#:~:text=Description%20The%20function%20computes%20the%20dissimilarity%20matrix%20of,that%20represents%20the%20average%20of%20multiple%20subsampling%20iterations.
soil_dist_Bsal <- avgdist(bcdf_Bsal, sample = 50, iterations = 10)

soil_dist_Bsal %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(day_a = case_when(str_detect(sample, "non-inoc") ~ 0,
                           str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D7_") ~ 7,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(day_b = case_when(str_detect(name, "non-inoc") ~ 0,
                           str_detect(name, "_D1_") ~ 1,
                           str_detect(name, "_D7_") ~ 7,
                           str_detect(name, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_b = case_when(str_detect(name, "FD") ~ "FD",
                            str_detect(name, "NG") ~ "NG",
                            str_detect(name, "CD") ~ "CD",
                            str_detect(name, "MC") ~ "MC",
                            str_detect(name, "DG") ~ "DG")) %>%
  mutate(day = if_else(day_b > day_a, day_b, day_a)) %>%
  mutate(diff = abs(day_a-day_b),
         period = if_else(day < 10, "early", "late")) %>%
  filter(site_a == site_b & day_a != day_b) %>%
  filter(!row_number() %in% c(2,4,6,10,12,14,17,19,21,22,23,24,28,34,36,38,41,43,45,46,47,48,53,55,57,61,63,65,66,67,71,76,78,80,82,83,84,88,93,95,97,99,100,101,105,110,112,114,117,119,121,123,124,125,129,134,136,138,141,143,145,147,148,149,153,159,161,163,165,167,169,171,172,173,177)) %>% #Remove rows that were comparing a day 0 sample (non-inoc) to a day 14 sample. Didn't know how else to do this besides viewing the tibble (pipe to View()) and removing those specific rows. Filtering by both conditions removed all non-inoc samples in name column which I didn't want.
  mutate(site_a = factor(site_a, levels = c("FD","DG","MC","NG","CD"))) %>%
  ggplot(aes(x=day, y=value, color=site_a, label = sample)) +
  geom_jitter(width=1.7, height=0, shape = 21, colour = "black", size = 5, stroke = 1, aes(fill = site_a)) +
  scale_fill_manual(values = nord(12)) +
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1))) + #https://stackoverflow.com/questions/40600824/how-to-apply-geom-smooth-for-every-group
  scale_x_continuous(breaks = c(1,7,14)) +
  scale_y_continuous(limits = c (0.84,0.99)) +
  scale_color_manual(values = nord(12)) +
  #scale_color_manual(values = lacroix_palette("PommeBaya")) +
  labs(x="Days after pathogen introduction",
       y="Bray-Curtis distance to previous day",
       color="Site") +
  #geom_text() +
  #geom_line() +
  #facet_wrap(~period, scales="free_x") +
  theme_bw() +
  theme_classic() +
  theme(text=element_text(size=20))

#ggsave("Figures/Bsal_bray_curtis_x_day.tiff", width=5, height=5)
ggsave("Figures/Bsal_bray_curtis_x_day.eps")

##For presentation
soil_dist_Bsal %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  mutate(day_a = case_when(str_detect(sample, "non-inoc") ~ 0,
                           str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(day_b = case_when(str_detect(name, "non-inoc") ~ 0,
                           str_detect(name, "_D1_") ~ 1,
                           str_detect(name, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_b = case_when(str_detect(name, "FD") ~ "FD",
                            str_detect(name, "NG") ~ "NG",
                            str_detect(name, "CD") ~ "CD",
                            str_detect(name, "MC") ~ "MC",
                            str_detect(name, "DG") ~ "DG")) %>%
  mutate(day = if_else(day_b > day_a, day_b, day_a)) %>%
  mutate(diff = abs(day_a-day_b),
         period = if_else(day < 10, "early", "late")) %>%
  filter(site_a == site_b & day_a != day_b) %>%
  filter(!row_number() %in% c(2,4,6,10,12,14,17,19,21,22,23,24,28,34,36,38,41,43,45,46,47,48,53,55,57,61,63,65,66,67,71,76,78,80,82,83,84,88,93,95,97,99,100,101,105,110,112,114,117,119,121,123,124,125,129,134,136,138,141,143,145,147,148,149,153,159,161,163,165,167,169,171,172,173,177)) %>% #Remove rows that were comparing a day 0 sample (non-inoc) to a day 14 sample. Didn't know how else to do this besides viewing the tibble (pipe to View()) and removing those specific rows. Filtering by both conditions removed all non-inoc samples in name column which I didn't want.
  mutate(site_a = factor(site_a, levels = c("FD","DG","MC","NG","CD"))) %>%
  ggplot(aes(x=day, y=value, color=site_a, label = sample)) +
  geom_point(size = 3,alpha = 0.2) +
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1))) + #https://stackoverflow.com/questions/40600824/how-to-apply-geom-smooth-for-every-group
  scale_x_discrete(limits = c(1,14)) +
  scale_y_continuous(limits = c (0.84,0.99)) +
  scale_color_manual(values = nord(12)) +
  #scale_color_manual(values = lacroix_palette("PommeBaya")) +
  labs(x="Days after pathogen introduction",
       y="Bray-Curtis distance to previous day",
       color="Site") +
  #geom_text() +
  #geom_line() +
  #facet_wrap(~period, scales="free_x") +
  theme_bw() +
  theme_classic()

ggsave("Figures/Bsal_bray_curtis_x_day_presentation.tiff", width=5, height=5)

#devtools::install_github('thomasp85/gganimate')
#install.packages("remotes")
#remotes::install_github("r-rust/hellorust")
#devtools::install_github("r-rust/gifski")
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "devel")
#BiocManager::install("microbiome")

library(microbiome)
library(ggplot2)
library(gganimate)
install.packages("av")
library(av)
library(gifski)
library(BiocManager)
#install.packages("nord")
library(nord)
#install.packages("gifski")

theme_set(theme_bw())

soil_dist_Bd %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  

BD_Bd <- subset_samples(BD_w_d0, Pathogen_type == "Bd")
sample_data(BD_Bd)

Bdseq <- microbiome::transform(BD_Bd, "compositional")
# Ordinate the data
set.seed(4235421)
# proj <- get_ordination(pseq, "MDS", "bray")
ord <- ordinate(Bdseq, "MDS", "bray")
p <- plot_ordination(Bdseq, ord, color = "Site") +
  geom_point(size = 9, alpha = .7) +
  scale_color_manual(values = nord(12))
p
p + transition_time(Day) +
  labs(title = "Day: {frame_time}") +
  shadow_mark(alpha = 0.3, size = 0.5)
animate(p)

BD_Bsal <- subset_samples(BD_w_d0, Pathogen_type == "Bsal")
sample_data(BD_Bsal)

Bsalseq <- microbiome::transform(BD_Bsal, "compositional")
# Ordinate the data
set.seed(4235421)
# proj <- get_ordination(pseq, "MDS", "bray")
ordBsal <- ordinate(Bsalseq, "MDS", "bray")
pBsal <- plot_ordination(Bsalseq, ordBsal, color = "Site") +
  geom_point(size = 10, alpha = .8) +
  scale_color_manual(values = nord(12))
pBsal
pBsal + transition_time(Day) +
  labs(title = "Day: {frame_time}") +
  shadow_mark(alpha = 0.3, size = 1)
animate(pBsal)
anim_save("Videos/Bsal_BC_x_day.mp4", animation = last_animation(), width = 5, height = 5)

Bdvid <- soil_dist_Bd %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  filter(!sample==str_detect(sample, "non-inoc")) %>%
  mutate(day_a = case_when(str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_a = factor(site_a, levels = c("FD","DG","MC","NG","CD"))) %>%
  ggplot(aes(x=day_a, y=value, color=site_a, label = sample)) +
  geom_point(shape = 21, colour = "black", size = 8, stroke = 0.5, alpha = .7, aes(fill=site_a)) +
  scale_color_manual(values = nord(12)) +
  scale_fill_manual(values = nord(12)) +
  scale_x_discrete(breaks = c(1,14)) +
  scale_y_continuous(limits = c(0.85,1)) +
  labs(fill = "Site",x="Time",y="Bray-Curtis Dissimilarity") +
  transition_time(day_a)

animate(Bdvid,renderer = gifski_renderer(file = "Videos/Bd_BC.gif"),width = 1200, height = 900, res = 300,fps = 9)

Bsalvid <- soil_dist_Bsal %>%
  as.matrix() %>%
  as_tibble(rownames = "sample") %>%
  pivot_longer(-sample) %>%
  filter(sample < name) %>%
  filter(!sample==str_detect(sample, "non-inoc")) %>%
  mutate(day_a = case_when(str_detect(sample, "_D1_") ~ 1,
                           str_detect(sample, "_D14_") ~ 14)) %>%
  mutate(site_a = case_when(str_detect(sample, "FD") ~ "FD",
                            str_detect(sample, "NG") ~ "NG",
                            str_detect(sample, "CD") ~ "CD",
                            str_detect(sample, "MC") ~ "MC",
                            str_detect(sample, "DG") ~ "DG")) %>%
  mutate(site_a = factor(site_a, levels = c("FD","DG","MC","NG","CD"))) %>%
  ggplot(aes(x=day_a, y=value, color=site_a, label = sample)) +
  geom_point(shape = 21, colour = "black", size = 8, stroke = 0.5, alpha = .7, aes(fill=site_a)) +
  scale_color_manual(values = nord(12)) +
  scale_fill_manual(values = nord(12)) +
  scale_x_discrete(breaks = c(1,14)) +
  labs(fill = "Site",x="Time",y="Bray-Curtis Dissimilarity") +
  transition_time(day_a)


animate(Bsalvid,renderer = gifski_renderer(file = "Videos/Bsal_BC.gif"),width = 1200, height = 900, res = 300,fps = 9)

BD_Bsal
distance.matrix <- vegdist(BD_Bsal)

sseq <- microbiome::transform(BD_w_d0, "compositional")
# Ordinate the data
set.seed(4235421)
# proj <- get_ordination(pseq, "MDS", "bray")
ord <- ordinate(sseq, "MDS", "bray")
p <- plot_ordination(sseq, ord, shape = "Pathogen_type", color = "Site") +
  geom_point(size = 5) +
  facet_wrap(~Pathogen_type)
p
p + transition_time(Day) +
  labs(title = "Day: {frame_time}")
animate(p, renderer = av_renderer('animation.mp4'), 
              width = 1280, height = 720, res = 104, fps = 25)

