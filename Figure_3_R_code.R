#Sarah McGrath-Blaser

#Set working directory
#setwd("C:/Users/Sarah/Dropbox (UFL)/Shared_Workspace/Sarah/qPCR/")

#install.packages("janitor")
library(janitor)

################################ Bd ###############################################################
#Read in data from qPCR results
data_D1 <- read.csv("Inputs/012921_Bd_D1_All_Live_Results.csv")
data_D1 <- row_to_names(data_D1, row_number = 5, remove_rows_above = T)
data_D14 <- read.csv("Inputs/012921_Bd_D14_All_Live_Results.csv")
data_D14 <- row_to_names(data_D14, row_number = 5, remove_rows_above = T)
data <- rbind(data_D1, data_D14)
data

#Structure says all columns are characters
str(data)
#Change column Quantity to numeric instead of character
data$Quantity <- as.numeric(data$Quantity)

#Create new column and multiply all values in Quantity by 200 (dilution factor)
data$Quant_mult <- data$Quantity * 200

#Make NAs 0
data[is.na(data)] <- 0
data

#Make new column Day and specify which samples belong to Days 1 and 14.
data$Day <- ifelse(grepl("D1_", data$`Sample Name`),"1", "14")
data
#Specify which are Standards
data$Day[data$Task == "STANDARD"] <- "Standard"

#Make new column and specify what the different treatments are (i.e., Bd pathogen load).
data$Treatment[grepl("1.1", data$`Sample Name`)] <- "10^3"
data$Treatment[grepl("1.2", data$`Sample Name`)] <- "10^3"
data$Treatment[grepl("1.3", data$`Sample Name`)] <- "10^3"
data$Treatment[grepl("2.1", data$`Sample Name`)] <- "10^5"
data$Treatment[grepl("2.2", data$`Sample Name`)] <- "10^5"
data$Treatment[grepl("2.3", data$`Sample Name`)] <- "10^5"

#Make new column and specify soil type.
data$Type[grepl("B", data$`Sample Name`)] <- "Broth"
data$Type[grepl("BdGF", data$`Sample Name`)] <- "Standard"
data$Type[grepl("neg", data$`Sample Name`)] <- "Standard"
data$Type[grepl("FD", data$`Sample Name`)] <- "Fontana Dam"
data$Type[grepl("NG", data$`Sample Name`)] <- "Newfound Gap"
data$Type[grepl("CD", data$`Sample Name`)] <- "Clingman's Dome"
data$Type[grepl("MC", data$`Sample Name`)] <- "Mt Cammerer"
data$Type[grepl("DG", data$`Sample Name`)] <- "Davenport Gap"
data

#Make new column to specify whether soil was autoclaved or non-autoclaved.
data$Soil_state[grepl("_A_", data$`Sample Name`)] <- "Autoclaved"
data$Soil_state[grepl("_NA_", data$`Sample Name`)] <- "Non-autoclaved"

#write to csv file so can more easily look for outliers
write.csv(data, file = "Outputs/qPCR_data_Bd.csv")

#remove rows containing outliers (found 4 that showed no amplification in one of the triplicates)
'%ni%' <- Negate('%in%')
data  <- subset(data, data$`Sample Name`%ni% c(
  "D1_NA_B1.3","D14_A_DG1.2","D14_NA_FD2.2","D14_NA_FD1.3")) 

library(tidyverse)
data
str(data)
data_avg <- data %>%
  as_tibble() %>%
  select(c('Sample Name', Quant_mult, Day, Treatment, Soil_state, Type)) %>%
  na.omit() %>%
  group_by(Type, Soil_state, Day, Treatment) %>%
  summarise(Quantity = mean(Quant_mult)) %>%
  unite(treatment, c("Soil_state","Treatment")) 
data_avg
data_avg$Quantity <- log10(data_avg$Quantity + 1)

data_avg <- as.data.frame(data_avg)

data_avg
colnames(data_avg)[1] <- "Soil_type"

########################### Bd qPCR STATS #####################################
data_Bd_103_NA <- subset(data, Soil_state == "Non-autoclaved" & Treatment == "10^3" & Type != "Broth")
data_Bd_103_NA <- subset(data, Soil_state == "Non-autoclaved" & Treatment == "10^3")
library(car)
# Levene's test with one independent variable
leveneTest(Quant_mult ~ Day, data = data_Bd_103_NA)
#Try to log transform and then analyze
data_Bd_103_NA$logQuant <- log(data_Bd_103_NA$Quant_mult + 1)
leveneTest(logQuant ~ Day, data = data_Bd_103_NA)

hist(data_Bd_103_NA$Quant_mult)
hist(data_Bd_103_NA$logQuant)
qqPlot(data_Bd_103_NA$Quant_mult)
qqPlot(data_Bd_103_NA$logQuant)
shapiro.test(data_Bd_103_NA$Quant_mult)
shapiro.test(data_Bd_103_NA$logQuant)
#Not normal and log transforming does not convince me that it coerces to a normal distribution

#Going with a non-parametric Mann-Whitney U/Wilcoxon-Mann-Whitney test
wilcox.test(Quant_mult ~ Day, data=data_Bd_103_NA) 
#Significantly different
#At a 0.05 signficiance level, we can conclude that the Bd loads between days 1 and 14 are nonidentical (i.e., significantly different)

#Repeat for the other variables
data_Bd_105_NA <- subset(data, Soil_state == "Non-autoclaved" & Treatment == "10^5" & Type != "Broth")
wilcox.test(Quant_mult ~ Day, data=data_Bd_105_NA) 
#Significantly different

data_Bd_103_A <- subset(data, Soil_state == "Autoclaved" & Treatment == "10^3" & Type != "Broth")
hist(data_Bd_103_A$Quant_mult)
shapiro.test(data_Bd_103_A$Quant_mult)
wilcox.test(Quant_mult ~ Day, data=data_Bd_103_A) 
#NOT significantly different

data_Bd_105_A <- subset(data, Soil_state == "Autoclaved" & Treatment == "10^5" & Type != "Broth")
hist(data_Bd_105_A$Quant_mult)
shapiro.test(data_Bd_105_A$Quant_mult)
wilcox.test(Quant_mult ~ Day, data=data_Bd_105_A) 
#NOT significantly different

data_Bd_broth <- subset(data, Type == "Broth")
hist(data_Bd_broth$Quant_mult)
shapiro.test(data_Bd_broth$Quant_mult)
wilcox.test(Quant_mult ~ Day, data=data_Bd_broth) 
#NOT significantly different

##Read in data from live and dead cell counts
library(tidyverse)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(ggplot2)
library(scales)

#Read in files
Bd_D4 <- read.csv("Inputs/Bd_D4_count_data.csv",header = T)
Bd_D7 <- read.csv("Inputs/Bd_D7_count_data.csv",header = T)
Bd_D14 <- read.csv("Inputs/Bd_D14_count_data.csv",header = T)
Bd_D4 <- select(Bd_D4, -starts_with("X") )
Bd_D7 <- select(Bd_D7, -starts_with("X"))
Bd_D14 <- select(Bd_D14, -starts_with("X"))

#Bd_D7 has lowercase green and red under Flour_type so need to change it to match the other days.
Bd_D7$Fluor_type[Bd_D7$Fluor_type == "green"] <- "Green"
Bd_D7$Fluor_type[Bd_D7$Fluor_type == "red"] <- "Red"

all_Bd <- rbind(Bd_D4,Bd_D7,Bd_D14)
all_Bd

#total_all <- merge(data,all_Bd,by=c("Soil_type", "treatment","Day"),all=T)

count_avg <- all_Bd %>%
  as_tibble() %>%
  group_by(Fluor_type,Soil_type, Sterilization_status, Inoc_type, Day) %>%
  summarise(Cell_count = mean(Count)) %>%
  unite(treatment, c("Sterilization_status","Inoc_type"))
count_avg

total <- merge(data_avg,count_avg,by=c("Soil_type", "treatment","Day"),all=T)
total
total$Day <- as.numeric(as.character(total$Day))
str(total)
total <- total %>%
  as_tibble() %>%
  pivot_wider(names_from = Fluor_type, values_from = Cell_count) %>%
  #spread(Fluor_type, Cell_count) %>%
  as.data.frame()
#unite(Cell_count, c("Fluor_type","Cell_count")) %>%
total

NA103 <- subset(total, treatment == "Non-autoclaved_10^3")
NA103

NA103$Day <- factor(NA103$Day, levels = c("1","4","7","14"))
NA103$Soil_type <- factor(NA103$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))
facet.labs <- c("Broth", "FD", "DG", "MC", "NG", "CD")
names(facet.labs) <- c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome")


Bd1 <- ggplot(NA103) +
  geom_col(aes(x=Day, y=Red),
           fill = "darkslateblue", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="cyan4",width = 0.6) +
  geom_point(data = NA103[!is.na(NA103$Quantity), ],aes(x=Day, 
                                                        y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = NA103[!is.na(NA103$Quantity), ],aes(x=Day, 
                                                       y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y.right = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))

Bd1
ggsave("Figures/Bd1_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)


A103 <- subset(total, treatment == "Autoclaved_10^3")
A103
A103$Day <- factor(A103$Day, levels = c("1","4","7","14"))
A103$Soil_type <- factor(A103$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))

Bd2 <- ggplot(A103) +
  geom_col(aes(x=Day, y=Red),
           fill = "darkslateblue", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="cyan4",width = 0.6) +
  geom_point(data = A103[!is.na(A103$Quantity), ],aes(x=Day, 
                                                      y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = A103[!is.na(A103$Quantity), ],aes(x=Day, 
                                                     y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
Bd2
ggsave("Figures/Bd2_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

NA105 <- subset(total, treatment == "Non-autoclaved_10^5")
NA105
NA105$Day <- factor(NA105$Day, levels = c("1","4","7","14"))
NA105$Soil_type <- factor(NA105$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))

Bd3 <- ggplot(NA105) +
  geom_col(aes(x=Day, y=Red),
           fill = "darkslateblue", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="cyan4",width = 0.6) +
  geom_point(data = NA105[!is.na(NA105$Quantity), ],aes(x=Day, 
                                                        y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = NA105[!is.na(NA105$Quantity), ],aes(x=Day, 
                                                       y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.y.right = element_blank(),
        strip.text.x = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
Bd3
ggsave("Figures/Bd3_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

A105 <- subset(total, treatment == "Autoclaved_10^5")
A105
A105$Day <- factor(A105$Day, levels = c("1","4","7","14"))
A105$Soil_type <- factor(A105$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))

Bd4 <- ggplot(A105) +
  geom_col(aes(x=Day, y=Red),
           fill = "darkslateblue", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="cyan4",width = 0.6) +
  geom_point(data = A105[!is.na(A105$Quantity), ],aes(x=Day, 
                                                      y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = A105[!is.na(A105$Quantity), ],aes(x=Day, 
                                                     y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.y.left = element_blank(),
        strip.text.x = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
Bd4
ggsave("Figures/Bd4_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

############################ BSAL #############################################################
#Read in data from qPCR results
data_Bsal_D1 <- read.csv("Inputs/031221_Bsal_D1_All_Live.csv")
data_Bsal_D1 <- row_to_names(data_Bsal_D1, row_number = 5, remove_rows_above = T)
data_Bsal_D7 <- read.csv("Inputs/031221_Bsal_D7_All_Live.csv")
data_Bsal_D7 <- row_to_names(data_Bsal_D7, row_number = 5, remove_rows_above = T)
data_Bsal_D14 <- read.csv("Inputs/040221_Bsal_D14_All_Live.csv")
data_Bsal_D14 <- row_to_names(data_Bsal_D14, row_number = 5, remove_rows_above = T)
data_Bsal <- rbind(data_Bsal_D1, data_Bsal_D7, data_Bsal_D14)
data_Bsal

#Structure says all columns are characters
str(data_Bsal)
#Change column Quantity to numeric instead of character
data_Bsal$Quantity <- as.numeric(data_Bsal$Quantity)

#Create new column and multiply all values in Quantity by 200 (dilution factor)
data_Bsal$Quant_mult <- data_Bsal$Quantity * 200

#Make NAs 0
data_Bsal[is.na(data_Bsal)] <- 0
data_Bsal

#Make new column Day and specify which samples belong to Days 1 and 7.
data_Bsal$Day <- ifelse(grepl("D1_", data_Bsal$`Sample Name`),"1",
                        ifelse(grepl("D7_", data_Bsal$`Sample Name`),"7","14"))
data_Bsal
#Specify which are Standards
data_Bsal$Day[data_Bsal$Task == "STANDARD"] <- "Standard"

#Make new column and specify what the different treatments are (i.e., Bd pathogen load).
data_Bsal$Treatment[grepl("1.1", data_Bsal$`Sample Name`)] <- "10^3"
data_Bsal$Treatment[grepl("1.2", data_Bsal$`Sample Name`)] <- "10^3"
data_Bsal$Treatment[grepl("1.3", data_Bsal$`Sample Name`)] <- "10^3"
data_Bsal$Treatment[grepl("2.1", data_Bsal$`Sample Name`)] <- "10^4"
data_Bsal$Treatment[grepl("2.2", data_Bsal$`Sample Name`)] <- "10^4"
data_Bsal$Treatment[grepl("2.3", data_Bsal$`Sample Name`)] <- "10^4"

#Make new column and specify soil type.
data_Bsal$Type[grepl("B", data_Bsal$`Sample Name`)] <- "Broth"
data_Bsal$Type[grepl("Bsal", data_Bsal$`Sample Name`)] <- "Standard"
data_Bsal$Type[grepl("neg", data_Bsal$`Sample Name`)] <- "Standard"
data_Bsal$Type[grepl("FD", data_Bsal$`Sample Name`)] <- "Fontana Dam"
data_Bsal$Type[grepl("NG", data_Bsal$`Sample Name`)] <- "Newfound Gap"
data_Bsal$Type[grepl("CD", data_Bsal$`Sample Name`)] <- "Clingman's Dome"
data_Bsal$Type[grepl("MC", data_Bsal$`Sample Name`)] <- "Mt Cammerer"
data_Bsal$Type[grepl("DG", data_Bsal$`Sample Name`)] <- "Davenport Gap"
data_Bsal

#Make new column to specify whether soil was autoclaved or non-autoclaved.
data_Bsal$Soil_state[grepl("_A_", data_Bsal$`Sample Name`)] <- "Autoclaved"
data_Bsal$Soil_state[grepl("_NA_", data_Bsal$`Sample Name`)] <- "Non-autoclaved"

#write to csv file so can more easily look for outliers
write.csv(data_Bsal, file = "Outputs/qPCR_data_Bsal.csv")

#Remove outliers
'%ni%' <- Negate('%in%')
data_Bsal <- subset(data_Bsal, data$`Sample Name`%ni% c(
  "D14_NA_B1.1","D14_NA_B2.1","D14_NA_CD1.3","D7_NA_DG1.1","D7_NA_FD1.1",
  "D14_NA_MC1.1","D14_A_MC1.3"
))

data_Bsal
str(data_Bsal)
data_Bsal_avg <- data_Bsal %>%
  as_tibble() %>%
  select(c('Sample Name', Quant_mult, Day, Treatment, Soil_state, Type)) %>%
  na.omit() %>%
  group_by(Type, Soil_state, Day, Treatment) %>%
  summarise(Quantity = mean(Quant_mult)) %>%
  unite(treatment, c("Soil_state","Treatment")) 
data_Bsal_avg

data_Bsal_avg$Quantity <- log10(data_Bsal_avg$Quantity + 1)

data_Bsal_avg <- as.data.frame(data_Bsal_avg)

data_Bsal_avg
colnames(data_Bsal_avg)[1] <- "Soil_type"

########################### Bsal qPCR STATS #####################################
data_Bsal_103_NA <- subset(data_Bsal, Soil_state == "Non-autoclaved" & Treatment == "10^3" & Type != "Broth")
data_Bsal_103_NA

# Levene's test with one independent variable
leveneTest(Quant_mult ~ Day, data = data_Bsal_103_NA)
#Try to log transform and then analyze
data_Bsal_103_NA$logQuant <- log(data_Bsal_103_NA$Quant_mult + 1)
leveneTest(logQuant ~ Day, data = data_Bsal_103_NA)

hist(data_Bsal_103_NA$Quant_mult)
hist(data_Bsal_103_NA$logQuant)
qqPlot(data_Bsal_103_NA$Quant_mult)
qqPlot(data_Bsal_103_NA$logQuant)
shapiro.test(data_Bsal_103_NA$Quant_mult)
shapiro.test(data_Bsal_103_NA$logQuant)
#Not normal and log transforming does not convince me that it coerces to a normal distribution

#Going with a non-parametric Kruskal-Wallis test because for Bsal we have 3 levels/days instead of two
kruskal.test(Quant_mult ~ Day, data=data_Bsal_103_NA) 
#Significantly different
library(FSA)
dunnTest(Quant_mult ~ Day, data=data_Bsal_103_NA, method="bh")
#Day 1 is sig diff from days 7 and 14 but 7 and 14 aren't sig diff from one another

#Repeat for the other variables
data_Bsal_104_NA <- subset(data_Bsal, Soil_state == "Non-autoclaved" & Treatment == "10^4" & Type != "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_104_NA) 
#Significantly different
dunnTest(Quant_mult ~ Day, data=data_Bsal_104_NA, method="bh")
#Day 1 is sig diff from days 7 and 14 but 7 and 14 aren't sig diff from one another


data_Bsal_103_A <- subset(data_Bsal, Soil_state == "Autoclaved" & Treatment == "10^3" & Type != "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_103_A) 
#NOT significantly different

data_Bsal_104_A <- subset(data_Bsal, Soil_state == "Autoclaved" & Treatment == "10^4" & Type != "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_104_A) 
#Significantly different
dunnTest(Quant_mult ~ Day, data=data_Bsal_104_A, method="bh")
#Day 1 is sig diff from days 7 and 14 but 7 and 14 aren't sig diff from one another

data_Bsal_broth <- subset(data_Bsal, Type == "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_broth) 
#Significantly different
dunnTest(Quant_mult ~ Day, data=data_Bsal_broth, method="bh")
######################################################################

#Read in files
Bsal_D1 <- read.csv("Inputs/Bsal_D1_count_data.csv",header = T)
Bsal_D4 <- read.csv("Inputs/Bsal_D4_count_data.csv",header = T)
Bsal_D14 <- read.csv("Inputs/Bsal_D14_count_data.csv",header = T)

all_Bsal <- rbind(Bsal_D1,Bsal_D4,Bsal_D14)
all_Bsal

count_Bsal_avg <- all_Bsal %>%
  as_tibble() %>%
  group_by(Fluor_type,Soil_type, Sterilization_status, Inoc_type, Day) %>%
  summarise(Cell_count = mean(Count)) %>%
  unite(treatment, c("Sterilization_status","Inoc_type"))
count_Bsal_avg

total_Bsal <- merge(data_Bsal_avg,count_Bsal_avg,by=c("Soil_type", "treatment","Day"),all=T)
total_Bsal
total_Bsal$Day <- as.numeric(as.character(total_Bsal$Day))
str(total_Bsal)
total_Bsal <- total_Bsal %>%
  as_tibble() %>%
  spread(Fluor_type, Cell_count) %>%
  as.data.frame()
#unite(Cell_count, c("Fluor_type","Cell_count")) %>%

total_Bsal

Bsal_NA103 <- subset(total_Bsal, treatment == "Non-autoclaved_10^3")
Bsal_NA103
Bsal_NA103$Day <- factor(Bsal_NA103$Day, levels = c("1","4","7","14"))
Bsal_NA103$Soil_type <- factor(Bsal_NA103$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))
facet.labs <- c("Broth", "FD", "DG", "MC", "NG", "CD")
names(facet.labs) <- c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome")

Bsal1 <- ggplot(Bsal_NA103) +
  geom_col(aes(x=Day, y=Red),
           fill = "seagreen4", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="seagreen3",width = 0.6) +
  geom_point(data = Bsal_NA103[!is.na(Bsal_NA103$Quantity), ],aes(x=Day, 
                                                                  y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = Bsal_NA103[!is.na(Bsal_NA103$Quantity), ],aes(x=Day, 
                                                                 y=Quantity*150), group=1, 
            size=1.5,linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "",expand = c(0,0)) +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y.right = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
Bsal1
ggsave("Figures/Bsal1_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

Bsal_A103 <- subset(total_Bsal, treatment == "Autoclaved_10^3")
Bsal_A103
Bsal_A103$Day <- factor(Bsal_A103$Day, levels = c("1","4","7","14"))
Bsal_A103$Soil_type <- factor(Bsal_A103$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))

Bsal2 <- ggplot(Bsal_A103) +
  geom_col(aes(x=Day, y=Red),
           fill = "seagreen4", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="seagreen3",width = 0.6) +
  geom_point(data = Bsal_A103[!is.na(Bsal_A103$Quantity), ],aes(x=Day, 
                                                                y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = Bsal_A103[!is.na(Bsal_A103$Quantity), ],aes(x=Day, 
                                                               y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y.left = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
Bsal2
ggsave("Figures/Bsal2_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

Bsal_NA104 <- subset(total_Bsal, treatment == "Non-autoclaved_10^4")
Bsal_NA104
Bsal_NA104$Day <- factor(Bsal_NA104$Day, levels = c("1","4","7","14"))
Bsal_NA104$Soil_type <- factor(Bsal_NA104$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))

Bsal3 <- ggplot(Bsal_NA104) +
  geom_col(aes(x=Day, y=Red),
           fill = "seagreen4", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="seagreen3",width = 0.6) +
  geom_point(data = Bsal_NA104[!is.na(Bsal_NA104$Quantity), ],aes(x=Day, 
                                                                  y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = Bsal_NA104[!is.na(Bsal_NA104$Quantity), ],aes(x=Day, 
                                                                 y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.y.right = element_blank(),
        strip.text.x = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
Bsal3
ggsave("Figures/Bsal3_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

Bsal_A104 <- subset(total_Bsal, treatment == "Autoclaved_10^4")
Bsal_A104
Bsal_A104$Day <- factor(Bsal_A104$Day, levels = c("1","4","7","14"))
Bsal_A104$Soil_type <- factor(Bsal_A104$Soil_type, levels = c("Broth", "Fontana Dam","Davenport Gap","Mt Cammerer", "Newfound Gap", "Clingman's Dome"))

Bsal4 <- ggplot(Bsal_A104) +
  geom_col(aes(x=Day, y=Red),
           fill = "seagreen4", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="seagreen3",width = 0.6) +
  geom_point(data = Bsal_A104[!is.na(Bsal_A104$Quantity), ],aes(x=Day, 
                                                                y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = Bsal_A104[!is.na(Bsal_A104$Quantity), ],aes(x=Day, 
                                                               y=Quantity*150), group=1, 
            size = 1.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = ""), limits = c(0,2200),
                     name = "") +
  facet_grid(~ Soil_type,labeller = labeller(Soil_type=facet.labs)) +
  xlab("") +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=20),
        axis.text.y.left = element_blank(),
        strip.text.x = element_blank()) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))

Bsal4
ggsave("Figures/Bsal4_fig3.tiff",units = "in",width = 9.5,height = 5,dpi=300)

################ Plots combined with elevation #################################

Elevation <- c(0,656,744,1278,1577,2027)
Site <- c("Control","Fontana Dam","Davenport Gap","Mt Cammerer","Newfound Gap","Clingman's Dome")
df<- data.frame(Elevation, Site)
df
df$Site <- factor(df$Site, levels = c("Control","Fontana Dam","Davenport Gap","Mt Cammerer","Newfound Gap","Clingman's Dome"))

a <- ggplot(data=df, aes(x=Site, y=Elevation)) +
  geom_bar(stat = "identity",width= 0.8) +
  theme_bw() +
  theme_classic()
a

b <- ggplot(NA103) +
  geom_col(aes(x=Day, y=Red),
           fill = "darkslateblue", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="cyan4",width = 0.6) +
  geom_point(data = NA103[!is.na(NA103$Quantity), ],aes(x=Day, 
                                                        y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = NA103[!is.na(NA103$Quantity), ],aes(x=Day, 
                                                       y=Quantity*150), group=1, 
            size = 0.5, linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = "Zoospore quantity (log)"), limits = c(0,2200),
                     name = "Cell count") +
  facet_grid(~ Soil_type) +
  theme_bw() +
  theme_classic() +
  #theme(text = element_text(size=25)) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))
b

#Arrenge figures
library("cowplot")

ggdraw()+
  draw_plot(a,x=0,y=0.61,width=0.865,height = 0.3)+#less width than the other, to match columns
  draw_plot(b,x=0,y=0,width = 0.9,height = 0.7)+
  draw_plot_label(label=c("a","b"),
                  x=c(0.02,0.02),
                  y=c(0.98,0.72))
ggsave("Fig-try1.tiff",units = "in",width = 9,height = 5,dpi=300)

c <- ggplot(Bsal_NA103) +
  geom_col(aes(x=Day, y=Red),
           fill = "seagreen4", width = 0.2) +
  geom_col(aes(x=Day, y=Green),
           alpha=0.3, fill="seagreen3",width = 0.6) +
  geom_point(data = Bsal_NA103[!is.na(Bsal_NA103$Quantity), ],aes(x=Day, 
                                                                  y=Quantity*150), group=1, 
             size = 4) +
  geom_line(data = Bsal_NA103[!is.na(Bsal_NA103$Quantity), ],aes(x=Day, 
                                                                 y=Quantity*150), group=1, 
            size=0.5,linetype = "dashed", color="black") + 
  scale_y_continuous(sec.axis = sec_axis(~./150, name = "Zoospore quantity (log)"), limits = c(0,2200),
                     name = "Cell count",expand = c(0,0)) +
  facet_grid(~ Soil_type) +
  theme_bw() +
  theme_classic() +
  #theme(text = element_text(size=25)) +
  scale_color_manual(name="Type",values = c("darkslateblue","cyan4"))

ggdraw()+
  draw_plot(a,x=0,y=0.61,width=0.865,height = 0.3)+#less width than the other, to match columns
  draw_plot(c,x=0,y=0,width = 0.9,height = 0.7)+
  draw_plot_label(label=c("a","b"),
                  x=c(0.02,0.02),
                  y=c(0.98,0.72))
ggsave("Fig-try2.tiff",units = "in",width = 9,height = 5,dpi=300)


###################### STATS ##########################
#Is there a significant difference between days for qPCR data?
#Is there a significant difference between days for live/dead cell count data?

########################### Live/Dead STATS #####################################
data_Bsal_103_NA <- subset(data_Bsal, Soil_state == "Non-autoclaved" & Treatment == "10^3" & Type != "Broth")
data_Bsal_103_NA

# Levene's test with one independent variable
leveneTest(Quant_mult ~ Day, data = data_Bsal_103_NA)
#Try to log transform and then analyze
data_Bsal_103_NA$logQuant <- log(data_Bsal_103_NA$Quant_mult + 1)
leveneTest(logQuant ~ Day, data = data_Bsal_103_NA)

hist(data_Bsal_103_NA$Quant_mult)
hist(data_Bsal_103_NA$logQuant)
qqPlot(data_Bsal_103_NA$Quant_mult)
qqPlot(data_Bsal_103_NA$logQuant)
shapiro.test(data_Bsal_103_NA$Quant_mult)
shapiro.test(data_Bsal_103_NA$logQuant)
#Not normal and log transforming does not convince me that it coerces to a normal distribution

#Going with a non-parametric Kruskal-Wallis test because for Bsal we have 3 levels/days instead of two
kruskal.test(Quant_mult ~ Day, data=data_Bsal_103_NA) 
#Significantly different
library(FSA)
dunnTest(Quant_mult ~ Day, data=data_Bsal_103_NA, method="bh")
#Day 1 is sig diff from days 7 and 14 but 7 and 14 aren't sig diff from one another

#Repeat for the other variables
data_Bsal_104_NA <- subset(data_Bsal, Soil_state == "Non-autoclaved" & Treatment == "10^4" & Type != "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_104_NA) 
#Significantly different
dunnTest(Quant_mult ~ Day, data=data_Bsal_104_NA, method="bh")
#Day 1 is sig diff from days 7 and 14 but 7 and 14 aren't sig diff from one another


data_Bsal_103_A <- subset(data_Bsal, Soil_state == "Autoclaved" & Treatment == "10^3" & Type != "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_103_A) 
#NOT significantly different

data_Bsal_104_A <- subset(data_Bsal, Soil_state == "Autoclaved" & Treatment == "10^4" & Type != "Broth")
kruskal.test(Quant_mult ~ Day, data=data_Bsal_104_A) 
#Significantly different
dunnTest(Quant_mult ~ Day, data=data_Bsal_104_A, method="bh")
#Day 1 is sig diff from days 7 and 14 but 7 and 14 aren't sig diff from one another
####################### Bd LIVE/DEAD cell STATS ################################
# Levene's test with one independent variable
leveneTest(Green ~ Day, data = NA103)
#p = 0.9703, > 0.05 so groups are equally variable

hist(NA103$Green)
qqPlot(NA103$Green)
shapiro.test(NA103$Green)
#Normal distribution

Live_anova <- aov(Green ~ Day, data = NA103)
summary(Live_anova)
#p = 0.645, the number of live cells (as a proxy of microbial activity) not significantly different between days

#Just checking to see what a kruskal-wallis test shows
kruskal.test(Green ~ Day, data = NA103)
#Also not significant

# Levene's test with one independent variable
leveneTest(Red ~ Day, data = NA103)
#p = 0.8496, > 0.05 so groups are equally variable

hist(NA103$Red)
qqPlot(NA103$Red)
shapiro.test(NA103$Red)
#Not normally distributed (left skewed)

kruskal.test(Red ~ Day, data = NA103) 
#Not significantly different

# Levene's test with one independent variable
leveneTest(Green ~ Day, data = NA105)
#p = 0.9187, > 0.05 so groups are equally variable

hist(NA105$Green)
qqPlot(NA105$Green)
shapiro.test(NA105$Green)
#Not normally distributed

kruskal.test(Green ~ Day, data = NA105)
#Not significant

# Levene's test with one independent variable
leveneTest(Red ~ Day, data = NA105)
#p = 0.2426, > 0.05 so groups are equally variable

hist(NA105$Red)
qqPlot(NA105$Red)
shapiro.test(NA105$Red)
#Not normally distributed (left skewed)

kruskal.test(Red ~ Day, data = NA105) 
#Not significantly different

#Should see a difference between microbially active and inactive soils
kruskal.test(Red ~ treatment, data = total)
#Not different
kruskal.test(Green ~ treatment, data = total)
dunnTest(Green ~ treatment, data = total, method="bh")
#Yes, sig diff between all autoclaved-non-autoclaved comparisons!


#### Bsal LIVE/DEAD cell STATS ####
# Levene's test with one independent variable
leveneTest(Green ~ Day, data = Bsal_NA103)
#p = 0.8335, > 0.05 so groups are equally variable

hist(Bsal_NA103$Green)
qqPlot(Bsal_NA103$Green)
shapiro.test(Bsal_NA103$Green)
#Not normally distributed

kruskal.test(Green ~ Day, data = Bsal_NA103)
#Not significant

# Levene's test with one independent variable
leveneTest(Red ~ Day, data = Bsal_NA103)
#p = 0.3883, > 0.05 so groups are equally variable

hist(Bsal_NA103$Red)
qqPlot(Bsal_NA103$Red)
shapiro.test(Bsal_NA103$Red)
#Not normally distributed (left skewed)

kruskal.test(Red ~ Day, data = Bsal_NA103) 
#Not significantly different

# Levene's test with one independent variable
leveneTest(Green ~ Day, data = Bsal_NA104)
#p = 0.4479, > 0.05 so groups are equally variable

hist(Bsal_NA104$Green)
qqPlot(Bsal_NA104$Green)
shapiro.test(Bsal_NA104$Green)
#Not normally distributed

kruskal.test(Green ~ Day, data = Bsal_NA104)
#Not significant

# Levene's test with one independent variable
leveneTest(Red ~ Day, data = Bsal_NA104)
#p = 0.4137, > 0.05 so groups are equally variable

hist(Bsal_NA104$Red)
qqPlot(Bsal_NA104$Red)
shapiro.test(Bsal_NA104$Red)
#Not normally distributed (left skewed)

kruskal.test(Red ~ Day, data = Bsal_NA104) 
#Not significantly different

#Should see a difference between microbially active and inactive soils
kruskal.test(Red ~ treatment, data = total_Bsal)
#Not different
kruskal.test(Green ~ treatment, data = total_Bsal)
dunnTest(Green ~ treatment, data = total_Bsal, method="bh")
#Yes, sig diff between all autoclaved-non-autoclaved comparisons!




