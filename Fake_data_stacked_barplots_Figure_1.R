library(reshape2)
library(ggplot2)
library(RColorBrewer)

#Note: plots saved as .eps files and colors changed to reflect native microbial community picture in inkscape

data <-  data.frame(A = c(5,3,2,2,2),        # Create example data
                    B = c(8,2,1,2,1),
                    C = c(6, 4, 3, 1, 1))
data                               # Print example data

data_long <- data                  # Converting data from wide to long format
data_long$subgroup <- as.factor(1:nrow(data_long))
data_long <- melt(data_long, id.vars = "subgroup")
data_long    

Inhibition <- ggplot(data_long,                  # Stacked barplot using ggplot2
       aes(x = variable,
           y = value,
           fill = subgroup)) +
  geom_bar(stat = "identity",position = "fill",color = "black") +
  scale_fill_brewer(palette = "Greens") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position = "none")

Inhibition
ggsave(file = "Figures/MicrobialInvasion_Inhibition.eps") 


data2 <- data.frame(A = c(5,3,2,2,2),        # Create example data
                   B = 6:2,
                   C = c(5, 1, 1, 1, 3))
data2                               # Print example data

data_long2 <- data2                  # Converting data from wide to long format
data_long2$subgroup <- as.factor(1:nrow(data_long2))
data_long2 <- melt(data_long2, id.vars = "subgroup")
data_long2 

Persistence <- ggplot(data_long2,                  # Stacked barplot using ggplot2
                      aes(x = variable,
                          y = value,
                          fill = subgroup)) +
  geom_bar(stat = "identity",position = "fill",color = "black") +
  scale_fill_brewer(palette = "BuPu") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position = "none")
Persistence
ggsave(file = "Figures/MicrobialInvasion_Persistence.eps") 

data3 <- data.frame(A = c(5,3,2,2,2),        # Create example data
                    B = c(2,1,2,3,2),
                    C = c(2, 2, 1, 1, 2))
data3                               # Print example data

data_long3 <- data3                  # Converting data from wide to long format
data_long3$subgroup <- as.factor(1:nrow(data_long3))
data_long3 <- melt(data_long3, id.vars = "subgroup")
data_long3 

Facilitation <- ggplot(data_long3,                  # Stacked barplot using ggplot2
                      aes(x = variable,
                          y = value,
                          fill = subgroup)) +
  geom_bar(stat = "identity",position = "fill",color = "black") +
  scale_fill_brewer(palette = "YlOrRd") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position = "none")
Facilitation
ggsave(file = "Figures/MicrobialInvasion_Facilitation.eps") 
