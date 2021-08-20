## FRASER ET AL, FIGURE 1

## LOAD LIBRARIES

library(tidyverse)

## FIGURE 1A ##

#### IMPORT DATA TABLE ####
Figure_1A <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Fraser_et_al_Figure1A.rds")
View(Figure_1A)

#### BAR PLOT WITH ERROR BARS ####
Figure_1A_Plot <- ggplot(Figure_1A, aes(x = reorder(gene.id, -gene.proportion), y = gene.proportion)) +
  geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper), width = 0.2) +
  geom_bar(stat = 'identity', fill = "black", width = 0.8) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.65), breaks = c(0, 0.2, 0.4, 0.6)) +
  ylab("Proportion of Patients") +
  theme(
    axis.text.x = element_text(face = "bold", 
                               size = 12, 
                               color = "black", 
                               angle = 90, 
                               hjust = 1, 
                               vjust = 0.5, 
                               margin=margin(10,0,0,0)),
    axis.text.y = element_text(face = "bold",
                               size = 30,
                               color = "black"),
    axis.title = element_text(face = 'bold', 
                              size = 36,
                              margin = margin(t=0,r=10,b=0,l=0)),
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black",
                                fill=NA,
                                size=2)
  )
Figure_1A_Plot

## FIGURE 1B

Figure_1B <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Fraser_et_al_Figure1B.rds")
View(Figure_1B)

#### DEFINE COLOURS BASED ON PCAWG TEMPLATE ####

# Define all colours 
# Coding SNV mutation subtypes & consequences 
nonsynonymous <- 'pink'
synonymous <- '#FFD700'
stop.gain <- '#8B4789'
stop.loss <- '#DA70D6'
indel.frameshift <- '#FF8C00'
indel.nonframeshift <- '#003366'
splicing <- '#00CED1'
# Non-coding SNV mutation subtypes, consequences & gene types
non.coding <- 'dodgerblue4'
promoter <- '#4C191E'
enhancer <- '#7F000F'
operator <- '#A84955'
silencer <- '#E78A96'
insulator <- '#FFC1C9'
lncRNA <- '#331900'
sncRNA <- '#594027'
tRNA <- '#A87849'
rRNA <- '#E7B98A'
miRNA <- '#FFE0C1'
utr5.utr3 <- '#1A1A1A'
intronic <- '#4D4D4D'
intergenic <- '#7F7F7F'
telomeres <- '#B3B3B3'
# Structral variant mutation subtypes
cna.gain <- '#FF0000'
cna.loss <- '#0000FF'
inversion <- '#FFA500'
transposition <- '#7300E7'
translocation <- '#458B00'

#### BAR PLOT WITH ERROR BARS ####

Figure_1B_Plot <- ggplot(Figure_1B, 
                                     aes(x = reorder(Abb.Code, -proportion.overall), y = proportion.overall)) +
  geom_errorbar(aes(ymin = proportion.overall.ci.lower, ymax = proportion.overall.ci.upper), width = 0.2) +
  geom_bar(aes(fill = Abberation.Type), position = "dodge", stat = 'identity', width = 0.8) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.65), breaks = c(0, 0.2, 0.4, 0.6)) +
  ylab("Proportion of Patients") +
  scale_x_discrete(labels = Figure_1B$gene.id) +
  scale_fill_manual(values = c(transposition,  nonsynonymous, cna.gain, cna.loss,  translocation, inversion, non.coding)) + # Added PCAWG colours
  guides(fill = guide_legend("MUTATION CLASS")) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10, 
                               color = "black", angle = 90, hjust = 1, vjust = 0.5, margin=margin(10,0,0,0)),
    axis.text.y = element_text(face = "bold", size = 30, color = "black"),
    axis.title = element_text(face = 'bold', size = 36, margin = margin(t=0,r=10,b=0,l=0)),
    axis.title.x = element_blank(),
    legend.title = element_text(face = "bold", size = 24, hjust = 0.5),
    legend.text = element_text(size = 18),
    legend.position = c(0.8, 0.4),
    legend.spacing.y = unit(1, units = "cm"), 
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
Figure_1B_Plot
