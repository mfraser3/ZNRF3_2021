#### FRASER ET AL - FIGURE S03 ####

#### LOAD LIBRARIES ####
library(tidyverse)
library(cowplot)


#### DEFINE DIRECTORIES ####
input_dir <- file.path("/Users/michaelfraser/OneDrive/Work/Projects/ZNRF3/FINAL DATA TABLES/Figures/Figure 1/")
output_dir <- file.path("/Users/michaelfraser/OneDrive/Work/Projects/ZNRF3/FINAL FIGURES/Draft_11/")

Figure_1B_Local <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Figure_1B_Local.rds")
Figure_1B_Mets <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Figure_1B_Mets.rds")

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

#### BAR PLOT WITH ERROR BARS  - LOCALIZED ####

DRIVERS.BY.ABERRATION.LOCAL.PLOT <- ggplot(Figure_1B_Local, 
                                           aes(x = reorder(Abb.Code, -proportion.localized), y = proportion.localized)) +
  geom_errorbar(aes(ymin = proportion.localized.ci.lower, ymax = proportion.localized.ci.upper), width = 0.2) +
  geom_bar(aes(fill = Abberation.Type), position = "dodge", stat = 'identity', width = 0.8) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.65), breaks = c(0, 0.2, 0.4, 0.6)) +
  ylab("Proportion of Patients") +
  ggtitle("Localized Prostate Cancer") +
  scale_x_discrete(labels = Figure_1B_Local$gene.id) +
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
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    plot.title = element_text(size = 30, color = "black", face = "bold")
  )
DRIVERS.BY.ABERRATION.LOCAL.PLOT

#### BAR PLOT WITH ERROR BARS - METASTATIC ####

DRIVERS.BY.ABERRATION.METS.PLOT <- ggplot(Figure_1B_Mets, 
                                          aes(x = reorder(Abb.Code, -proportion.mets), y = proportion.mets)) +
  geom_errorbar(aes(ymin = proportion.mets.ci.lower, ymax = proportion.mets.ci.upper), width = 0.2) +
  geom_bar(aes(fill = Abberation.Type), position = "dodge", stat = 'identity', width = 0.8) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.83), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  ylab("Proportion of Patients") +
  ggtitle("Metastatic Prostate Cancer") +
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
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    plot.title = element_text(size = 30, color = "black", face = "bold")
  )
DRIVERS.BY.ABERRATION.METS.PLOT

## MERGE PLOTS ####
ZNRF3_FigureS03_by_Cohort_Merged <- plot_grid(DRIVERS.BY.ABERRATION.LOCAL.PLOT, DRIVERS.BY.ABERRATION.METS.PLOT, axis = c("l"), ncol = 1)
ZNRF3_FigureS03_by_Cohort_Merged
