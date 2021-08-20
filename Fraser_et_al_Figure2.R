## FRASER ET AL, FIGURE 2

## LOAD LIBRARIES ####

library(easyGgplot2)
library(tidyverse)
library(cowplot)
library(ggrepel)

#### DEFINE DIRECTORIES ####
input_dir <- file.path("/Users/michaelfraser/OneDrive/Work/Projects/ZNRF3/FINAL DATA TABLES/")
output_dir <- file.path("/Users/michaelfraser/OneDrive/Work/Projects/ZNRF3/FINAL FIGURES/Draft_11/")

# IMPORT DATA
Figure2 <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Fraser_et_al_Figure2.rds")
Figure2 <- Figure2[1:113,]

#### DEFINE COLOURS BASED ON PCAWG ####

# Define all colours 
# Coding SNV mutation subtypes & consequences 
nonsynonymous <- '#698B69'
synonymous <- '#FFD700'
stop.gain <- '#8B4789'
stop.loss <- '#DA70D6'
indel.frameshift <- '#FF8C00'
indel.nonframeshift <- '#003366'
splicing <- '#00CED1'
# Non-coding SNV mutation subtypes, consequences & gene types
non.coding <- '#A80015'
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

Mutation_Colours <- c("purple", "chartreuse3", "orange", "dodgerblue4")

#### BARPLOT OF DELTA PROPORTIONS ####

ADJUSTED.DELTA.PROPORTIONS <- ggplot(data = Figure2, aes(x=reorder(Abb.Code, adjusted.delta.proportion), y = adjusted.delta.proportion)) +
  ylab("Adjusted \U0394 Proportion") +
  labs(fill = "MUTATION TYPE") +
  scale_x_discrete(labels = c('CHD1',
                              'RB1',
                              'NKX3-1',
                              'TTN',
                              'CDK12',
                              'KDM6A',
                              'chr3:125 Mbp',
                              'MSH2',
                              'MED12',
                              'PCDH15',
                              'chr3:129 Mbp',
                              'chr6:97 Mbp',
                              'SPOP',
                              'CTNNB1',
                              'KDM6A',
                              'chr7:145019604',
                              'BRCA1',
                              'PRKDC',
                              'chr22:31826396',
                              'PIK3CA',
                              'chr16:64910033',
                              'chr2:131394079',
                              'chr13:58055742',
                              'MLH1',
                              'SYNE1',
                              'BRAF',
                              'chr4:148 Mbp',
                              'BRCA1',
                              'GNAS',
                              'ASXL2',
                              'chr4:39684557',
                              'MYC',
                              'NCOR2',
                              'BRAF',
                              'ERCC2',
                              'MSH2',
                              'NCOR1',
                              'AKT1',
                              'MSH6',
                              'ATM',
                              'IDH1',
                              'AXL',
                              'PTEN',
                              'ETV5',
                              'MSH2',
                              'MED12',
                              'CDKN2A',
                              'ERCC2',
                              'CDKN2A',
                              'PRB4',
                              'HRAS',
                              'ASXL2',
                              'TPTE',
                              'RB1',
                              'ADAMTS20',
                              'PRKDC',
                              'ATM',
                              'ZFHX3',
                              'RYR2',
                              'FAT3',
                              'MSH6',
                              'HDAC4',
                              'QRICH2',
                              'NCOR1',
                              'SPTA1',
                              'CHD1',
                              'ASXL2',
                              'ETV4',
                              'CSMD3',
                              'KDM6A',
                              'ERG',
                              'chr17:25 Mbp',
                              'FOXA1',
                              'KMT2C',
                              'APC',
                              'PTEN',
                              'ZBTB16',
                              'BRCA2',
                              'AHNAK2',
                              'ETV1',
                              'KMT2D',
                              'ATM',
                              'ETV1',
                              'FOXA1',
                              'HDAC4',
                              'TP53',
                              'ZNRF3',
                              'chr10:89 Mbp',
                              'CDKN2A',
                              'APC',
                              'MUC16',
                              'ZNRF3',
                              'AKT1',
                              'CTNNB1',
                              'chr21:42 Mbp',
                              'CDH1',
                              'AR',
                              'chr7:61 Mbp',
                              'CCND1',
                              'ETV5',
                              'NKX3-1',
                              'FOXA1',
                              'BRAF',
                              'NCOR1',
                              'BRCA2',
                              'PRKDC',
                              'PTEN',
                              'TP53',
                              'ZFHX3',
                              'RB1',
                              'MYC',
                              'TP53',
                              'AR')) +
  ylim(-0.25,0.75) +
  geom_errorbar(inherit.aes = TRUE, stat = "identity",
                aes(ymin = adjusted.delta.proportion.lower, ymax = adjusted.delta.proportion.ci.upper)) +
  geom_bar(fill = "black", stat = "identity", color = "white", width = 1) +
  coord_flip() +
  geom_bar(fill = "black", stat = "identity", color = "white", width = 1) +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold",
                                size = 28),
    axis.line = element_line(),
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    axis.text.y = element_text(face = "bold", size= 16, color = "black"),
    legend.position = c(0,0),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(0.1,0.5,0,0.1), "cm")
  )
ADJUSTED.DELTA.PROPORTIONS

## DELTA DELTA PROPORTIONS ####
DELTA.DELTA.PROPORTIONS <- ggplot(data = Figure2, aes(x=reorder(Abb.Code, proportion.delta), y = observed.delta.minus.expected.delta)) +
  ylab("Observed - Expected\n\U0394 Proportion") +
  labs(fill = "MUTATION TYPE") +
  scale_x_discrete(labels = c('CHD1',
                              'RB1',
                              'NKX3-1',
                              'chr3:125 Mbp',
                              'chr3:129 Mbp',
                              'chr6:97 Mbp',
                              'chr7:145019604',
                              'PRKDC',
                              'chr22:31826396',
                              'PIK3CA',
                              'MSH2',
                              'chr13:58055742',
                              'chr16:64910033',
                              'chr2:131394079',
                              'chr4:148 Mbp',
                              'SPOP',
                              'chr4:39684557',
                              'MYC',
                              'BRAF',
                              'ERCC2',
                              'MSH2',
                              'AXL',
                              'PTEN',
                              'ETV5',
                              'MED12',
                              'IDH1',
                              'CDKN2A',
                              'PRB4',
                              'CDKN2A',
                              'ERCC2',
                              'MLH1',
                              'BRAF',
                              'HRAS',
                              'AKT1',
                              'ATM',
                              'TPTE',
                              'MSH6',
                              'MED12',
                              'NCOR1',
                              'NCOR2',
                              'KDM6A',
                              'RB1',
                              'BRCA1',
                              'ADAMTS20',
                              'NCOR1',
                              'ASXL2',
                              'ASXL2',
                              'ETV4',
                              'GNAS',
                              'PCDH15',
                              'QRICH2',
                              'CDK12',
                              'SPTA1',
                              'PRKDC',
                              'KDM6A',
                              'ERG',
                              'KDM6A',
                              'CTNNB1',
                              'FAT3',
                              'ZFHX3',
                              'ATM',
                              'chr17:25 Mbp',
                              'RYR2',
                              'FOXA1',
                              'PTEN',
                              'ETV1',
                              'BRCA1',
                              'HDAC4',
                              'HDAC4',
                              'KMT2D',
                              'FOXA1',
                              'APC',
                              'CSMD3',
                              'BRCA2',
                              'TP53',
                              'ZNRF3',
                              'AHNAK2',
                              'KMT2C',
                              'ASXL2',
                              'chr10:89 Mbp',
                              'SYNE1',
                              'MSH2',
                              'chr21:42 Mbp',
                              'CHD1',
                              'MSH6',
                              'chr7:61 Mbp',
                              'MUC16',
                              'AR',
                              'TTN',
                              'ZBTB16',
                              'CDKN2A',
                              'ETV1',
                              'CCND1',
                              'ATM',
                              'AKT1',
                              'NKX3-1',
                              'ZNRF3',
                              'APC',
                              'CTNNB1',
                              'CDH1',
                              'FOXA1',
                              'ETV5',
                              'BRAF',
                              'NCOR1',
                              'TP53',
                              'BRCA2',
                              'PRKDC',
                              'ZFHX3',
                              'PTEN',
                              'MYC',
                              'RB1',
                              'TP53',
                              'AR')) +
  ylim(-0.27,0.7) +
  geom_errorbar(inherit.aes = TRUE, stat = "identity",
                aes(ymin = observed.delta.minus.expected.delta.ci.lower, ymax = observed.delta.minus.expected.delta.ci.upper)) +
  geom_bar(fill = "black", stat = "identity", color = "white", width = 1) +
  coord_flip() +
  geom_bar(fill = "black", stat = "identity", color = "white", width = 1) +
  coord_flip() +
  theme(
    panel.background = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold",
                                size = 26),
    axis.line = element_line(),
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    axis.text.y = element_blank(),
    legend.position = c(0,0),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(0.1,0.5,0,0.1), "cm")
  )
DELTA.DELTA.PROPORTIONS

##### SCATTER PLOT #####

Figure2.Scatter <- Figure2 %>%
  ggplot(aes(x = proportion.mets, y = proportion.localized)) +
  geom_point(aes(size = log10.qval, color = type), alpha = 0.5) +
  scale_size_continuous(expression(bold(paste("\n-log"[1*0]~"q-value")))) + 
  scale_color_manual(expression(bold(paste("Mutation Type"))),
                     values = Mutation_Colours) +
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  geom_text_repel(data = Figure2 %>%
                    filter(gene.id == "SPOP" | 
                             gene.id == "AR" & type == "CNA" | 
                             gene.id == "ZNRF3" & type == "CNA" |
                             gene.id == "TP53" & type == "SNV" |
                             gene.id == "MYC" & type == "CNA" |
                             gene.id == "CCND1" & type == "CNA"),
                  aes(label = gene.id),
                  nudge_y = 0.1,
                  size = 10,
                  fontface = "italic",
                  segment.colour = "dodgerblue4") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75)) +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  ylab("Proportion (Localized)") +
  xlab("Proportion (Metastatic)") +
  expand_limits(x = c(0,0.8), y = c(0,0.8)) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA,
                                size = 2,
                                color = "black"),
    axis.title.x = element_text(size = 36,
                                face = "bold",
                                color = "black",
                                margin = margin(t = 10, b= 0, l = 0, r = 0)),
    axis.title.y = element_text(size = 36,
                                face = "bold",
                                color = "black",
                                margin = margin(t = 0, b= 0, l = 0, r = 10)),
    axis.text.x = element_text(size = 24,
                               face = "bold",
                               color = "black",
                               margin = margin(t = 10, b = 0, l = 0, r = 0)),
    axis.text.y = element_text(size = 24,
                               face = "bold",
                               color = "black",
                               margin = margin(t = 0, b = 0, l = 0, r = 10)),
    legend.title = element_text(size = 22,
                                color = "black",
                                face = "bold",
                                hjust = 1),
    legend.text = element_text(size = 18,
                               color = "black",
                               face = "bold"),
    legend.key = element_rect(size = 15,
                              fill = "NA"),
    legend.position = c(0.2, 0.65)
  )

Figure2.Scatter
#### COVARIATE BAR ####
DELTA.PROPORTIONS.TYPE <- ggplot(data = Figure2, aes(x=reorder(Abb.Code, adjusted.delta.proportion))) +
  coord_flip(ylim = c(0,0.01)) +
  geom_bar(aes(fill = type), color = "white", width = 1) +
  scale_fill_manual(values = Mutation_Colours) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  )

DELTA.PROPORTIONS.TYPE
#### LEGEND ####

DELTA.PROPORTIONS.LEGEND <- ggplot(data = Figure2, aes(x = type)) +
  coord_flip(ylim = c(0,0.01)) +
  geom_bar(aes(fill = type)) +
  scale_fill_manual(values = Mutation_Colours) +
  scale_y_continuous(expand = c(0,0))+
  theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.background = element_blank()
  )

DELTA.PROPORTIONS.LEGEND

#### BARPLOT OF Q-VALUES ####

ADJUSTED.DELTA.PROPORTIONS.Q_VALUES <- ggplot(data = Figure2, aes(x=reorder(Abb.Code, adjusted.delta.proportion), y = log10.qval)) +
  xlab("Gene\n") +
  ylab(expression(bold(paste("\n-log"[1*0]~"q-value")))) +
  labs(fill = "MUTATION TYPE") +
  scale_x_discrete(labels = c('CHD1',
                              'RB1',
                              'NKX3-1',
                              'TTN',
                              'CDK12',
                              'KDM6A',
                              'chr3:125 Mbp',
                              'MSH2',
                              'MED12',
                              'PCDH15',
                              'chr3:129 Mbp',
                              'chr6:97 Mbp',
                              'SPOP',
                              'CTNNB1',
                              'KDM6A',
                              'chr7:145019604',
                              'BRCA1',
                              'PRKDC',
                              'chr22:31826396',
                              'PIK3CA',
                              'chr16:64910033',
                              'chr2:131394079',
                              'chr13:58055742',
                              'MLH1',
                              'SYNE1',
                              'BRAF',
                              'chr4:148 Mbp',
                              'BRCA1',
                              'GNAS',
                              'ASXL2',
                              'chr4:39684557',
                              'MYC',
                              'NCOR2',
                              'BRAF',
                              'ERCC2',
                              'MSH2',
                              'NCOR1',
                              'AKT1',
                              'MSH6',
                              'ATM',
                              'IDH1',
                              'AXL',
                              'PTEN',
                              'ETV5',
                              'MSH2',
                              'MED12',
                              'CDKN2A',
                              'ERCC2',
                              'CDKN2A',
                              'PRB4',
                              'HRAS',
                              'ASXL2',
                              'TPTE',
                              'RB1',
                              'ADAMTS20',
                              'PRKDC',
                              'ATM',
                              'ZFHX3',
                              'RYR2',
                              'FAT3',
                              'MSH6',
                              'HDAC4',
                              'QRICH2',
                              'NCOR1',
                              'SPTA1',
                              'CHD1',
                              'ASXL2',
                              'ETV4',
                              'CSMD3',
                              'KDM6A',
                              'ERG',
                              'chr17:25 Mbp',
                              'FOXA1',
                              'KMT2C',
                              'APC',
                              'PTEN',
                              'ZBTB16',
                              'BRCA2',
                              'AHNAK2',
                              'ETV1',
                              'KMT2D',
                              'ATM',
                              'ETV1',
                              'FOXA1',
                              'HDAC4',
                              'TP53',
                              'ZNRF3',
                              'chr10:89 Mbp',
                              'CDKN2A',
                              'APC',
                              'MUC16',
                              'ZNRF3',
                              'AKT1',
                              'CTNNB1',
                              'chr21:42 Mbp',
                              'CDH1',
                              'AR',
                              'chr7:61 Mbp',
                              'CCND1',
                              'ETV5',
                              'NKX3-1',
                              'FOXA1',
                              'BRAF',
                              'NCOR1',
                              'BRCA2',
                              'PRKDC',
                              'PTEN',
                              'TP53',
                              'ZFHX3',
                              'RB1',
                              'MYC',
                              'TP53',
                              'AR')) +
  scale_y_continuous(breaks = c(0, 2.5, 5), expand = c(0,0)) +
  geom_bar(fill = "black", stat = "identity", color = "white", width = 1) +
  geom_hline(yintercept = 1.3, color = "red") +
  coord_flip(ylim=c(0,5.1)) +
  theme(
    panel.background = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold",
                                size = 28),
    axis.line = element_line(),
    panel.border = element_rect(fill = NA),
    axis.text.x = element_text(face = "bold", size = 24, color = "black"),
    axis.text.y = element_blank(),
    legend.position = c(0.65,0.1),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(0.1,0.6,0,.6), "cm")
  )
ADJUSTED.DELTA.PROPORTIONS.Q_VALUES
#### MERGE PLOTS TO FINAL FIGURE ####

DELTA.DELTA.PROPORTIONS.WITH.INSET <- ggdraw() +
  draw_plot(DELTA.DELTA.PROPORTIONS) + 
  draw_plot(DELTA.PROPORTIONS.LEGEND, x = 0.62, y = -0.3, width = 1, height = 1, scale = 0.075, hjust = 0.5) + 
  draw_label("MUTATION CLASS", fontface = "bold", size = 18, x = 0.6, y = 0.25, hjust = 0) +
  draw_label("SV", fontface = "bold", size = 12, x = 0.66, y = 0.229, hjust = 0, vjust = 1) +
  draw_label("SNV", fontface = "bold", size = 12, x = 0.66, y = 0.213, hjust = 0, vjust = 1) +
  draw_label("ncSNV", fontface = "bold", size = 12, x = 0.66, y = 0.199, hjust = 0, vjust = 1) +
  draw_label("CNA", fontface = "bold", size = 12, x = 0.66, y = 0.183, hjust = 0, vjust = 1)

DELTA.DELTA.PROPORTIONS.WITH.INSET

ZNRF3_Figure2_Merged <- plot_grid(ADJUSTED.DELTA.PROPORTIONS, DELTA.PROPORTIONS.TYPE, ADJUSTED.DELTA.PROPORTIONS.Q_VALUES, align = "h", ncol = 3, rel_widths = c(27,2,16))
ZNRF3_Figure2_Merged <- ZNRF3_Figure2_Merged + draw_plot(DELTA.PROPORTIONS.LEGEND, x = 0.35, y = -0.35, width = 1, height = 1, scale = 0.06, hjust = 0.5) + 
  draw_label("MUTATION CLASS", fontface = "bold", size = 24, x = 0.305, y = 0.19, hjust = 0) +
  draw_label("SV", fontface = "bold", size = 14, x = 0.39, y = 0.172, hjust = 0, vjust = 1) +
  draw_label("SNV", fontface = "bold", size = 14, x = 0.39, y = 0.16, hjust = 0, vjust = 1) +
  draw_label("ncSNV", fontface = "bold", size = 14, x = 0.39, y = 0.148, hjust = 0, vjust = 1) +
  draw_label("CNA", fontface = "bold", size = 14, x = 0.39, y = 0.135, hjust = 0, vjust = 1)
ZNRF3_Figure2_Merged
