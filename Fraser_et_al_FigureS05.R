# FRASER ET AL - FIGURE S05 #

# Libraries ####

library(tidyverse)

# Directories ####

input_dir <- file.path("~/OneDrive/Work/Projects/ZNRF3/")

#
# Import data - CPCGENE ####
saveRDS(ZNRF3.CNA.RNA.TABLE, "/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCG.ZNRF3.CNA.RNA.TABLE.rds")
CPCG.ZNRF3.CNA.RNA.TABLE <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCG.ZNRF3.CNA.RNA.TABLE.rds")
# DENSITY PLOT - CPCG #####
CPCG.ZNRF3.CNA.RNA.TABLE %>%
  ggplot(aes(x = ZNRF3.RNA)) +
  geom_density(size = 2) +
  labs(x = expression(bold('ZNRF3 RNA Abundance ('*log["2"]*')')),
       y = "Density") +
  xlim(6,8) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     limits = c(0,2.1)) +
  theme(
    axis.text.x = element_text(size = 36,
                               face = "bold",
                               color = "black"),
    axis.text.y = element_text(size = 50,
                               face = "bold",
                               color = "black"),
    axis.title.x = element_text(size = 44,
                                face = "bold",
                                color = "black",
                                margin = margin(t=10,b=0,r=0,l=0)),
    axis.title.y = element_text(size = 50,
                                face = "bold",
                                color = "black",
                                margin = margin(t=0,b=0,r=10,l=0)),
    panel.border = element_rect(fill = "NA",
                                size = 2),
    panel.background = element_blank(),
    legend.title = element_text(size = 18,
                                face = "bold",
                                color = "black"),
    legend.text = element_text(size = 12,
                               face = "italic",
                               color = "black"),
    aspect.ratio = 1
  )

## BOX AND WHISKERS PLOT - CPCG ####

ZNRF3.RNA.BOXPLOT <- ggplot2.boxplot(
  data = CPCG.ZNRF3.CNA.RNA.TABLE,
  xName = 'ZNRF3.CNA',
  yName = 'ZNRF3.RNA',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ytitle = "ZNRF3 RNA Abundance\n",
  ytitleFont = c(24, "bold","black"),
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA
) + theme(
  plot.background = element_blank()
)

#Add p-value as annotation
wilcox.test(ZNRF3.RNA ~ ZNRF3.CNA, data = CPCG.ZNRF3.CNA.RNA.TABLE)

ZNRF3.RNA.BOXPLOT <- ZNRF3.RNA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 7,
  label = "P == 2.11 %*% 10^-2",
  size = 7,
  fontface = "bold",
  parse=TRUE
)

ZNRF3.RNA.BOXPLOT <- ZNRF3.RNA.BOXPLOT + 
  theme(axis.text.y = element_text(face = "bold", size = 24), 
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(size = 0.5, color='black'),
        panel.border = element_rect(fill = NA,
                                    size = 2,
                                    colour = "black"),
        aspect.ratio = 1) + 
  stat_boxplot(geom='errorbar', width = 0.5,) +
  scale_x_discrete(labels = c("1" = "\nZNRF3 Loss", "0" = "\nZNRF3 Neutral"))
#View box plot
ZNRF3.RNA.BOXPLOT


#
# IMPORT DATA - TCGA ####

TCGA_ZNRF3_CNA_RNA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA_ZNRF3_CNA_RNA.rds")

# DENSITY PLOT - TCGA ####
TCGA_ZNRF3_CNA_RNA %>%
  ggplot(aes(x = ZNRF3.RNA)) +
  geom_density(size = 2) +
  labs(x = expression(bold('ZNRF3 RNA Abundance ('*log["2"]*')')),
       y = "Density") +
  xlim(6,12) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4),
                     limits = c(0,0.6)) +
  theme(
    axis.text.x = element_text(size = 36,
                               face = "bold",
                               color = "black"),
    axis.text.y = element_text(size = 50,
                               face = "bold",
                               color = "black"),
    axis.title.x = element_text(size = 44,
                                face = "bold",
                                color = "black",
                                margin = margin(t=10,b=0,r=0,l=0)),
    axis.title.y = element_text(size = 50,
                                face = "bold",
                                color = "black",
                                margin = margin(t=0,b=0,r=10,l=0)),
    panel.border = element_rect(fill = "NA",
                                size = 2),
    panel.background = element_blank(),
    legend.title = element_text(size = 18,
                                face = "bold",
                                color = "black"),
    legend.text = element_text(size = 12,
                               face = "italic",
                               color = "black"),
    aspect.ratio = 1
  )

# BOX AND WHISKERS PLOT - TCGA ####
TCGA.ZNRF3.RNA.BOXPLOT <- ggplot2.boxplot(
  data = TCGA_ZNRF3_CNA_RNA,
  xName = 'ZNRF3.CNA',
  yName = 'ZNRF3.RNA',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ytitle = "ZNRF3 RNA Abundance\n",
  ytitleFont = c(24, "bold","black"),
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA
) + theme(
  plot.background = element_blank()
)

#Add p-value as annotation
wilcox.test(ZNRF3.RNA ~ ZNRF3.CNA, data = TCGA_ZNRF3_CNA_RNA)
TCGA.ZNRF3.RNA.BOXPLOT <- TCGA.ZNRF3.RNA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 11.5,
  label = "P == 4.06 %*% 10^-6",
  size = 7,
  fontface = "bold",
  parse=TRUE
)

TCGA.ZNRF3.RNA.BOXPLOT <- TCGA.ZNRF3.RNA.BOXPLOT + 
  theme(axis.text.y = element_text(face = "bold", size = 24), 
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(size = 0.5, color='black'),
        panel.border = element_rect(fill = NA,
                                    size = 2,
                                    colour = "black"),
        aspect.ratio = 1) + 
  stat_boxplot(geom='errorbar', width = 0.5,) +
  scale_x_discrete(labels = c("1" = "\nZNRF3 Loss", "0" = "\nZNRF3 Neutral"))
#View box plot
TCGA.ZNRF3.RNA.BOXPLOT

#
#
# IMPORT DATA - EOPC ####

DKFZ_ZNRF3_CNA_RNA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/DKFZ_ZNRF3_CNA_RNA.rds")

# DENSITY PLOT - EOPC ####

DKFZ_ZNRF3_CNA_RNA %>%
  ggplot(aes(x = ZNRF3.RNA)) +
  geom_density(size = 2) +
  labs(x = expression(bold('ZNRF3 RNA Abundance ('*log["2"]*')')),
       y = "Density") +
  xlim(2,11.5) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1, 0.15, 0.2),
                     limits = c(0,0.25)) +
  theme(
    axis.text.x = element_text(size = 36,
                               face = "bold",
                               color = "black"),
    axis.text.y = element_text(size = 50,
                               face = "bold",
                               color = "black"),
    axis.title.x = element_text(size = 44,
                                face = "bold",
                                color = "black",
                                margin = margin(t=10,b=0,r=0,l=0)),
    axis.title.y = element_text(size = 50,
                                face = "bold",
                                color = "black",
                                margin = margin(t=0,b=0,r=10,l=0)),
    panel.border = element_rect(fill = "NA",
                                size = 2),
    panel.background = element_blank(),
    legend.title = element_text(size = 18,
                                face = "bold",
                                color = "black"),
    legend.text = element_text(size = 12,
                               face = "italic",
                               color = "black"),
    aspect.ratio = 1
  )

# BOX AND WHISKERS PLOT - EOPC ####

EOPC.ZNRF3.RNA.BOXPLOT <- ggplot2.boxplot(
  data = DKFZ_ZNRF3_CNA_RNA,
  xName = 'ZNRF3.CNA',
  yName = 'ZNRF3.RNA',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ytitle = "ZNRF3 RNA Abundance\n",
  ytitleFont = c(24, "bold","black"),
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA
) + theme(
  plot.background = element_blank()
)

#Add p-value as annotation
wilcox.test(ZNRF3.RNA ~ ZNRF3.CNA, data = DKFZ_ZNRF3_CNA_RNA)
EOPC.ZNRF3.RNA.BOXPLOT <- EOPC.ZNRF3.RNA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 9,
  label = "P = 0.113",
  size = 7,
  fontface = "bold",
  parse=FALSE
)

EOPC.ZNRF3.RNA.BOXPLOT <- EOPC.ZNRF3.RNA.BOXPLOT + 
  theme(axis.text.y = element_text(face = "bold", size = 24), 
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(size = 0.5, color='black'),
        panel.border = element_rect(fill = NA,
                                    size = 2,
                                    colour = "black"),
        aspect.ratio = 1) + 
  stat_boxplot(geom='errorbar', width = 0.5,) +
  scale_x_discrete(labels = c("1" = "\nZNRF3 Loss", "0" = "\nZNRF3 Neutral"))
#View box plot
EOPC.ZNRF3.RNA.BOXPLOT

#
# IMPORT DATA - MSKCC ####

MSKCC_ZNRF3_CNA_RNA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/MSKCC_ZNRF3_CNA_RNA.rds")

# DENSITY PLOT - MSKCC ####

MSKCC_ZNRF3_CNA_RNA %>%
  ggplot(aes(x = ZNRF3.RNA)) +
  geom_density(size = 2) +
  labs(x = expression(bold('ZNRF3 RNA Abundance ('*log["2"]*')')),
       y = "Density") +
  xlim(11,13) +
  scale_y_continuous(breaks = c(0, 1,2,3,4),
                     limits = c(0,4.5)) +
  theme(
    axis.text.x = element_text(size = 36,
                               face = "bold",
                               color = "black"),
    axis.text.y = element_text(size = 50,
                               face = "bold",
                               color = "black"),
    axis.title.x = element_text(size = 44,
                                face = "bold",
                                color = "black",
                                margin = margin(t=10,b=0,r=0,l=0)),
    axis.title.y = element_text(size = 50,
                                face = "bold",
                                color = "black",
                                margin = margin(t=0,b=0,r=10,l=0)),
    panel.border = element_rect(fill = "NA",
                                size = 2),
    panel.background = element_blank(),
    legend.title = element_text(size = 18,
                                face = "bold",
                                color = "black"),
    legend.text = element_text(size = 12,
                               face = "italic",
                               color = "black"),
    aspect.ratio = 1
  )

# BOX AND WHISKERS PLOT - MSKCC ####
MSKCC.ZNRF3.RNA.BOXPLOT <- ggplot2.boxplot(
  data = MSKCC_ZNRF3_CNA_RNA,
  xName = 'ZNRF3.CNA',
  yName = 'ZNRF3.RNA',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ytitle = "ZNRF3 RNA Abundance\n",
  ytitleFont = c(24, "bold","black"),
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA
) + theme(
  plot.background = element_blank()
)

#Add p-value as annotation
wilcox.test(ZNRF3.RNA ~ ZNRF3.CNA, data = MSKCC_ZNRF3_CNA_RNA)
MSKCC.ZNRF3.RNA.BOXPLOT <- MSKCC.ZNRF3.RNA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 12.1,
  label = "P == 1.13 %*% 10^-2",
  size = 7,
  fontface = "bold",
  parse=TRUE
)

MSKCC.ZNRF3.RNA.BOXPLOT <- MSKCC.ZNRF3.RNA.BOXPLOT + 
  theme(axis.text.y = element_text(face = "bold", size = 24), 
        legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(size = 0.5, color='black'),
        panel.border = element_rect(fill = NA,
                                    size = 2,
                                    colour = "black"),
        aspect.ratio = 1) + 
  stat_boxplot(geom='errorbar', width = 0.5,) +
  scale_x_discrete(labels = c("1" = "\nZNRF3 Loss", "0" = "\nZNRF3 Neutral"))
#View box plot
MSKCC.ZNRF3.RNA.BOXPLOT

