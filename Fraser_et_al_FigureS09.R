# FRASER ET AL - FIGURE S09 #

# LOAD LIBRARIES ####
library(tidyverse)
library(survminer)
library(easyGgplot2)

# LOAD DATA - CPCG ####

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCGENE.OUTCOME.ADJUSTED.PGA.rds")

CPCG.ZNRF3.CNA.AdjPGA.SUMMARY <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  group_by(ZNRF3.CNA) %>%
  summarise(MEAN.PGA= mean(pga.no.chr22), 
            SD.PGA = sd(pga.no.chr22),
            N.PGA = n()) %>%
  mutate(SE.PGA = SD.PGA/sqrt(N.PGA),
         lower.ci.PGA = MEAN.PGA - qt(1 - (0.05 / 2), N.PGA - 1) * SE.PGA,
         upper.ci.PGA = MEAN.PGA + qt(1 - (0.05 / 2), N.PGA - 1) * SE.PGA)
CPCG.ZNRF3.CNA.AdjPGA.SUMMARY

## BOX AND WHISKERS PLOT - CPCG ####

ZNRF3.ADJ.PGA.BOXPLOT <- ggplot2.boxplot(
  data = CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA,
  xName = 'ZNRF3.CNA',
  yName = 'adjusted.pga',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ytitle = "Percentage Genome Alteration (PGA)",
  ytitleFont = c(24, "bold","black"),
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA)

ZNRF3.ADJ.PGA.BOXPLOT

#Add p-value as annotation
ZNRF3.ADJ.PGA.BOXPLOT <- ZNRF3.ADJ.PGA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 41,
  label = "P == 4.57 %*% 10^-3",
  size = 10,
  fontface = "bold",
  parse=TRUE
) +
  labs(y = expression(bold('Adjusted PGA')))

ZNRF3.ADJ.PGA.BOXPLOT <- ZNRF3.ADJ.PGA.BOXPLOT + 
  theme(axis.title.y = element_text(size = 36,
                                    face = "bold",
                                    color = "black",
                                    margin = margin(t=0,b=0,r=10,l=0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 30,
                                   color = "black",
                                   margin = margin(t=0,b=0,l=0,r=10)),
        axis.text.x = element_text(face = "bold",
                                   size = 30,
                                   color = "black",
                                   margin = margin(t=10,b=0,l=0,r=0)),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    size = 2),
        aspect.ratio = 1) +
  stat_boxplot(geom='errorbar',
               width = 0.5)
ZNRF3.ADJ.PGA.BOXPLOT <- ZNRF3.ADJ.PGA.BOXPLOT +
  scale_x_discrete(labels = c("1" = "ZNRF3 Loss", "0" = "ZNRF3 Neutral"))
#View box plot
ZNRF3.ADJ.PGA.BOXPLOT



#
# LOAD AND TIDY DATA - TCGA ####

saveRDS(TCGA_CNA_RNA_PGA_PFS_, "/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA_CNA_RNA_PGA_PFS_.rds")
TCGA_CNA_RNA_PGA_PFS_ <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA_CNA_RNA_PGA_PFS_.rds")
View(TCGA_CNA_RNA_PGA_PFS_)

TCGA <- TCGA_CNA_RNA_PGA_PFS_

TCGA <- TCGA %>%
  filter(Adjusted.PGA != "#N/A") %>%
  mutate(median.adjusted.pga = ifelse(Adjusted.PGA > median(Adjusted.PGA), 1, 0))

TCGA$Adjusted.PGA <- as.numeric(TCGA$Adjusted.PGA)
TCGA$Adjusted.PGA <- (TCGA$Adjusted.PGA * 100)

View(TCGA)
TCGA.PFS.MVCOX <- coxph(Surv(time, bcr) ~ ZNRF3.CNA + median.adjusted.pga, data = TCGA)
summary(TCGA.PFS.MVCOX)

TCGA.ZNRF3.CNA.AdjPGA.SUMMARY <- TCGA %>%
  group_by(ZNRF3.CNA) %>%
  summarise(MEAN.PGA= mean(Adjusted.PGA), 
            SD.PGA = sd(Adjusted.PGA),
            N.PGA = n()) %>%
  mutate(SE.PGA = SD.PGA/sqrt(N.PGA),
         lower.ci.PGA = MEAN.PGA - qt(1 - (0.05 / 2), N.PGA - 1) * SE.PGA,
         upper.ci.PGA = MEAN.PGA + qt(1 - (0.05 / 2), N.PGA - 1) * SE.PGA)
TCGA.ZNRF3.CNA.AdjPGA.SUMMARY

wilcox.test(Adjusted.PGA ~ ZNRF3.CNA, data = TCGA)
## BOX AND WHISKERS PLOT - TCGA ####
TCGA.ZNRF3.ADJ.PGA.BOXPLOT <- ggplot2.boxplot(
  data = TCGA,
  xName = 'ZNRF3.CNA',
  yName = 'Adjusted.PGA',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA)

#Add p-value as annotation
TCGA.ZNRF3.ADJ.PGA.BOXPLOT <- TCGA.ZNRF3.ADJ.PGA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 95,
  label = "P == 2.11 %*% 10^-10",
  size = 10,
  fontface = "bold",
  parse=TRUE
) +
  labs(y = expression(bold('Adjusted PGA')))

TCGA.ZNRF3.ADJ.PGA.BOXPLOT <- TCGA.ZNRF3.ADJ.PGA.BOXPLOT + 
  theme(axis.title.y = element_text(size = 36,
                                    face = "bold",
                                    color = "black",
                                    margin = margin(t=0,b=0,r=10,l=0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 30,
                                   color = "black",
                                   margin = margin(t=0,b=0,l=0,r=10)),
        axis.text.x = element_text(face = "bold",
                                   size = 30,
                                   color = "black",
                                   margin = margin(t=10,b=0,l=0,r=0)),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    size = 2),
        aspect.ratio = 1) +
  stat_boxplot(geom='errorbar',
               width = 0.5)
TCGA.ZNRF3.ADJ.PGA.BOXPLOT <- TCGA.ZNRF3.ADJ.PGA.BOXPLOT +
  scale_x_discrete(labels = c("1" = "ZNRF3 Loss", "0" = "ZNRF3 Neutral"))
TCGA.ZNRF3.ADJ.PGA.BOXPLOT

#
# LOAD AND TIDY DATA - MSKCC ####

Taylor <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Taylor.ZNRF3.CNA.Adjusted.PGA.rds")

Taylor <- Taylor %>%
  filter(Adjusted.PGA != "#N/A") %>%
  mutate(median.adjusted.pga = ifelse(Adjusted.PGA > median(Adjusted.PGA), 1, 0))

Taylor$Adjusted.PGA <- as.numeric(Taylor$Adjusted.PGA * 100)

View(Taylor)

Taylor.ZNRF3.CNA.AdjPGA.SUMMARY <- Taylor %>%
  group_by(ZNRF3.CNA) %>%
  summarise(MEAN.PGA= mean(Adjusted.PGA), 
            SD.PGA = sd(Adjusted.PGA),
            N.PGA = n()) %>%
  mutate(SE.PGA = SD.PGA/sqrt(N.PGA),
         lower.ci.PGA = MEAN.PGA - qt(1 - (0.05 / 2), N.PGA - 1) * SE.PGA,
         upper.ci.PGA = MEAN.PGA + qt(1 - (0.05 / 2), N.PGA - 1) * SE.PGA)
Taylor.ZNRF3.CNA.AdjPGA.SUMMARY

wilcox.test(as.numeric(Adjusted.PGA) ~ ZNRF3.CNA, Taylor)

## BOX AND WHISKERS PLOT - MSKCC ####
library(easyGgplot2)
Taylor.ZNRF3.ADJ.PGA.BOXPLOT <- ggplot2.boxplot(
  data = Taylor,
  xName = 'ZNRF3.CNA',
  yName = 'Adjusted.PGA',
  groupName = 'ZNRF3.CNA',
  groupColors =c ('#FF0000', '#104E8B'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ytitle = "Percentage Genome Alteration (PGA)",
  ytitleFont = c(24, "bold","black"),
  xTickLabelFont = c(24,"bold", "black"),
  outlier.shape = NA)

#Add p-value as annotation
Taylor.ZNRF3.ADJ.PGA.BOXPLOT <- Taylor.ZNRF3.ADJ.PGA.BOXPLOT + annotate(
  "text",
  x = 2,
  y = 85,
  label = "P == 6.4 %*% 10^-4",
  size = 10,
  fontface = "bold",
  parse=TRUE
) +
  labs(y = expression(bold('Adjusted PGA')))

Taylor.ZNRF3.ADJ.PGA.BOXPLOT <- Taylor.ZNRF3.ADJ.PGA.BOXPLOT + 
  theme(axis.title.y = element_text(size = 36,
                                    face = "bold",
                                    color = "black",
                                    margin = margin(t=0,b=0,r=10,l=0)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 30,
                                   color = "black",
                                   margin = margin(t=0,b=0,l=0,r=10)),
        axis.text.x = element_text(face = "bold",
                                   size = 30,
                                   color = "black",
                                   margin = margin(t=10,b=0,l=0,r=0)),
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,
                                    size = 2),
        aspect.ratio = 1) +
  stat_boxplot(geom='errorbar',
               width = 0.5)
Taylor.ZNRF3.ADJ.PGA.BOXPLOT <- Taylor.ZNRF3.ADJ.PGA.BOXPLOT +
  scale_x_discrete(labels = c("1" = "ZNRF3 Loss", "0" = "ZNRF3 Neutral"))
#View box plot
Taylor.ZNRF3.ADJ.PGA.BOXPLOT
