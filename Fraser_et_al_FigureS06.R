# FRASER ET AL - FIGURE S06 #

# LOAD LIBRARIES ####
library(tidyverse)

# IMPORT DATA - TCGA ####

TCGA_ZNRF3_RNA_ISUP <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA_ZNRF3_RNA_ISUP.rds")

## CALCULATE MEANS, 95% CI FOR GROUPS BASED ON ZNRF3 CNA STATUS - TCGA ####
TCGA_ZNRF3_RNA_ISUP %>%
  group_by(isup) %>%
  summarise(mean = mean(ZNRF3.RNA.ZScore),
            std = sqrt(var(ZNRF3.RNA.ZScore)),
            lower = mean(ZNRF3.RNA.ZScore) - qnorm(.975)*std/sqrt(n()),
            upper = mean(ZNRF3.RNA.ZScore) + qnorm(.975)*std/sqrt(n()))



## BINARIZE ISUP INTO GRADE 1-2 AND GRADE 3 AND ABOVE - TCGA ####
TCGA <- TCGA_ZNRF3_RNA_ISUP %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 1, 2)))

## STATISTICAL TESTING OF GROUP MEANS - TCGA ####

isup3 = as.factor(TCGA$isup)
isup3.means <- tapply(TCGA$ZNRF3.RNA.ZScore, isup3, mean)
isup3.means
isup.nn <- tapply(TCGA$ZNRF3.RNA.ZScore, TCGA$isup, length); isup.sd <- 
  tapply(TCGA$ZNRF3.RNA.ZScore, TCGA$isup, sd); isup.nn; isup.sd
tcga.isup.anova <- aov(TCGA$ZNRF3.RNA.ZScore ~ TCGA$isup)
summary(tcga.isup.anova, intercept = F)
tcga.isup.anova <- aov(TCGA$ZNRF3.RNA.ZScore ~ isup3, data = TCGA)
pairwise.t.test(TCGA$ZNRF3.RNA.ZScore, TCGA$isup, p.adjust="BH")
TukeyHSD(tcga.isup.anova, conf.level = 0.95)

## CREATE BASE PLOT - TCGA ####

TCGA_ZNRF3_RNA_ISUP <- TCGA_ZNRF3_RNA_ISUP %>%
  drop_na(ZNRF3.CNA)
View(TCGA_ZNRF3_RNA_ISUP)

TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- ggplot2.boxplot(
  data = TCGA,
  xName = 'isup',
  yName = 'ZNRF3.RNA.ZScore',
  groupName = 'isup',
  groupColors =c ('#FF0000', 'dodgerblue2', 'pink'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ylim = c(-2.5, 3.5),
  ytitle = "ZNRF3 RNA Abundance (Z-Score)",
  outlier.shape = NA)

#Add p-value as annotation
TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + annotate(
  "text",
  x = 1.05,
  y = -2.5,
  label = "P == 2.90 %*% 10^-6", # Add p-value from one-way ANOVA above
  size = 10,
  fontface = "bold",
  parse=TRUE
)

TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT +
  theme(axis.title.y = element_text(size = 30,
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

TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + 
  scale_x_discrete(labels = c("1" = "Grade 1", "2" = "Grade 2", "3" = "\U2265Grade 3"))

#View box plot
TCGA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT

#
#
# IMPORT DATA - CPCG ####

CPCGENE <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCG.ZNRF3.RNA.ISUP.rds")
## CALCULATE MEANS, 95% CI FOR GROUPS BASED ON ZNRF3 CNA STATUS ####

CPCGENE %>%
  group_by(as.factor(isup)) %>%
  summarise(mean = mean(ZNRF3.RNA.ZScore),
            std = sqrt(var(ZNRF3.RNA.ZScore)),
            lower = mean(ZNRF3.RNA.ZScore) - qnorm(.975)*std/sqrt(n()),
            upper = mean(ZNRF3.RNA.ZScore) + qnorm(.975)*std/sqrt(n()))

## BINARIZE ISUP INTO GRADE 1-2 AND GRADE 3 AND ABOVE - CPCG ####
CPCGENE <- CPCGENE %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 1, 2)))

## STATISTICAL TESTING OF GROUP MEANS - CPCG ####

isup3 <- as.factor(CPCGENE$isup)
isup3.means <- tapply(CPCGENE$ZNRF3.RNA.ZScore, isup3, mean)
isup3.means
isup.nn <- tapply(CPCGENE$ZNRF3.RNA.ZScore, CPCGENE$isup, length); isup.sd <- 
  tapply(CPCGENE$ZNRF3.RNA.ZScore, CPCGENE$isup, sd); isup.nn; isup.sd
cpcg.isup.anova <- aov(CPCGENE$ZNRF3.RNA.ZScore~CPCGENE$isup)
summary(cpcg.isup.anova, intercept = F)
cpcg.isup.anova <- aov(CPCGENE$ZNRF3.RNA.ZScore ~ isup3, data = CPCGENE)
pairwise.t.test(CPCGENE$ZNRF3.RNA.ZScore, CPCGENE$isup, p.adjust="BH")
TukeyHSD(cpcg.isup.anova, conf.level = 0.95)


## CREATE BASE PLOT - CPCG ####

CPCGENE <- CPCGENE %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 1, 2)))

CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- ggplot2.boxplot(
  data = CPCGENE,
  xName = 'isup',
  yName = 'ZNRF3.RNA.ZScore',
  groupName = 'isup',
  groupColors =c ('#FF0000', 'dodgerblue2', 'pink'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ylim = c(-2.5, 3.5),
  ytitle = "ZNRF3 RNA Abundance (Z-Score)",
  outlier.shape = NA)

#Add p-value as annotation
CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + annotate(
  "text",
  x = 1.05,
  y = -2.5,
  label = "P == 2.43 %*% 10^-3", # Add p-value from one-way ANOVA above
  size = 10,
  fontface = "bold",
  parse=TRUE
)

CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT +
  theme(axis.title.y = element_text(size = 30,
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

CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + 
  scale_x_discrete(labels = c("1" = "Grade 1", "2" = "Grade 2", "3" = "\U2265Grade 3"))

#View box plot
CPCG.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT

#
#
# IMPORT DATA - EOPC ####

DKFZ <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/EOPC.ZNRF3.RNA.ISUP.rds")

## CALCULATE MEANS, 95% CI FOR GROUPS BASED ON ZNRF3 CNA STATUS - EOPC ####
DKFZ %>%
  group_by(isup) %>%
  summarise(mean = mean(ZNRF3.RNA.ZScore),
            std = sqrt(var(ZNRF3.RNA.ZScore)),
            lower = mean(ZNRF3.RNA.ZScore) - qnorm(.975)*std/sqrt(n()),
            upper = mean(ZNRF3.RNA.ZScore) + qnorm(.975)*std/sqrt(n()))

## BINARIZE ISUP INTO GRADE 1-2 AND GRADE 3 AND ABOVE - EOPC ####
DKFZ <- DKFZ %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 1, 2)))

## STATISTICAL TESTING OF GROUP MEANS - EOPC ####
isup2 <- as.factor(DKFZ$isup)
isup2.means <- tapply(DKFZ$ZNRF3.RNA.ZScore, isup2, mean)
isup2.means
isup.nn <- tapply(DKFZ$ZNRF3.RNA.ZScore, isup2, length); isup.sd <- 
  tapply(DKFZ$ZNRF3.RNA.ZScore, isup2, sd); isup.nn; isup.sd
dkfz.isup.anova <- aov(DKFZ$ZNRF3.RNA.ZScore ~ isup2)
summary(dkfz.isup.anova, intercept = F)
dkfz.isup.anova <- aov(DKFZ$ZNRF3.RNA.ZScore ~ isup2, data = DKFZ)
pairwise.t.test(DKFZ$ZNRF3.RNA.ZScore, isup2, p.adjust="BH")
TukeyHSD(dkfz.isup.anova, conf.level = 0.95)

# CREATE BASE PLOT - EOPC ####

DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- ggplot2.boxplot(
  data = DKFZ,
  xName = 'isup',
  yName = 'ZNRF3.RNA.ZScore',
  groupName = 'isup',
  groupColors =c ('#FF0000', 'dodgerblue2', 'pink'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ylim = c(-2.5, 3.5),
  ytitle = "ZNRF3 RNA Abundance (Z-Score)",
  outlier.shape = NA)

#Add p-value as annotation
DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + annotate(
  "text",
  x = 1.05,
  y = -2.5,
  label = "P == 3.57 %*% 10^-6", # Add p-value from one-way ANOVA above
  size = 10,
  fontface = "bold",
  parse=TRUE
)


DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + 
  theme(axis.title.y = element_text(size = 30,
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

DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + 
  scale_x_discrete(labels = c("1" = "Grade 1", "2" = "Grade 2", "3" = "\U2265Grade 3"))

#View box plot
DKFZ.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT

#
# IMPORT DATA - LTRI ####

ZLOTTA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/LTRI.ZNRF3.RNA.ISUP.rds")

## BINARIZE ISUP INTO GRADE 1-2 AND GRADE 3 AND ABOVE - LTRI ####
mean.znrf3 <- mean(ZLOTTA$znrf3.rna)
sd.znrf3 <- sd(ZLOTTA$znrf3.rna)

ZLOTTA <- ZLOTTA %>%
  drop_na(znrf3.rna) %>%
  group_by(patient.number) %>%
  mutate(ZNRF3.ZSCORE = ((znrf3.rna - mean.znrf3)/sd.znrf3))

## STATISTICAL TESTING OF GROUP MEANS - LTRI ####
isup <- as.factor(ZLOTTA$isup)
isup.means <- tapply(ZLOTTA$ZNRF3.ZSCORE, ZLOTTA$isup, mean)
isup.means
isup.nn <- tapply(ZLOTTA$ZNRF3.ZSCORE, ZLOTTA$isup, length); isup.sd <- tapply(ZLOTTA$ZNRF3.ZSCORE, ZLOTTA$isup, sd); isup.nn; isup.sd
zlotta.isup.anova <- aov(ZLOTTA$ZNRF3.ZSCORE~isup)
summary(zlotta.isup.anova, intercept = F)
zlotta.isup.anova <- aov(ZLOTTA$ZNRF3.ZSCORE ~ as.factor(ZLOTTA$isup), data = ZLOTTA)
pairwise.t.test(ZLOTTA$ZNRF3.ZSCORE, isup, p.adjust="BH")
TukeyHSD(zlotta.isup.anova, conf.level = 0.95)

# CREATE BASE PLOT - LTRI ####
ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- ggplot2.boxplot(
  data = ZLOTTA,
  xName = 'isup',
  yName = 'ZNRF3.ZSCORE',
  groupName = 'isup',
  groupColors =c ('#FF0000', 'dodgerblue', 'pink'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ylim = c(-2.5, 3.5),
  ytitle ="ZNRF3 RNA Abundance (Z-Score)",
  outlier.shape = NA)

#Add p-value as annotation
ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + annotate(
  "text",
  x = 1.05,
  y = -2.5,
  label = "P == 2.04 %*% 10^-2", # Add p-value from one-way ANOVA above
  size = 10,
  fontface = "bold",
  parse=TRUE
)

ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT  <- ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT  + 
  theme(axis.title.y = element_text(size = 30,
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
ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT  <- ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT  +
  scale_x_discrete(labels = c("1" = "Grade 2", "0" = "Grade 1", "2" = "\U2265 Grade 3"))
#View box plot
ZLOTTA.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT 

#
# IMPORT DATA - MSKCC ####
MSKCC <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/MSKCC.ZNRF3.RNA.ISUP.rds")

MSKCC <- MSKCC %>%
  drop_na(isup)
View(MSKCC)

## STATISTICAL TESTING OF GROUP MEANS - MSKCC ####

isup2 = as.factor(MSKCC$isup)
isup2.means <- tapply(MSKCC$ZNRF3.RNA.ZScore, isup2, mean)
isup2.means
isup.nn <- tapply(MSKCC$ZNRF3.RNA.ZScore, isup2, length); isup.sd <- 
  tapply(MSKCC$ZNRF3.RNA.ZScore, isup2, sd); isup.nn; isup.sd
mskcc.isup.anova <- aov(MSKCC$ZNRF3.RNA.ZScore ~ isup2)
summary(mskcc.isup.anova, intercept = F)
mskcc.isup.anova <- aov(MSKCC$ZNRF3.RNA.ZScore ~ isup2, data = MSKCC)
pairwise.t.test(MSKCC$ZNRF3.RNA.ZScore, isup2, p.adjust="BH")
TukeyHSD(mskcc.isup.anova, conf.level = 0.95)

## CALCULATE MEANS, 95% CI FOR GROUPS BASED ON ISUP GROUP - MSKCC ####
MSKCC %>%
  group_by(isup) %>%
  summarise(mean = mean(ZNRF3.RNA.ZScore),
            std = sqrt(var(ZNRF3.RNA.ZScore)),
            lower = mean(ZNRF3.RNA.ZScore) - qnorm(.975)*std/sqrt(n()),
            upper = mean(ZNRF3.RNA.ZScore) + qnorm(.975)*std/sqrt(n()))


MSKCC <- MSKCC %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 1,2)))
# CREATE BASE PLOT - MSKCC ####
MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- ggplot2.boxplot(
  data = MSKCC,
  xName = 'isup',
  yName = 'ZNRF3.RNA.ZScore',
  groupName = 'isup',
  groupColors =c ('#FF0000', 'dodgerblue2', 'pink'),
  addDot = TRUE,
  dotSize = 1,
  dotPosition = "jitter",
  xShowTitle = FALSE,
  ylim = c(-2.5, 3.5),
  ytitle = "ZNRF3 RNA Abundance (Z-Score)",
  outlier.shape = NA)

MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT


#Add p-value as annotation
MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT + annotate(
  "text",
  x = 1.05,
  y = -2.5,
  label = "P == 0.821", # Add p-value from one-way ANOVA above
  size = 10,
  fontface = "bold",
  parse=TRUE
)

MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT <- MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT +
  theme(axis.title.y = element_text(size = 30,
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
MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT  <- MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT  +
  scale_x_discrete(labels = c("1" = "Grade 1",
                              "2" = "Grade 2",
                              "3" = "\U2265Grade 3"))

MSKCC.ZNRF3.RNA.ZSCORE.ISUP.BOXPLOT
