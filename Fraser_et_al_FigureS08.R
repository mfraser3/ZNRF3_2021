# FRASER ET AL - FIGURE S08 #

# LOAD LIBRARIES ####
library(tidyverse)
library(survminer)
library(survival)

# LOAD DATA - LTRI ####

ZLOTTA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/LTRI.ZNRF3.RNA.ISUP.rds")

ZLOTTA <- ZLOTTA %>%
  group_by(patient.id) %>%
  mutate(time_to_bcr_months = time_to_bcr * 12)

# COX PROPORTIONAL HAZARDS MODELS - LTRI ####
ZLOTTA.ZNRF3.RNA.BCR.COX <- coxph(Surv(time_to_bcr_months, bcr.bin)~
                                    znrf3.rna.status, data = ZLOTTA)
summary(ZLOTTA.ZNRF3.RNA.BCR.COX)

ZLOTTA.ZNRF3.RNA.METS.COX <- coxph(Surv(time_to_mets, mets)~
                                     znrf3.rna.status, data = ZLOTTA)
summary(ZLOTTA.ZNRF3.RNA.METS.COX)


#### KM curve of biochemical relapse in ZLOTTA patients, stratified by ZNRF3 RNA abundance (median dichotomized) - LTRI ####

fit <- survfit(Surv(time_to_bcr_months, bcr.bin) ~ znrf3.rna.status, data = ZLOTTA)
ZLOTTA.ZNRF3.RNA.BCR.KM <- ggsurvplot(
  fit,
  data = ZLOTTA,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time Post-Treatment (Months)",
  ylab = "BCR-Free Survival",
  xlim = c(0,102),
  ylim = c(0,1.01),
  break.time.by = 24,
  axes.offset = FALSE,
  font.legend =
    c(16),
  font.x =
    c(26, "bold"),
  font.y =
    c(26,"bold"),
  legend =
    c(.18, 0.35),
  legend.labs =
    c("High", "Low"),
  legend.title = "ZNRF3 RNA",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 6
)

ZLOTTA.ZNRF3.RNA.BCR.KM$table <- ZLOTTA.ZNRF3.RNA.BCR.KM$table +
  theme(
    plot.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,1,1,1,"cm")
  )

ZLOTTA.ZNRF3.RNA.BCR.KM$plot <- ZLOTTA.ZNRF3.RNA.BCR.KM$plot +
  ggplot2::annotate("text",
                    x=0.5, 
                    y=0.18,
                    label = "HR = 4.85 (2.01, 11.7)", 
                    hjust=0, 
                    size=10)
ZLOTTA.ZNRF3.RNA.BCR.KM$plot <- ZLOTTA.ZNRF3.RNA.BCR.KM$plot +
  ggplot2::annotate("text",
                    x= 0.5, 
                    y=0.08,
                    label = "P: 4.56 %*% 10^-4",
                    hjust=0,
                    size=10,
                    parse = TRUE)

ZLOTTA.ZNRF3.RNA.BCR.KM$plot <- ZLOTTA.ZNRF3.RNA.BCR.KM$plot + theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA,
                              size = 2),
  axis.line = element_line(size = 0.5, 
                           color='black'),
  axis.text.x = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.text.y = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.title.x = element_text(size = 30,
                              face = "bold",
                              color = "black",
                              margin = margin(t=10,b=0,r=0,l=0)),
  axis.title.y = element_text(size = 36,
                              face = "bold",
                              color = "black",
                              margin = margin(t=0,b=0,r=10,l=10)),
  legend.title = element_text(size = 24,
                              face = "bold",
                              color = "black"),
  legend.text = element_text(size = 22,
                             face = "italic",
                             color = "black"),
  legend.background = element_blank(),
  plot.margin = margin(1,1,1,1, "cm")
)

ZLOTTA.ZNRF3.RNA.BCR.KM

#### KM curve of metastatic relapse in ZLOTTA patients, stratified by ZNRF3 RNA abundance (median dichotomized) - LTRI ####
fit.mets <- survfit(Surv(time_to_mets, mets) ~ znrf3.rna.status, data = ZLOTTA)
ZLOTTA.ZNRF3.RNA.METS.KM <- ggsurvplot(
  fit.mets,
  data = ZLOTTA,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time Post-Treatment (Years)",
  ylab = "Metastasis-Free Survival",
  xlim = c(0,9),
  ylim = c(0,1.01),
  break.time.by = 2,
  axes.offset = FALSE,
  font.legend =
    c(16),
  font.x =
    c(26, "bold"),
  font.y =
    c(26,"bold"),
  legend =
    c(.21, 0.35),
  legend.labs =
    c("High", "Low"),
  legend.title = "ZNRF3 RNA",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 6
)

ZLOTTA.ZNRF3.RNA.METS.KM$table <- ZLOTTA.ZNRF3.RNA.METS.KM$table +
  theme(
    plot.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,1,1,1,"cm")
  )

ZLOTTA.ZNRF3.RNA.METS.KM$plot <- ZLOTTA.ZNRF3.RNA.METS.KM$plot +
  ggplot2::annotate("text",
                    x=0.18, 
                    y=0.18,
                    label = "HR = 6.63 (1.24, 35.5)", 
                    hjust=0, 
                    size=10)
ZLOTTA.ZNRF3.RNA.METS.KM$plot <- ZLOTTA.ZNRF3.RNA.METS.KM$plot +
  ggplot2::annotate("text",
                    x= 0.18, 
                    y=0.08,
                    label = "P: 0.027",
                    hjust=0,
                    size=10,
                    parse = TRUE)

ZLOTTA.ZNRF3.RNA.METS.KM$plot <- ZLOTTA.ZNRF3.RNA.METS.KM$plot + theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA,
                              size = 2),
  axis.line = element_line(size = 0.5, 
                           color='black'),
  axis.text.x = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.text.y = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.title.x = element_text(size = 30,
                              face = "bold",
                              color = "black",
                              margin = margin(t=10,b=0,r=0,l=0)),
  axis.title.y = element_text(size = 28,
                              face = "bold",
                              color = "black",
                              margin = margin(t=0,b=0,r=10,l=10)),
  legend.title = element_text(size = 28,
                              face = "bold",
                              color = "black"),
  legend.text = element_text(size = 24,
                             face = "italic",
                             color = "black"),
  legend.background = element_blank(),
  plot.margin = margin(1,1,1,1, "cm")
)

ZLOTTA.ZNRF3.RNA.METS.KM

#
# LOAD DATA - TCGA ####
saveRDS(TCGA.ZNRF3.RNA, "/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA.ZNRF3.RNA.rds")
TCGA.ZNRF3.RNA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA.ZNRF3.RNA.rds")

# COX PROPORTIONAL HAZARDS MODELS - TCGA ####
TCGA.ZNRF3.RNA.PFS.COX <- coxph(Surv(time, PFS)~
                                  RNA.BIN, data = TCGA.ZNRF3.RNA)
summary(TCGA.ZNRF3.RNA.PFS.COX)

#KM curve of 36-month biochemical relapse in TCGA patients, stratified by ZNRF3 RNA abundance (upper quartile vs lower three quartiles) - TCGA ####

fit <- survfit(Surv(time, PFS) ~ RNA.BIN, data = TCGA.ZNRF3.RNA)
TCGA.ZNRF3.RNA.PFS.KM <- ggsurvplot(
  fit,
  data = TCGA.ZNRF3.RNA,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time Post-Treatment (Months)",
  ylab = "Progression-Free Survival",
  xlim = c(0,125),
  ylim = c(0,1.01),
  break.time.by = 24,
  axes.offset = FALSE,
  font.legend =
    c(16),
  font.x =
    c(26, "bold"),
  font.y =
    c(26,"bold"),
  legend =
    c(.23, 0.4),
  legend.labs =
    c("High", "Low"),
  legend.title = "ZNRF3 RNA",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 6
)

TCGA.ZNRF3.RNA.PFS.KM$table <- TCGA.ZNRF3.RNA.PFS.KM$table +
  theme(
    plot.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,1,1,1,"cm")
  )

TCGA.ZNRF3.RNA.PFS.KM$plot <- TCGA.ZNRF3.RNA.PFS.KM$plot +
  ggplot2::annotate("text",
                    x=0.5, 
                    y=0.2,
                    label = "HR = 1.83 (1.18, 2.82)", 
                    hjust=0, 
                    size=10)
TCGA.ZNRF3.RNA.PFS.KM$plot <- TCGA.ZNRF3.RNA.PFS.KM$plot +
  ggplot2::annotate("text",
                    x= 0.5, 
                    y=0.1,
                    label = "P: 6.63 %*% 10^-3",
                    hjust=0,
                    size=10,
                    parse = TRUE)

TCGA.ZNRF3.RNA.PFS.KM$plot <- TCGA.ZNRF3.RNA.PFS.KM$plot + theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA,
                              size = 2),
  axis.line = element_line(size = 0.5, 
                           color='black'),
  axis.text.x = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.text.y = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.title.x = element_text(size = 30,
                              face = "bold",
                              color = "black",
                              margin = margin(t=10,b=0,r=0,l=0)),
  axis.title.y = element_text(size = 28,
                              face = "bold",
                              color = "black",
                              margin = margin(t=0,b=0,r=10,l=10)),
  legend.title = element_text(size = 30,
                              face = "bold",
                              color = "black"),
  legend.text = element_text(size = 24,
                             face = "italic",
                             color = "black"),
  legend.background = element_blank(),
  plot.margin = margin(1,1,1,1, "cm")
)

TCGA.ZNRF3.RNA.PFS.KM

#
# LOAD DATA - EOPC ####

DKFZ <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/DKFZ.ZNRF3.RNA.Outcome.rds")

# COX PROPORTIONAL HAZARDS NODEL - EOPC ####

DKFZ.ZNRF3.RNA.PFS.COX <- coxph(Surv(time, bcr)~
                                  ZNRF3.RNA.BIN, data = DKFZ)
summary(DKFZ.ZNRF3.RNA.PFS.COX)

#KM curve of 36-month biochemical relapse in TCGA patients, stratified by ZNRF3 RNA abundance (upper quartile vs lower three quartiles) ####

fit <- survfit(Surv(time, bcr) ~ ZNRF3.RNA.BIN, data = DKFZ)
DKFZ.ZNRF3.RNA.BCR.KM <- ggsurvplot(
  fit,
  data = DKFZ,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time Post-Treatment (Months)",
  ylab = "BCR-Free Survival",
  xlim = c(0,75),
  ylim = c(0,1.01),
  break.time.by = 24,
  axes.offset = FALSE,
  font.legend =
    c(16),
  font.x =
    c(26, "bold"),
  font.y =
    c(26,"bold"),
  legend =
    c(.23, 0.4),
  legend.labs =
    c("High", "Low"),
  legend.title = "ZNRF3 RNA",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 6
)

DKFZ.ZNRF3.RNA.BCR.KM$table <- DKFZ.ZNRF3.RNA.BCR.KM$table +
  theme(
    plot.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0,1,1,1,"cm")
  )

DKFZ.ZNRF3.RNA.BCR.KM$plot <- DKFZ.ZNRF3.RNA.BCR.KM$plot +
  ggplot2::annotate("text",
                    x=0.5, 
                    y=0.2,
                    label = "HR = 3.37 (1.20, 9.51)", 
                    hjust=0, 
                    size=10)
DKFZ.ZNRF3.RNA.BCR.KM$plot <- DKFZ.ZNRF3.RNA.BCR.KM$plot +
  ggplot2::annotate("text",
                    x= 0.5, 
                    y=0.1,
                    label = "P: 0.022",
                    hjust=0,
                    size=10,
                    parse = TRUE)

DKFZ.ZNRF3.RNA.BCR.KM$plot <- DKFZ.ZNRF3.RNA.BCR.KM$plot + theme(
  panel.background = element_blank(),
  panel.border = element_rect(fill = NA,
                              size = 2),
  axis.line = element_line(size = 0.5, 
                           color='black'),
  axis.text.x = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.text.y = element_text(size = 24,
                             face = "bold",
                             color = "black"),
  axis.title.x = element_text(size = 30,
                              face = "bold",
                              color = "black",
                              margin = margin(t=10,b=0,r=0,l=0)),
  axis.title.y = element_text(size = 36,
                              face = "bold",
                              color = "black",
                              margin = margin(t=0,b=0,r=10,l=10)),
  legend.title = element_text(size = 30,
                              face = "bold",
                              color = "black"),
  legend.text = element_text(size = 24,
                             face = "italic",
                             color = "black"),
  legend.background = element_blank(),
  plot.margin = margin(1,1,1,1, "cm")
)

DKFZ.ZNRF3.RNA.BCR.KM
