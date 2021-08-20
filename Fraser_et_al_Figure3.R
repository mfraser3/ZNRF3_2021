# FRASER ET AL - FIGURE 3A


# Load Libraries ####
library(ggplot2)
library(readr)
library(survminer)
library(survival)
library(forestmodel)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplotify)

# Define directories ####
input_dir <- file.path("~/OneDrive/Work/Projects/ZNRF3")
input_dir

# Import data ####
CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCGENE.OUTCOME.ADJUSTED.PGA.rds")

## FIGURE 3A ####
## Univariable Cox Proportional Hazards Model ####
summary(coxph(Surv(bcr.time, bcr.bin) ~ ZNRF3, data = znrf3.mvcox.table))

## Kaplan-Meier Curve ####
fit.bcr <- survfit(Surv(bcr.time, bcr.bin) ~ ZNRF3, data = znrf3.mvcox.table)
CPCG.ZNRF3.BCR.KM <- ggsurvplot(
  fit.bcr,
  data = znrf3.mvcox.table,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time Post-Treatment (Months)",
  ylab = "BCR-Free Survival",
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
    c(.25, 0.3),
  legend.labs =
    c("Neutral", "Loss"),
  legend.title = "ZNRF3 CNA Status",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 7
)

CPCG.ZNRF3.BCR.KM$table <- CPCG.ZNRF3.BCR.KM$table +
  theme(
    plot.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(1,1,1,1,"cm")
  )

CPCG.ZNRF3.BCR.KM$plot <- CPCG.ZNRF3.BCR.KM$plot +
  ggplot2::annotate("text",
                    x=1, 
                    y=0.15,
                    label = "HR = 2.18 (1.31, 3.64)", 
                    hjust=0, 
                    size=10)
CPCG.ZNRF3.BCR.KM$plot <- CPCG.ZNRF3.BCR.KM$plot +
  ggplot2::annotate("text",
                    x= 1, 
                    y=0.08,
                    label = "P: 2.87 %*% 10^-3",
                    hjust=0,
                    size=10,
                    parse = TRUE)

CPCG.ZNRF3.BCR.KM$plot <- CPCG.ZNRF3.BCR.KM$plot + theme(
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
  axis.title.x = element_text(size = 36,
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

CPCG.ZNRF3.BCR.KM


## FIGURE 3B ####
## Forest Plot - BCR####
CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  drop_na(mets.time, isup2) 

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA  %>%
  mutate(ZNRF3.CNA = recode(ZNRF3.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(psa.bin = recode(psa.bin, `0` = "< 10ng/ml", `1` = "\U2265 10ng/ml")) %>%
  mutate(tcat.bin = recode(tcat.bin, `0` = "cT1", `1` = "cT2a/b", `3` = "cT2c"))

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA2 <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  transmute(bcr.bin,
            bcr.time,
            ZNRF3 = factor(ZNRF3.CNA),
            PSA = factor(psa.bin),
            TCat = factor(tcat.bin),
            PGA = adjusted.pga,
            ISUP = factor(isup2))

znrf3.mvcox.table <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  transmute(bcr.time,
            bcr.bin,
            `ZNRF3`= factor(ZNRF3.CNA),
            `ISUP Grade`= factor(isup2),
            `PSA (Continuous)` = psa,
            `Clinical T-Category` = factor(tcat.bin))

znrf3.mvcox.table$ZNRF3 <- relevel(znrf3.mvcox.table$ZNRF3, ref = "Neutral")

median(CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA$adjusted.pga)
panels <- list(
  list(width = 0.03),
  list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.1, display = ~level),
  list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.55, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(
    width = 0.05,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)


ZNRF3.CPCG.FOREST <- forest_model(coxph(Surv(bcr.time, bcr.bin) ~ ., znrf3.mvcox.table), panels)

ZNRF3.CPCG.FOREST +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black"),
    aspect.ratio = 0.5
  )

## FIGURE 3C ####
# Univariable CoxPH Model ####
summary(coxph(Surv(mets.time, mets.bin) ~ ZNRF3, data = znrf3.mvcox.table))

# Kaplan-Meier Curve ####
fit.mets <- survfit(Surv(mets.time, mets.bin) ~ ZNRF3, data = znrf3.mvcox.table)
CPCG.ZNRF3.METS.KM <- ggsurvplot(
  fit.mets,
  data = znrf3.mvcox.table,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  xlab = "Time Post-Treatment (Years)",
  ylab = "Metastasis-Free Survival",
  xlim = c(0,11),
  ylim = c(0,1.01),
  break.time.by = 2.5,
  axes.offset = FALSE,
  font.legend =
    c(16),
  font.x =
    c(26, "bold"),
  font.y =
    c(26,"bold"),
  legend =
    c(.25, 0.3),
  legend.labs =
    c("Neutral", "Loss"),
  legend.title = "ZNRF3 CNA Status",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  fontsize = 7
)

CPCG.ZNRF3.METS.KM$table <- CPCG.ZNRF3.METS.KM$table +
  theme(
    plot.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(1,1,1,1,"cm")
  )

CPCG.ZNRF3.METS.KM$plot <- CPCG.ZNRF3.METS.KM$plot +
  ggplot2::annotate("text",
                    x=0.1, 
                    y=0.15,
                    label = "HR = 4.57 (2.12, 9.84)", 
                    hjust=0, 
                    size=10)
CPCG.ZNRF3.METS.KM$plot <- CPCG.ZNRF3.METS.KM$plot +
  ggplot2::annotate("text",
                    x= 0.1, 
                    y=0.08,
                    label = "P: 1.03 %*% 10^-4",
                    hjust=0,
                    size=10,
                    parse = TRUE)

CPCG.ZNRF3.METS.KM$plot <- CPCG.ZNRF3.METS.KM$plot + theme(
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
  axis.title.x = element_text(size = 36,
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

CPCG.ZNRF3.METS.KM

## FIGURE 3D ####
## Forest Plot - Metastasis ####
CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  drop_na(mets.time, isup2) 

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA  %>%
  mutate(ZNRF3.CNA = recode(ZNRF3.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(psa.bin = recode(psa.bin, `0` = "< 10ng/ml", `1` = "\U2265 10ng/ml")) %>%
  mutate(tcat.bin = recode(tcat.bin, `0` = "cT1", `1` = "cT2a/b", `3` = "cT2c"))

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA2 <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  transmute(mets.time,
            mets.bin,
            ZNRF3 = factor(ZNRF3.CNA),
            PSA = factor(psa.bin),
            TCat = factor(tcat.bin),
            PGA = adjusted.pga,
            ISUP = factor(isup2))

znrf3.mvcox.table <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  transmute(mets.time,
            mets.bin,
            `ZNRF3`= factor(ZNRF3.CNA),
            `ISUP Grade`= factor(isup2),
            `PSA (Continuous)` = psa,
            `Clinical T-Category` = factor(tcat.bin))

znrf3.mvcox.table$ZNRF3 <- relevel(znrf3.mvcox.table$ZNRF3, ref = "Neutral")

median(CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA$adjusted.pga)
panels <- list(
  list(width = 0.03),
  list(width = 0.1, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.1, display = ~level),
  list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.55, item = "forest", hjust = 0.5, heading = "Hazard ratio", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA),
  list(
    width = 0.05,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)


ZNRF3.CPCG.FOREST <- forest_model(coxph(Surv(mets.time, mets.bin) ~ ., znrf3.mvcox.table), panels)

ZNRF3.CPCG.FOREST +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black"),
    aspect.ratio = 0.5
  )
