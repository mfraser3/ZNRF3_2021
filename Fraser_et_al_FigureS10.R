# FRASER ET AL - FIGURE S10 #

## LOAD LIBRARIES ####
library(tidyverse)
library(survival)
library(survminer)

## LOAD AND TIDY DATA - CPCG ####

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCGENE.OUTCOME.ADJUSTED.PGA.rds")

x <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA$ZNRF3.RNA
y <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA$XBP1.RNA

# SPEARMAN TEST - CPCG ####
cor.test(x, y, alternative = "two.sided", method = "spearman")

# SCATTERPLOT - CPCG ####
ZNRF3.XBP1.RNA <- ggplot(data = CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA, aes(x = ZNRF3.RNA, y = XBP1.RNA)) +
  geom_point() +
  geom_smooth(colour = "black", method = "lm", se = FALSE) +
  xlab("ZNRF3 RNA Abundance") +
  ylab("XPB1 RNA Abundance") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 2),
    axis.text = element_text(face = "bold",
                             size = 32,
                             colour = "black"),
    axis.title.x = element_text(face = "bold", size = 45,
                                margin = margin(t= 10, r = 0, b= 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = 45,
                                margin = margin(t=0,r=10,b=0,l=0))
  ) +
  annotate("text",
           label = "\U03C1 = 0.372",
           x = 7,
           y = 12.25,
           size = 14,
           hjust = 0) +
  annotate("text",
           label = "P: 3.18 %*% 10^-8",
           x = 7,
           y = 12.1,
           size = 14,
           hjust = 0,
           parse = TRUE)

ZNRF3.XBP1.RNA

## LOAD DATA - TCGA ####
TCGA_ZNRF3_XBP1_RNA_PFS <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/TCGA_ZNRF3_XBP1_RNA_PFS.rds")

# SPEARMAN TEST - TCGA ####

x <- TCGA_ZNRF3_XBP1_RNA_PFS$ZNRF3.RNA
y <- TCGA_ZNRF3_XBP1_RNA_PFS$XBP1.RNA

cor.test(x, y, alternative = "two.sided", method = "spearman")

## SCATTERPLOT - TCGA ####
ZNRF3.XBP1.RNA <- ggplot(data = TCGA_ZNRF3_XBP1_RNA_PFS, aes(x = ZNRF3.RNA, y = XBP1.RNA)) +
  geom_point() +
  geom_smooth(colour = "black", method = "lm", se = FALSE) +
  xlab("ZNRF3 RNA Abundance") +
  ylab("XPB1 RNA Abundance") +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 2),
    axis.text = element_text(face = "bold",
                             size = 32,
                             colour = "black"),
    axis.title.x = element_text(face = "bold", size = 45,
                                margin = margin(t= 10, r = 0, b= 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = 45,
                                margin = margin(t=0,r=10,b=0,l=0))
  ) +
  annotate("text",
           label = "\U03C1 = 0.171",
           x = 7,
           y = 12.25,
           size = 14,
           hjust = 0) +
  annotate("text",
           label = "P: 1.37 %*% 10^-4",
           x = 7,
           y = 12.1,
           size = 14,
           hjust = 0,
           parse = TRUE)

ZNRF3.XBP1.RNA
## COX PROPORTIONAL HAZARDS MODELS - TCGA ####

TCGA_ZNRF3_XBP1_RNA_PFS <- TCGA_ZNRF3_XBP1_RNA_PFS %>%
  mutate(ZNRF3.XBP1.RNA.BIN = as.numeric(ifelse(ZNRF3.RNA.BIN == 0 & XBP1.RNA.BIN == 0, 0,
                                                ifelse(ZNRF3.RNA.BIN == 1 & XBP1.RNA.BIN == 0,1,
                                                       ifelse(ZNRF3.RNA.BIN == 0 & XBP1.RNA.BIN == 1, 2,
                                                              ifelse(ZNRF3.RNA.BIN == 1 & XBP1.RNA.BIN == 1, 3, NA))))))

summary(coxph(Surv(pfs.time, pfs.bin) ~ ZNRF3.RNA.BIN + XBP1.RNA.BIN, data = TCGA_ZNRF3_XBP1_RNA_PFS))

## LOG RANK TEST - TCGA ####
survdiff(Surv(pfs.time, pfs.bin) ~ ZNRF3.XBP1.RNA.BIN, data = TCGA_ZNRF3_XBP1_RNA_PFS)

# KAPLAN MEIER CURVE - TCGA ####

fit <- survfit(Surv(pfs.time, pfs.bin) ~ ZNRF3.XBP1.RNA.BIN, data = TCGA_ZNRF3_XBP1_RNA_PFS)
TCGA.ZNRF3.XBP1.RNA.KM <- ggsurvplot(
  fit,
  data = TCGA_ZNRF3_XBP1_RNA_PFS,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B", "black", "chartreuse3"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  risk.table.title = "Number At Risk",
  xlab = "Time Post-Treatment (Months)",
  ylab = "Progression-Free Survival",
  xlim = c(0,122),
  ylim = c(0,1.01),
  break.time.by = 24,
  axes.offset = FALSE,
  font.legend =
    c(16),
  font.x =
    c(20, "bold"),
  font.y =
    c(20,"bold"),
  legend =
    c(.25, 0.38),
  legend.labs =
    c("Both High",
      "ZNRF3 Low",
      "XBP1 Low",
      "Both Low"
    ),
  legend.title = "ZNRF3/XBP1 RNA",
  fontsize = 6,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE
)

TCGA.ZNRF3.XBP1.RNA.KM$table <- TCGA.ZNRF3.XBP1.RNA.KM$table +
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

TCGA.ZNRF3.XBP1.RNA.KM$plot <- TCGA.ZNRF3.XBP1.RNA.KM$plot +
  ggplot2::annotate("text",
                    x=0.3, 
                    y=0.17,
                    label = "P: 3.0 %*% 10^-4",
                    hjust=0, 
                    size=8,
                    color = "black",
                    parse = TRUE) +
  ggplot2::annotate("text",
                    x=0.3, 
                    y=0.07,
                    label = "(Logrank test)",
                    hjust=0, 
                    size=8,
                    color = "black")

TCGA.ZNRF3.XBP1.RNA.KM$plot <- TCGA.ZNRF3.XBP1.RNA.KM$plot + 
  theme(
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
    axis.title.y = element_text(size = 26,
                                face = "bold",
                                color = "black",
                                margin = margin(t=0,b=0,r=10,l=10)),
    legend.title = element_text(size = 22,
                                face = "bold",
                                color = "black"),
    legend.text = element_text(size = 20,
                               face = "italic",
                               color = "black"),
    legend.background = element_blank(),
    plot.margin = margin(1,1,1,1, "cm")
  )

TCGA.ZNRF3.XBP1.RNA.KM

#
# FOREST PLOT - TCGA ####
znrf3.xbp1.mvcox.table <- TCGA_ZNRF3_XBP1_RNA_PFS %>%
  transmute(pfs.time,
            pfs.bin,
            `ZNRF3`= factor(ZNRF3.RNA.BIN),
            `XBP1` = factor(XBP1.RNA.BIN))

znrf3.xbp1.mvcox.table <- znrf3.xbp1.mvcox.table  %>%
  mutate(ZNRF3 = recode(ZNRF3, `0` = "High", `1` = "Low")) %>%
  mutate(XBP1 = recode(XBP1, `0` = "High", `1` = "Low"))


median(CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA$adjusted.pga)
panels <- list(
  list(width = 0.03),
  list(width = 0.07, display = ~variable, fontface = "bold", heading = "Variable"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.1, display = ~level, heading = "Level", fontface = "italic", hjust = 0.5),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.05, display = ~n, hjust = 0.5, heading = "N"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.75, item = "forest", hjust = 0.5, heading = "Hazard Ratio", linetype = "dashed",
    line_x = 0
  ),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(width = 0.12, display = ~ ifelse(reference, "Reference", sprintf(
    "%0.2f (%0.2f, %0.2f)",
    trans(estimate), trans(conf.low), trans(conf.high)
  )), display_na = NA, heading = "HR (95% CI)"),
  list(width = 0.03, item = "vline", hjust = 0.5),
  list(
    width = 0.05,
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
    display_na = NA, hjust = 0.5, heading = "p-value"
  ),
  list(width = 0.03)
)


TCGA.ZNRF3.XBP1.RNA.PFS.FOREST <- forest_model(coxph(Surv(pfs.time, pfs.bin) ~ ZNRF3 * XBP1, znrf3.xbp1.mvcox.table), panels,
                                               format_options = forest_model_format_options(text_size = 8))

TCGA.ZNRF3.XBP1.RNA.PFS.FOREST +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black",
                               angle = 90,
                               hjust = 1,
                               vjust = 0.5),
    aspect.ratio = 0.9
  )

