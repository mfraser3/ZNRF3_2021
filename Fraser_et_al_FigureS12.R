# FRASER ET AL - FIGURE S12 ##

## LOAD LIBRARIES ####
library(tidyverse)
library(survival)
library(survminer)

## LOAD DATA - EOPC ####

ZNRF3_Loss_GO_UpReg_DKFZ <- read_csv(paste0(input_dir,"/FINAL DATA TABLES/ZNRF3_Loss_GO_UpReg_DKFZ.csv"))
saveRDS(ZNRF3_Loss_GO_UpReg_DKFZ, "/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/ZNRF3_Loss_GO_UpReg_DKFZ.rds")

## BARPLOT - EOPC ####

ZNRF3_Loss_GO_UpReg_DKFZ %>%
  ggplot(aes(x = reorder(Gene.Set, -q.value), y = log.qval, fill = log.qval)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_gradient(low = "orange", high = "red") +
  scale_y_continuous(expand = c(0,1)) +
  labs(fill = expression("-Log"[10]*" Q-Value")) +
  ylab(expression("\n-Log"[10]*" Q-Value")) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 18,
                               face = "bold",
                               color = "black"),
    axis.text.x = element_text(size = 24,
                               face = "bold",
                               color = "black"),
    axis.title.x = element_text(size = 32,
                                face = "bold",
                                color = "black"),
    panel.border = element_rect(fill = NA,
                                size = 2,
                                color = "black"),
    legend.position = "none",
    aspect.ratio = 1
  )

#
## LOAD DATA - ABIDA ####

ZNRF3_Loss_GO_UpReg_Abida <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/ZNRF3_Loss_GO_UpReg_Abida.rds")

## BARPLOT - ABIDA ####
ZNRF3_Loss_GO_UpReg_Abida %>%
  ggplot(aes(x = reorder(Gene.Set, -q.value), y = log.qval, fill = log.qval)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_gradient(low = "orange", high = "red") +
  scale_y_continuous(expand = c(0,1)) +
  labs(fill = expression("-Log"[10]*" Q-Value")) +
  ylab(expression("\n-Log"[10]*" Q-Value")) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 18,
                               face = "bold",
                               color = "black"),
    axis.text.x = element_text(size = 24,
                               face = "bold",
                               color = "black"),
    axis.title.x = element_text(size = 32,
                                face = "bold",
                                color = "black"),
    panel.border = element_rect(fill = NA,
                                size = 2,
                                color = "black"),
    legend.position = "none",
    aspect.ratio = 1
  )

#
## LOAD OUTCOME DATA ####

CPCG_CRPC_ENRICHED_OUTCOME <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCG_CRPC_ENRICHED_OUTCOME.rds")

## LOG RANK TESTS ####
ZNRF3.CCND1.BCR.LOGRANK <- survdiff(Surv(bcr.time, bcr.bin) ~ ZNRF3.CCND1.BIN, data = CPCG_CRPC_ENRICHED_OUTCOME, rho = 0)
ZNRF3.CCND1.BCR.LOGRANK

## KAPLAN-MEIER CURVE ####
fit <- survfit(Surv(bcr.time, bcr.bin) ~ ZNRF3.CCND1.BIN, data = CPCG_CRPC_ENRICHED_OUTCOME)
CPCG.BCR.ZNRF3.CCND1.CNA.KM <- ggsurvplot(
  fit,
  data = CPCG_CRPC_ENRICHED_OUTCOME,
  size = 1,
  palette = 
    c("#FF0000", "#104E8B", "black", "orange"),
  conf.int = FALSE,
  pval = FALSE,
  risk.table = TRUE,
  risk.table.title = "Number At Risk",
  xlab = "Time (Months)",
  ylab = "BCR-Free Survival",
  xlim = c(0,32),
  ylim = c(0,1.01),
  break.time.by = 6,
  axes.offset = FALSE,
  font.legend =
    c(14),
  font.x =
    c(20, "bold"),
  font.y =
    c(20,"bold"),
  legend =
    c(.7, 0.705),
  legend.labs =
    c("Both Neutral",
      "CCND1 Gain",
      "ZNRF3 Loss",
      "ZNRF3 Loss, CCND1 Gain"),
  legend.title = "Status",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE
)

CPCG.BCR.ZNRF3.CCND1.CNA.KM$table <- CPCG.BCR.ZNRF3.CCND1.CNA.KM$table +
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

CPCG.BCR.ZNRF3.CCND1.CNA.KM$plot <- CPCG.BCR.ZNRF3.CCND1.CNA.KM$plot +
  ggplot2::annotate("text",
                    x=0.3, 
                    y=0.17,
                    label = "P: 5.0 %*% 10^-9",
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


CPCG.BCR.ZNRF3.CCND1.CNA.KM$plot <- CPCG.BCR.ZNRF3.CCND1.CNA.KM$plot + theme(
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
  legend.title = element_blank(),
  legend.text = element_text(size = 18,
                             face = "italic",
                             color = "black"),
  legend.background = element_blank(),
  plot.margin = margin(1,1,1,1, "cm")
)

CPCG.BCR.ZNRF3.CCND1.CNA.KM
