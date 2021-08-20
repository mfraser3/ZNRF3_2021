## FRASER ET AL - FIGURE S04 ##

## LOAD LIBRARIES ####

library(readr)
library(survminer)
library(survival)
library(forestmodel)
library(tidyverse)

CPCG_CRPC_ENRICHED_OUTCOME <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCG_CRPC_ENRICHED_OUTCOME.rds")

CPCG_FOUR_GENES.OUTCOME <- CPCG_CRPC_ENRICHED_OUTCOME %>%
  select(cpcg.id, mets.time, mets.bin, ZNRF3.CNA, MYC.CNA, CCND1.CNA, PRKDC.CNA, TP53.CNA, ETV5.CNA, CDK12.CNA) %>%
  drop_na(mets.time, mets.bin)

CPCG_FOUR_GENES.OUTCOME <- CPCG_FOUR_GENES.OUTCOME %>%
  mutate(ZNRF3.CNA = recode(ZNRF3.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(MYC.CNA = recode(MYC.CNA, `0` = "Neutral", `1` = "Gain")) %>%
  mutate(CCND1.CNA = recode(CCND1.CNA, `0` = "Neutral", `1` = "Gain")) %>%
  mutate(PRKDC.CNA = recode(PRKDC.CNA, `0` = "Neutral", `1` = "Gain")) %>%
  mutate(TP53.CNA = recode(TP53.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(CDK12.CNA = recode(CDK12.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(ETV5.CNA = recode(ETV5.CNA, `0` = "Neutral", `1` = "Loss"))


forest.table <- CPCG_FOUR_GENES.OUTCOME %>%
  transmute(mets.time,
            mets.bin,
            `ZNRF3`= factor(ZNRF3.CNA),
            `MYC`= factor(MYC.CNA),
            `CCND1` = factor(CCND1.CNA),
            `PRKDC` = factor(PRKDC.CNA),
            `TP53` = factor(TP53.CNA),
            `CDK12` = factor(CDK12.CNA),
            `ETV5` = factor(ETV5.CNA))

forest.table$ZNRF3 <- relevel(forest.table$ZNRF3, ref = "Neutral")
forest.table$MYC <- relevel(forest.table$MYC, ref = "Neutral")
forest.table$CCND1 <- relevel(forest.table$CCND1, ref = "Neutral")
forest.table$PRKDC <- relevel(forest.table$PRKDC, ref = "Neutral")
forest.table$TP53 <- relevel(forest.table$TP53, ref = "Neutral")
forest.table$CDK12 <- relevel(forest.table$CDK12, ref = "Neutral")
forest.table$ETV5 <- relevel(forest.table$ETV5, ref = "Neutral")

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
    display = ~ ifelse(reference, "", format.pval(p.value, digits = 2, eps = 0.001)),
    display_na = NA, hjust = 1, heading = "p"
  ),
  list(width = 0.03)
)


ZNRF3.CPCG.FOREST <- forest_model(coxph(Surv(mets.time, mets.bin) ~ ., forest.table), panels,
                                  format_options = forest_model_format_options(text_size = 14))

ZNRF3.CPCG.FOREST +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black"),
    aspect.ratio = 0.5
  )

summary(ZNRF3.CPCG.FOREST)
ggsave(filename = "/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Communications/Supplementary Figures/Driversv12-FigureS04.tiff", width = 20, height = 10, plot = ZNRF3.CPCG.FOREST)
