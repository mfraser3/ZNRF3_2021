# FRASER ET AL - FIGURE S07 #

# LOAD LIBRARIES ####
library(factorial2x2)
library(tidyverse)
library(broom)
library(reporttools)
library(survminer)
library(easyGgplot2)
library(forestmodel)

# LOAD AND TIDY DATA _FIG S07A ####
CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCGENE.OUTCOME.ADJUSTED.PGA.rds")

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  drop_na(mets.time, isup, ZNRF3.RNA, fraser.sig) 

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 2,
                                  ifelse(isup >= 3, 3, "NA")))) %>%
  drop_na(isup.bin)

CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA  %>%
  mutate(ZNRF3.CNA = recode(ZNRF3.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(psa.bin = recode(psa.bin, `0` = "< 10ng/ml", `1` = "\U2265 10ng/ml")) %>%
  mutate(tcat.bin = recode(tcat.bin, `0` = "cT1", `1` = "cT2a/b", `3` = "cT2c"))



## CONSTRUCT FOREST PLOT - FIG S07A ####
znrf3.mvcox.table <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  transmute(bcr.time,
            bcr.bin,
            `ZNRF3 RNA Abundance`= ZNRF3.RNA,
            `ISUP Grade`= factor(isup.bin),
            `PSA` = factor(psa.bin),
            `Clinical T-Category` = factor(tcat.bin))

znrf3.mvcox.table$`Clinical T-Category` <- relevel(znrf3.mvcox.table$`Clinical T-Category`, ref = "cT1")

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
    display_na = NA, hjust = 0.5, heading = "p-value (Wald test)"
  ),
  list(width = 0.03)
)


ZNRF3.CPCG.BCR.FOREST <- forest_model(coxph(Surv(bcr.time, bcr.bin) ~ `ZNRF3 RNA Abundance` + `ISUP Grade` + `PSA` + `Clinical T-Category` , znrf3.mvcox.table), panels,
                                       format_options = forest_model_format_options(text_size = 8))

ZNRF3.CPCG.BCR.FOREST +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black",
                               angle = 90,
                               hjust = 1,
                               vjust = 0.5),
    aspect.ratio = 0.9
  )



## CONSTRUCT FOREST PLOT - FIG S07B ####

znrf3.mvcox.table <- CPCG_ZNRF3_OUTCOME_ADJUSTED_PGA %>%
  transmute(bcr.time,
            bcr.bin,
            `ZNRF3 RNA Abundance`= ZNRF3.RNA,
            `ZNRF3 Loss`= factor(ZNRF3.CNA))

znrf3.mvcox.table$`ZNRF3 Loss` <- relevel(znrf3.mvcox.table$`ZNRF3 Loss`, ref = "Neutral")

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
    display_na = NA, hjust = 0.5, heading = "p-value (Wald test)"
  ),
  list(width = 0.03)
)


ZNRF3.CPCG.BCR.FOREST <- forest_model(coxph(Surv(bcr.time, bcr.bin) ~ ., znrf3.mvcox.table), panels,
                                       format_options = forest_model_format_options(text_size = 8))

ZNRF3.CPCG.BCR.FOREST +
  theme(
    axis.text.x = element_text(size = 18,
                               face = "bold",
                               color = "black",
                               angle = 90,
                               hjust = 1,
                               vjust = 0.5),
    aspect.ratio = 0.9
  )
