# FRASER ET AL - FIGURE S11 #

## LOAD LIBRARIES ####
library(tidyverse)
library(forestmodel)

## LOAD AND TIDY DATA ####

ZNRF3.ALL <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/CPCGENE.OUTCOME.ADJUSTED.PGA.rds")
ZNRF3.ALL <- ZNRF3.ALL %>%
  mutate(isup.bin = ifelse(isup == 1, 1,
                           ifelse(isup == 2, 2,
                                  ifelse(isup == 3, 3,
                                         ifelse(isup >= 4, 4, "NA"))))) %>%
  drop_na(isup.bin)

ZNRF3.ALL <- ZNRF3.ALL %>%
  mutate(ZNRF3.CNA = recode(ZNRF3.CNA, `0` = "Neutral", `1` = "Loss")) %>%
  mutate(tcat.bin = recode(tcat.bin, `0` = "cT1", `1` = "cT2a/b", `3` = "cT2c")) %>%
  mutate(idc.ca = recode(idc.ca, `0` = "Absent", `1` = "Present"))

## CONSTRUCT FOREST PLOT ####
znrf3.mvcox.table <- ZNRF3.ALL %>%
  transmute(mets.time,
            mets.bin,
            `ZNRF3`= factor(ZNRF3.CNA),
            `ISUP Grade`= factor(isup.bin),
            `PSA` = psa,
            `Clinical T-Category` = factor(tcat.bin),
            `Adjusted PGA` = adjusted.pga,
            `IDC-P/CA` = factor(idc.ca))

znrf3.mvcox.table$ZNRF3 <- relevel(znrf3.mvcox.table$ZNRF3, ref = "Neutral")

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


ZNRF3.CPCG.BCR.FOREST <- forest_model(coxph(Surv(mets.time, mets.bin) ~ `ZNRF3` + `ISUP Grade` + `PSA` + `Clinical T-Category` + `Adjusted PGA` + `IDC-P/CA`, znrf3.mvcox.table), panels,
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
