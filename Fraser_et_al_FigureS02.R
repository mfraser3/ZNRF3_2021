## FRASER ET AL - FIGURE S2
#3 LOAD LIBRARIES ####
library(tidyverse)
library(readr)
library(reshape2)
library(tidyr)
library(stringr)
library(mudata2)
library(scales)
library(cowplot)
library(patchwork)

## LOAD AND TIDY DATA ####
Drivers_All_Mutations <- readRDS("/Users/michaelfraser/OneDrive/Work/Manuscripts/2020/ZNRF3/FINAL/Nature Cancer/Data and Code/Drivers_All_Mutations.rds")

Drivers_Long <- melt(Drivers_All_Mutations, id.vars = c("patient.id"), na.rm = FALSE)

Gene.Names <- str_sub(Drivers_Long$variable, 0, -6)

Mutation_type <- str_sub(Drivers_Long$variable, start = -4)

Drivers_Long <- Drivers_Long %>%
  mutate(gene.name = Gene.Names) %>%
  mutate(mutation.type = Mutation_type)

Drivers_Long$mutation.type <- recode(Drivers_Long$mutation.type, LOSS = "CNA",
                                     GAIN = "CNA",
                                     SNV1 = "SNV",
                                     SNV2 = "ncSNV",
                                     SV01 = "SV")
Drivers_Long$value <- recode(as.factor(Drivers_Long$value),
                             `2` = "1",
                             `3` = "1")

Drivers_Long[is.na(Drivers_Long)] <- "NA"

### PLOTS ####

Drivers_Heatmap <- Drivers_Long %>%
  ggplot(aes(y = patient.id, x = variable, color = as.factor(value))) +
  scale_color_manual(name = "Mutation",
                     labels = c("Absent", "Present", "NA"),
                     values = c("snow2", "dodgerblue4", "peachpuff")) +
  geom_point(shape = 15) +
  ylab("Patient") +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 45,
                                face = "bold",
                                margin = margin(t= 0, r = 10, b= 0, l = 0)),
    panel.border = element_rect(fill = NA,
                                size = 1,
                                color = "black"),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 36,
                                face = "bold"),
    legend.text = element_text(size = 24,
                               face = "bold")
  )
Drivers_Heatmap

Mutation_Colours <- c("purple", "chartreuse3", "orange", "dodgerblue4")


Drivers_Heatmap_Type_Key <- Drivers_Long %>%
  ggplot(aes(x = variable)) +
  geom_bar(aes(fill = mutation.type), alpha = 0.9) +
  scale_y_continuous(limits = c(1,1), oob = rescale_none) +
  scale_x_discrete(labels = c('NKX3-1',
                              'MYC',
                              'RB1',
                              'TP53',
                              'CDH1',
                              'PTEN',
                              'ZFHX3',
                              'AR',
                              'CHD1',
                              'ETV1',
                              'PRKDC',
                              'ETV5',
                              'BRCA2',
                              'ZNRF3',
                              'CCND1',
                              'NCOR1',
                              'TP53',
                              'ZBTB16',
                              'BRAF',
                              'APC',
                              'AKT1',
                              'CDKN2A',
                              'TTN',
                              'CTNNB1',
                              'FOXA1',
                              'ATM',
                              'SPOP',
                              'CDK12',
                              'BRCA1',
                              'MUC16',
                              'ERG',
                              'MSH6',
                              'chr21:42 Mbp',
                              'MSH2',
                              'AR',
                              'MSH2',
                              'HDAC4',
                              'CSMD3',
                              'PTEN',
                              'ASXL2',
                              'SYNE1',
                              'FOXA1',
                              'KMT2C',
                              'PTEN',
                              'AHNAK2',
                              'APC',
                              'RYR2',
                              'RB1',
                              'KMT2D',
                              'chr6:97 Mbp',
                              'BRCA2',
                              'ATM',
                              'FAT3',
                              'ZFHX3',
                              'KDM6A',
                              'RB1',
                              'SPTA1',
                              'chr7:61 Mbp',
                              'PRKDC',
                              'CTNNB1',
                              'PCDH15',
                              'MED12',
                              'chr3:125 Mbp',
                              'CHD1',
                              'chr17:25 Mbp',
                              'chr3:129 Mbp',
                              'NKX3-1',
                              'TP53',
                              'NCOR2',
                              'chr10:89 Mbp',
                              'PRB4',
                              'GNAS',
                              'ETV1',
                              'chr4:148 Mbp',
                              'PRKDC',
                              'ZNRF3',
                              'QRICH2',
                              'KDM6A',
                              'NCOR1',
                              'ATM',
                              'ADAMTS20',
                              'ASXL2',
                              'TPTE',
                              'ASXL2',
                              'BRCA1',
                              'NCOR1',
                              'ETV4',
                              'BRAF',
                              'KDM6A',
                              'MYC',
                              'HDAC4',
                              'PIK3CA',
                              'chr7:145019604',
                              'AKT1',
                              'BRAF',
                              'FOXA1',
                              'chr16:64910033',
                              'MSH6',
                              'ERCC2',
                              'MSH2',
                              'chr13:58055742',
                              'chr2:131394079',
                              'chr22:31826396',
                              'chr4:39684557',
                              'IDH1',
                              'ETV5',
                              'HRAS',
                              'MLH1',
                              'CDKN2A',
                              'AXL',
                              'ERCC2',
                              'MED12',
                              'CDKN2A')) +
  scale_fill_manual(name = "Mutation Type",
                    values = Mutation_Colours) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 10,
                               face = "bold",
                               angle = 90,
                               hjust = 1,
                               vjust = 0.5),
    axis.text.y = element_blank(),
    legend.title = element_text(size = 36,
                                face = "bold"),
    legend.text = element_text(size = 24,
                               face = "bold")
  )
Drivers_Heatmap_Type_Key


Drivers_Merged <- Drivers_Heatmap / Drivers_Heatmap_Type_Key + plot_layout(heights = c(75,1))
Drivers_Merged
