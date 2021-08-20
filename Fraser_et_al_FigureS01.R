############################################################################
#  VENN DIAGRAM OF DRIVER ABERRATION TYPES IN LOCALIZED AND METASTATIC PC  #
############################################################################

## LOAD LIBRARIES ####

library(VennDiagram)
library(BoutrosLab.plotting.general)


## PLOT TRIPLE VENN DIAGRAM #####

draw.triple.venn(
  area1 = 1834,
  area2 = 1759,
  area3 = 302,
  n12 = (1834-85),
  n13 = 298,
  n23 = 302,
  n123 = 298,
  category = c("CNA", "SNV", "SV"),
  fontfamily = "Arial",
  cat.cex = 4,
  cat.col = c("dodgerblue3", "chartreuse4", "pink"),
  cat.pos = c(0,0,0),
  cat.fontfamily = "Arial",
  cat.fontface = "bold",
  fill = c("dodgerblue3", "chartreuse4", "pink"),
  sep.dist = 0.03,
  scaled = TRUE,
  cex = 4
)
dev.off()
dev.new()
