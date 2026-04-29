# suppfig2_original.R
# Supplementary Figure 2: Effect of age difference at perfusion (standalone panel)
#
# Source: 2502_PerinatalEnrichment_Publication(2).Rmd
# Code copied directly from the original analysis notebook.
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Model:
#   Relative_Jacobians ~ Housing + Sex + Group  (P43 + P96 combined)
#   t-value extracted: tvalue-Groupperinatal, FDR 5%, low = 2.35
#
# Outputs:
#   figures/suppfig2_age_difference_at_perfusion.png

library(RMINC)
library(MRIcrotome)
library(dplyr)
library(here)

dataDir <- '/well/lerch/users/wwk430/P7_EE_MaternalBehaviour/'

inputFile <- file.path(dataDir, 'AdultPerinatLSQ6/DataBMRC.csv')
mydataTable2 <- read.csv(inputFile)
mydataTable2 <- mydataTable2 %>% filter(mydataTable2$Removed_from_EE == "no")

mydataTable2$Housing <- factor(mydataTable2$Housing)
mydataTable2$Housing <- relevel(mydataTable2$Housing, ref = "shoebox")
mydataTable2$Sex <- factor(mydataTable2$Sex)
mydataTable2$Sex <- relevel(mydataTable2$Sex, ref = "M")
mydataTable2$Group <- factor(mydataTable2$Group)
mydataTable2$Group <- relevel(mydataTable2$Group, ref = "adult")

lmPerinatAdultHousing <- mincLm(Relative_Jacobians ~ Housing + Sex + Group, mydataTable2,
                                  mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
print(lmPerinatAdultHousing)

qvPerinatAdultHousing <- mincFDR(lmPerinatAdultHousing,
                                   mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc",
                                   method = "FDR")
print(qvPerinatAdultHousing)

fileAdultPerinatAge <- file.path(tempdir(), "lmAdultPerinatAge.mnc")
mincWriteVolume(lmPerinatAdultHousing, fileAdultPerinatAge, "tvalue-Groupperinatal")

anatVol <- mincArray(mincGetVolume("AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc"))
statVolPerinatAdultAge <- mincArray(mincGetVolume(fileAdultPerinatAge))

png(here("figures", "suppfig2_age_difference_at_perfusion.png"),
    units = "cm", width = 18, height = 20, res = 150)
sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultAge, low=2.35, high=6, symmetric = T) %>%
  legend("t-stats FDR 5%") %>%
  addtitle('Effect of age difference at perfusion') %>%
  draw()
dev.off()

message("suppfig2_original.R complete")

sessionInfo()
