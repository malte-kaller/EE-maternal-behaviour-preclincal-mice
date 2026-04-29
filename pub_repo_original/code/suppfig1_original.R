# suppfig1_original.R
# Supplementary Figure 1: Housing effect at P43 — all subjects vs males only
#
# Source: 2502_PerinatalEnrichment_Publication(1).Rmd
# Code copied directly from the original analysis notebook.
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
# Data directory: /well/lerch/users/wwk430/P7_EE_MaternalBehaviour/
#
# Models:
#   lmPerinatHousing:      Relative_Jacobians ~ Housing + Sex   (all P43)
#   lmPerinatHousingMales: Relative_Jacobians ~ Housing         (P43 males only)
#
# Outputs:
#   figures/suppfig1_P43_housing_all_vs_males.png

library(RMINC)
library(MRIcrotome)
library(dplyr)
library(here)

dataDir <- '/well/lerch/users/wwk430/P7_EE_MaternalBehaviour/'

inputFile <- file.path(dataDir, 'AdultPerinatLSQ6/DataBMRC.csv')
mydataTable2 <- read.csv(inputFile)

mydataTable2 <- mydataTable2 %>% filter(mydataTable2$Removed_from_EE == "no")
dataPerinatal <- mydataTable2 %>% filter(mydataTable2$Group == "perinatal" & mydataTable2$Removed_from_EE == "no")

dataPerinatal$Housing <- factor(dataPerinatal$Housing)
dataPerinatal$Housing <- relevel(dataPerinatal$Housing, ref = "shoebox")
dataPerinatal$Sex <- factor(dataPerinatal$Sex)
dataPerinatal$Sex <- relevel(dataPerinatal$Sex, ref = "M")

# P43 all subjects
lmPerinatHousing <- mincLm(Relative_Jacobians ~ Housing + Sex, dataPerinatal,
                             mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
print(lmPerinatHousing)

qvPerinatHousing <- mincFDR(lmPerinatHousing,
                              mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc",
                              method = "FDR")
print(qvPerinatHousing)

# P43 males only
lmPerinatHousingMales <- mincLm(Relative_Jacobians ~ Housing,
                                  dataPerinatal %>% filter(dataPerinatal$Sex == "M"),
                                  mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
print(lmPerinatHousingMales)

qvPerinatHousingMales <- mincFDR(lmPerinatHousingMales,
                                   mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc",
                                   method = "FDR")
print(qvPerinatHousingMales)

filePerinat <- file.path(tempdir(), "lmPerinatHousing2.mnc")
mincWriteVolume(lmPerinatHousing, filePerinat, "tvalue-Housingenriched")

filePerinatMales <- file.path(tempdir(), "lmPerinatHousingMales2.mnc")
mincWriteVolume(lmPerinatHousingMales, filePerinatMales, "tvalue-Housingenriched")

anatVol            <- mincArray(mincGetVolume("AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc"))
statVolPerinat     <- mincArray(mincGetVolume(filePerinat))
statVolPerinatMales <- mincArray(mincGetVolume(filePerinatMales))

png(here("figures", "suppfig1_P43_housing_all_vs_males.png"),
    units = "cm", width = 36, height = 20, res = 150)
sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinat, low=2.75, high=6, symmetric = T) %>%
  legend("t-stats Housing effect at P43 (perinatally enriched) FDR 5%") %>%
  sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatMales, low=3.80, high=6, symmetric = T) %>%
  legend("t-stats Housing effect in males only at P43 (perinatally enriched) FDR 5%") %>%
  draw()
dev.off()

message("suppfig1_original.R complete")

sessionInfo()
