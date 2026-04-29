# suppfig4_original.R
# Supplementary Figure 4: Combined EE effect models — 4-panel overlay
#   Panel 1: Common EE effect P43+P96 (t=2.85, FDR 5%)
#   Panel 2: Effect of age difference at perfusion (t=2.35, FDR 5%)
#   Panel 3: Perinatal vs adulthood EE interaction (t=2.60, FDR 20%)
#   Panel 4: Perinatal EE effect at P7 (t=2.48, FDR 5%)
#
# Source: 2502_PerinatalEnrichment_Publication(2).Rmd
# Code copied directly from the original analysis notebook.
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Outputs:
#   figures/suppfig4_combined_EE_models.png

library(RMINC)
library(MRIcrotome)
library(dplyr)
library(here)

dataDir <- '/well/lerch/users/wwk430/P7_EE_MaternalBehaviour/'

inputFile   <- file.path(dataDir, 'AdultPerinatLSQ6/DataBMRC.csv')
inputFileP7 <- file.path(dataDir, 'Project_3_Data_Classification_DETBMRC.csv')

mydataTable2 <- read.csv(inputFile) %>% filter(Removed_from_EE == "no")
mydataTableP7 <- read.csv(inputFileP7) %>% filter(Image_Quality == "Good")

mydataTable2$Housing <- factor(mydataTable2$Housing)
mydataTable2$Housing <- relevel(mydataTable2$Housing, ref = "shoebox")
mydataTable2$Sex <- factor(mydataTable2$Sex)
mydataTable2$Sex <- relevel(mydataTable2$Sex, ref = "M")
mydataTable2$Group <- factor(mydataTable2$Group)
mydataTable2$Group <- relevel(mydataTable2$Group, ref = "adult")

mydataTableP7$Housing <- factor(mydataTableP7$Housing)
mydataTableP7$Housing <- relevel(mydataTableP7$Housing, ref = "SH")
mydataTableP7$Sex <- factor(mydataTableP7$Sex)
mydataTableP7$Sex <- relevel(mydataTableP7$Sex, ref = "M")
mydataTableP7$Litter <- factor(mydataTableP7$Litter)
mydataTableP7$Litter <- relevel(mydataTableP7$Litter, ref = "1")

# Combined P43+P96 models
lmPerinatAdultHousing <- mincLm(Relative_Jacobians ~ Housing + Sex + Group, mydataTable2,
                                  mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
lmPerinatAdultHousingInt <- mincLm(Relative_Jacobians ~ Housing*Group + Sex, mydataTable2,
                                    mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")

# P7 model
lmAllP7 <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex,
                   mydataTableP7, mask = "Pups_DWFSE_P7_2-nlin-3_mask.mnc")

fileAdultPerinatH   <- file.path(tempdir(), "lmAdultPerinatHousing.mnc")
fileAdultPerinatAge <- file.path(tempdir(), "lmAdultPerinatAge.mnc")
fileAdultPerinatInt <- file.path(tempdir(), "lmAdultPerinatInt.mnc")
fileAllP7           <- file.path(tempdir(), "lmAllP7Housing1.mnc")

mincWriteVolume(lmPerinatAdultHousing,    fileAdultPerinatH,   "tvalue-Housingenriched")
mincWriteVolume(lmPerinatAdultHousing,    fileAdultPerinatAge, "tvalue-Groupperinatal")
mincWriteVolume(lmPerinatAdultHousingInt, fileAdultPerinatInt, "tvalue-Housingenriched:Groupperinatal")
mincWriteVolume(lmAllP7,                  fileAllP7,           "tvalue-HousingEE")

anatVol   <- mincArray(mincGetVolume("AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc"))
anatVolP7 <- mincArray(mincGetVolume("Pups_DWFSE_P7_2-nlin-3.mnc"))

statVolPerinatAdultH   <- mincArray(mincGetVolume(fileAdultPerinatH))
statVolPerinatAdultAge <- mincArray(mincGetVolume(fileAdultPerinatAge))
statVolPerinatAdultInt <- mincArray(mincGetVolume(fileAdultPerinatInt))
statVolP7              <- mincArray(mincGetVolume(fileAllP7))

png(here("figures", "suppfig4_combined_EE_models.png"),
    units = "cm", width = 72, height = 20, res = 150)
sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultH, low=2.85, high=6, symmetric = T) %>%
  legend("t-stats FDR 5%") %>%
  addtitle('Exposure to EE common effects') %>%
  sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultAge, low=2.35, high=6, symmetric = T) %>%
  legend("t-stats FDR 5%") %>%
  addtitle('Effect of age difference at perfusion') %>%
  sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultInt, low=2.60, high=6, symmetric = T) %>%
  legend("t-stats FDR 20%") %>%
  addtitle('Perinatal vs Adulthood exposure to EE') %>%
  sliceSeries(nrow=9, ncol=2, begin=55, end=195) %>%
  anatomy(anatVolP7, low=400, high=3500) %>%
  overlay(statVolP7, low=2.48, high=6, symmetric = T) %>%
  legend("t-stats FDR 5%") %>%
  addtitle('Perinatal exposure to EE: effect at P7') %>%
  draw()
dev.off()

message("suppfig4_original.R complete")

sessionInfo()
