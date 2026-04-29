# suppfig3_original.R
# Supplementary Figure 3: P7 voxelwise model components
#           (Housing effect, Litter Size effect, Litter Order effect)
#
# Source: 2502_PerinatalEnrichment_Publication(1).Rmd
# Code copied directly from the original analysis notebook.
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Model: HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex
# t-values: HousingEE (FDR 1%, t=3.23), Size_Litter (FDR 1%, t=3.68), Litter2 (FDR 1%, t=4.32)
#
# Outputs:
#   figures/suppfig3_P7_model_components.png

library(RMINC)
library(MRIcrotome)
library(dplyr)
library(here)

dataDir <- '/well/lerch/users/wwk430/P7_EE_MaternalBehaviour/'

inputFileP7 <- file.path(dataDir, 'Project_3_Data_Classification_DETBMRC.csv')
mydataTableP7 <- read.csv(inputFileP7)
mydataTableP7 <- mydataTableP7 %>% filter(mydataTableP7$Image_Quality == "Good")

mydataTableP7$Housing <- factor(mydataTableP7$Housing)
mydataTableP7$Housing <- relevel(mydataTableP7$Housing, ref = "SH")
mydataTableP7$Sex <- factor(mydataTableP7$Sex)
mydataTableP7$Sex <- relevel(mydataTableP7$Sex, ref = "M")
mydataTableP7$Litter <- factor(mydataTableP7$Litter)
mydataTableP7$Litter <- relevel(mydataTableP7$Litter, ref = "1")

lmAllP7 <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex,
                   mydataTableP7, mask = "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
print(lmAllP7)

qAllP7 <- mincFDR(lmAllP7, mask = "Pups_DWFSE_P7_2-nlin-3_mask.mnc", method = "FDR")
print(qAllP7)

fileAllP7         <- file.path(tempdir(), "lmAllP7Housing1.mnc")
fileAllP7LitterSize  <- file.path(tempdir(), "lmAllP7LiSize1.mnc")
fileAllP7LitterOrder <- file.path(tempdir(), "lmAllP7Liorder.mnc")

mincWriteVolume(lmAllP7, fileAllP7,            "tvalue-HousingEE")
mincWriteVolume(lmAllP7, fileAllP7LitterSize,  "tvalue-Size_Litter")
mincWriteVolume(lmAllP7, fileAllP7LitterOrder, "tvalue-Litter2")

anatVolP7      <- mincArray(mincGetVolume("Pups_DWFSE_P7_2-nlin-3.mnc"))
statVolP7      <- mincArray(mincGetVolume(fileAllP7))
statVolP7Size  <- mincArray(mincGetVolume(fileAllP7LitterSize))
statVolP7Order <- mincArray(mincGetVolume(fileAllP7LitterOrder))

png(here("figures", "suppfig3_P7_model_components.png"),
    units = "cm", width = 54, height = 20, res = 150)
sliceSeries(nrow=9, ncol=2, begin=55, end=195) %>%
  anatomy(anatVolP7, low=400, high=3500) %>%
  overlay(statVolP7, low=3.23, high=6, symmetric = T) %>%
  legend("t-stats Housing effect at P7 FDR 1%") %>%
  sliceSeries() %>%
  anatomy() %>%
  overlay(statVolP7Size, low=3.68, high=6, symmetric = T) %>%
  legend("t-stats litter size effect at P7 FDR 1%") %>%
  sliceSeries() %>%
  anatomy() %>%
  overlay(statVolP7Order, low=4.32, high=6, symmetric = T) %>%
  legend("t-stats Litter order effect at P7 FDR 1%") %>%
  draw()
dev.off()

message("suppfig3_original.R complete")

sessionInfo()
