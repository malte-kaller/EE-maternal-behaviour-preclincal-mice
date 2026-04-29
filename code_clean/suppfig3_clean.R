# suppfig3_clean.R
# Supplementary Figure 3: P7 voxelwise model components
#           (Housing effect, Litter Size effect, Litter Order effect)
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Model: HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex
# t-values: HousingEE (FDR 1%, t=3.23), Size_Litter (FDR 1%, t=3.68), Litter2 (FDR 1%, t=4.32)
#
# Note: if fig2_clean.R has already been run and the lmAllP7 t-stat map is saved
# in stat_volumes/, load it directly instead of re-running the model.
#
# Outputs:
#   figures/suppfig3_P7_model_components.png

library(here)
library(dplyr)
library(RMINC)
library(MRIcrotome)

DATA_DIR <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_P7  <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
ANAT_P7  <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3.mnc")
OUT_VOLS <- here("data", "derived", "stat_volumes")

# ── Data ──────────────────────────────────────────────────────────────────────

p7 <- read.csv(file.path(DATA_DIR, "Project_3_Data_Classification_DETBMRC.csv")) |>
  filter(Image_Quality == "Good") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "SH"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Litter  = relevel(factor(Litter),  ref = "1")
  )

# ── Model ─────────────────────────────────────────────────────────────────────

lm_p7 <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex,
                 p7, mask = MASK_P7)
print(mincFDR(lm_p7, mask = MASK_P7, method = "FDR"))

file_housing <- file.path(OUT_VOLS, "suppfig3_lmP7_housing_tstat.mnc")
file_size    <- file.path(OUT_VOLS, "suppfig3_lmP7_littersize_tstat.mnc")
file_order   <- file.path(OUT_VOLS, "suppfig3_lmP7_litterorder_tstat.mnc")

mincWriteVolume(lm_p7, file_housing, "tvalue-HousingEE")
mincWriteVolume(lm_p7, file_size,    "tvalue-Size_Litter")
mincWriteVolume(lm_p7, file_order,   "tvalue-Litter2")

# ── Figure ────────────────────────────────────────────────────────────────────

anat_p7    <- mincArray(mincGetVolume(ANAT_P7))
stat_house <- mincArray(mincGetVolume(file_housing))
stat_size  <- mincArray(mincGetVolume(file_size))
stat_order <- mincArray(mincGetVolume(file_order))

png(here("figures", "suppfig3_P7_model_components.png"),
    units = "cm", width = 54, height = 20, res = 150)
sliceSeries(nrow = 9, ncol = 2, begin = 55, end = 195) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_house, low = 3.23, high = 6, symmetric = TRUE) |>
  legend("t-stats Housing effect at P7 FDR 1%") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 55, end = 195) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_size,  low = 3.68, high = 6, symmetric = TRUE) |>
  legend("t-stats Litter size effect at P7 FDR 1%") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 55, end = 195) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_order, low = 4.32, high = 6, symmetric = TRUE) |>
  legend("t-stats Litter order effect at P7 FDR 1%") |>
  draw()
dev.off()

message("suppfig3_clean.R complete")

sessionInfo()
