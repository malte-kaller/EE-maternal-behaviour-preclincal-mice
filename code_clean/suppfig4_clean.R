# suppfig4_clean.R
# Supplementary Figure 4: Combined EE effect models — 4-panel overlay
#   Panel 1: Common EE effect P43+P96 (~ Housing + Sex + Group, FDR 5%, t=2.85)
#   Panel 2: Effect of age difference at perfusion (Groupperinatal, FDR 5%, t=2.35)
#   Panel 3: Perinatal vs adulthood EE interaction (Housing:Group, FDR 20%, t=2.60)
#   Panel 4: Perinatal EE effect at P7 (HousingEE, FDR 5%, t=2.48)
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Outputs:
#   figures/suppfig4_combined_EE_models.png

library(here)
library(dplyr)
library(RMINC)
library(MRIcrotome)

DATA_DIR <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_PA  <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
ANAT_PA  <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc")
MASK_P7  <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
ANAT_P7  <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3.mnc")
OUT_VOLS <- here("data", "derived", "stat_volumes")

# ── Data ──────────────────────────────────────────────────────────────────────

raw <- read.csv(file.path(DATA_DIR, "AdultPerinatLSQ6/DataBMRC.csv")) |>
  filter(Removed_from_EE == "no") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "shoebox"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Group   = relevel(factor(Group),   ref = "adult")
  )

p7 <- read.csv(file.path(DATA_DIR, "Project_3_Data_Classification_DETBMRC.csv")) |>
  filter(Image_Quality == "Good") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "SH"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Litter  = relevel(factor(Litter),  ref = "1")
  )

# ── Models ────────────────────────────────────────────────────────────────────

lm_common   <- mincLm(Relative_Jacobians ~ Housing + Sex + Group,    raw, mask = MASK_PA)
lm_interact <- mincLm(Relative_Jacobians ~ Housing * Group + Sex,    raw, mask = MASK_PA)
lm_p7       <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex,
                       p7, mask = MASK_P7)

cat("Common model FDR:\n");    print(mincFDR(lm_common,   mask = MASK_PA, method = "FDR"))
cat("Interaction model FDR:\n"); print(mincFDR(lm_interact, mask = MASK_PA, method = "FDR"))

file_housing <- file.path(OUT_VOLS, "suppfig4_lmPA_common_housing.mnc")
file_age     <- file.path(OUT_VOLS, "suppfig4_lmPA_age_tstat.mnc")
file_interact <- file.path(OUT_VOLS, "suppfig4_lmPA_interaction.mnc")
file_p7      <- file.path(OUT_VOLS, "suppfig4_lmP7_housing.mnc")

mincWriteVolume(lm_common,   file_housing,  "tvalue-Housingenriched")
mincWriteVolume(lm_common,   file_age,      "tvalue-Groupperinatal")
mincWriteVolume(lm_interact, file_interact, "tvalue-Housingenriched:Groupperinatal")
mincWriteVolume(lm_p7,       file_p7,       "tvalue-HousingEE")

# ── Figure ────────────────────────────────────────────────────────────────────

anat_pa   <- mincArray(mincGetVolume(ANAT_PA))
anat_p7   <- mincArray(mincGetVolume(ANAT_P7))
stat_h    <- mincArray(mincGetVolume(file_housing))
stat_age  <- mincArray(mincGetVolume(file_age))
stat_int  <- mincArray(mincGetVolume(file_interact))
stat_p7   <- mincArray(mincGetVolume(file_p7))

png(here("figures", "suppfig4_combined_EE_models.png"),
    units = "cm", width = 72, height = 20, res = 150)
sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_h,   low = 2.85, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 5%") |>
  addtitle("Exposure to EE: common effect (P43 + P96)") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_age, low = 2.35, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 5%") |>
  addtitle("Effect of age difference at perfusion") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_int, low = 2.60, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 20%") |>
  addtitle("Perinatal vs adulthood EE (interaction)") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 55, end = 195) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_p7,  low = 2.48, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 5%") |>
  addtitle("Perinatal exposure to EE: effect at P7") |>
  draw()
dev.off()

message("suppfig4_clean.R complete")

sessionInfo()
