# suppfig2_clean.R
# Supplementary Figure 2: Effect of age difference at perfusion (standalone panel)
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Model: Relative_Jacobians ~ Housing + Sex + Group  (P43 + P96 combined)
# t-value: Groupperinatal, FDR 5%, low = 2.35
#
# Note: if fig1_clean.R has already been run and saved stat_volumes/lmPA_age_tstat.mnc,
# skip the model fitting below and load that file directly.
#
# Outputs:
#   figures/suppfig2_age_difference_at_perfusion.png

library(here)
library(dplyr)
library(RMINC)
library(MRIcrotome)

DATA_DIR <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_PA  <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
ANAT_PA  <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc")
OUT_VOLS <- here("data", "derived", "stat_volumes")

# ── Data ──────────────────────────────────────────────────────────────────────

raw <- read.csv(file.path(DATA_DIR, "AdultPerinatLSQ6/DataBMRC.csv")) |>
  filter(Removed_from_EE == "no") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "shoebox"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Group   = relevel(factor(Group),   ref = "adult")
  )

# ── Model ─────────────────────────────────────────────────────────────────────

lm_common <- mincLm(Relative_Jacobians ~ Housing + Sex + Group, raw, mask = MASK_PA)
print(mincFDR(lm_common, mask = MASK_PA, method = "FDR"))

file_age <- file.path(OUT_VOLS, "suppfig2_lmPA_age_tstat.mnc")
mincWriteVolume(lm_common, file_age, "tvalue-Groupperinatal")

# ── Figure ────────────────────────────────────────────────────────────────────

anat_pa  <- mincArray(mincGetVolume(ANAT_PA))
stat_age <- mincArray(mincGetVolume(file_age))

png(here("figures", "suppfig2_age_difference_at_perfusion.png"),
    units = "cm", width = 18, height = 20, res = 150)
sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_age, low = 2.35, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 5%") |>
  addtitle("Effect of age difference at perfusion") |>
  draw()
dev.off()

message("suppfig2_clean.R complete")

sessionInfo()
