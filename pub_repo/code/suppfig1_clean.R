# suppfig1_clean.R
# Supplementary Figure 1: Housing effect at P43 — all subjects vs males only
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Models:
#   ~ Relative_Jacobians ~ Housing + Sex   (all P43, FDR 5%, t = 2.75)
#   ~ Relative_Jacobians ~ Housing         (P43 males only, FDR 5%, t = 3.80)
#
# Outputs:
#   figures/suppfig1_P43_housing_all_vs_males.png

library(here)
library(dplyr)
library(RMINC)
library(MRIcrotome)

DATA_DIR <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_PA  <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
ANAT_PA  <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc")
OUT_VOLS <- here("data", "derived", "stat_volumes")

# ── Data ──────────────────────────────────────────────────────────────────────

data_p43 <- read.csv(file.path(DATA_DIR, "AdultPerinatLSQ6/DataBMRC.csv")) |>
  filter(Removed_from_EE == "no", Group == "perinatal") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "shoebox"),
    Sex     = relevel(factor(Sex),     ref = "M")
  )

data_p43_males <- filter(data_p43, Sex == "M")

# ── Models ────────────────────────────────────────────────────────────────────

lm_p43_all   <- mincLm(Relative_Jacobians ~ Housing + Sex, data_p43,       mask = MASK_PA)
lm_p43_males <- mincLm(Relative_Jacobians ~ Housing,       data_p43_males, mask = MASK_PA)

cat("All P43:\n");    print(mincFDR(lm_p43_all,   mask = MASK_PA, method = "FDR"))
cat("P43 males:\n"); print(mincFDR(lm_p43_males, mask = MASK_PA, method = "FDR"))

file_all   <- file.path(OUT_VOLS, "suppfig1_lmP43_all_tstat.mnc")
file_males <- file.path(OUT_VOLS, "suppfig1_lmP43_males_tstat.mnc")

mincWriteVolume(lm_p43_all,   file_all,   "tvalue-Housingenriched")
mincWriteVolume(lm_p43_males, file_males, "tvalue-Housingenriched")

# ── Figure ────────────────────────────────────────────────────────────────────

anat_pa     <- mincArray(mincGetVolume(ANAT_PA))
stat_all    <- mincArray(mincGetVolume(file_all))
stat_males  <- mincArray(mincGetVolume(file_males))

png(here("figures", "suppfig1_P43_housing_all_vs_males.png"),
    units = "cm", width = 36, height = 20, res = 150)
sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_all,   low = 2.75, high = 6, symmetric = TRUE) |>
  legend("t-stats Housing effect at P43 (perinatally enriched) FDR 5%") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_males, low = 3.80, high = 6, symmetric = TRUE) |>
  legend("t-stats Housing effect in males only at P43 (perinatally enriched) FDR 5%") |>
  draw()
dev.off()

message("suppfig1_clean.R complete")

sessionInfo()
