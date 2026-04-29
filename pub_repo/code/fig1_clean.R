# fig1_clean.R
# Figure 1: Voxelwise EE effects at P43/P96 (Fig 1C) and ROI z-scores (Fig 1D)
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc Jacobian/atlas files
# stored on the Oxford BMRC cluster.
#
# Datasets used:
#   P (perinatal, P43): n = 49 — perinatally enriched or standard housing
#   A (adult,    P96): n = 25 — enriched in adulthood or standard housing
#
# Models:
#   Fig 1C (left):  Relative_Jacobians ~ Housing + Sex + Group  (common EE effect)
#   Fig 1C (right): Relative_Jacobians ~ Housing * Group + Sex  (perinatal-specific)
#
# Fig 1D z-score: (roi/brain - mean_SH_P43) / SD_all_P43+P96
#   Males only for P43; all animals for P96 (matching original analysis).
#   SD is shared across P43 and P96 males to allow cross-age comparison.
#
# Outputs:
#   figures/fig1C_voxelwise_perinat_adult.png
#   figures/fig1D_roi_zscores_hippocampus.png
#   figures/fig1D_roi_zscores_visual_areas.png
#   figures/fig1D_roi_zscores_somatomotor.png
#   figures/fig1D_roi_zscores_orbital.png
#   figures/fig1D_roi_zscores_auditory.png

library(here)
library(RMINC)
library(MRIcrotome)
library(data.tree)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(cowplot)

# ── Cluster paths ─────────────────────────────────────────────────────────────
# Update if data directory changes on the cluster.

DATA_DIR    <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_PA     <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
ANAT_PA     <- file.path(DATA_DIR, "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc")
DEFS_PA     <- file.path(DATA_DIR, "AdultPerinatLSQ6/DSURQE_40micron_R_mapping.csv")
ALLEN_JSON  <- "/well/lerch/shared/tools/atlases/Allen_Brain/Allen_hierarchy_definitions.json"
DSURQE_MAP  <- "/well/lerch/shared/tools/atlases/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/mappings/DSURQE_40micron_R_mapping.csv"

# ── 1. Load and filter data ───────────────────────────────────────────────────

raw <- read.csv(file.path(DATA_DIR, "AdultPerinatLSQ6/DataBMRC.csv")) |>
  filter(Removed_from_EE == "no") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "shoebox"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Group   = relevel(factor(Group),   ref = "adult"),
    Cage    = factor(Cage)
  )

data_p43 <- filter(raw, Group == "perinatal")
data_p96 <- filter(raw, Group == "adult")

# Male-only and SH-only subsets for atlas z-scores
data_p43_males    <- filter(raw, Group == "perinatal", Sex == "M")
data_p43_sh_males <- filter(raw, Group == "perinatal", Sex == "M", Housing == "shoebox")
data_all_males    <- filter(raw, Sex == "M")
data_p96_sh       <- filter(raw, Group == "adult",     Housing == "shoebox")

# ── 2. Voxelwise models (Fig 1C) ──────────────────────────────────────────────

# Common EE effect across P43 and P96
lm_common <- mincLm(Relative_Jacobians ~ Housing + Sex + Group, raw, mask = MASK_PA)
fdr_common <- mincFDR(lm_common, mask = MASK_PA, method = "FDR")
print(fdr_common)  # t ≈ 2.85 at FDR 5%

# Perinatal-specific interaction
lm_interact <- mincLm(Relative_Jacobians ~ Housing * Group + Sex, raw, mask = MASK_PA)
fdr_interact <- mincFDR(lm_interact, mask = MASK_PA, method = "FDR")
print(fdr_interact)  # t ≈ 2.60 at FDR 20%

# Write t-stat volumes
file_common  <- file.path(tempdir(), "fig1_common_housing.mnc")
file_age     <- file.path(tempdir(), "fig1_age_group.mnc")
file_interact <- file.path(tempdir(), "fig1_perinat_interact.mnc")

mincWriteVolume(lm_common,   file_common,   "tvalue-Housingenriched")
mincWriteVolume(lm_common,   file_age,      "tvalue-Groupperinatal")
mincWriteVolume(lm_interact, file_interact, "tvalue-Housingenriched:Groupperinatal")

# Load anatomy and stat arrays
anat_pa      <- mincArray(mincGetVolume(ANAT_PA))
stat_common  <- mincArray(mincGetVolume(file_common))
stat_age     <- mincArray(mincGetVolume(file_age))
stat_interact <- mincArray(mincGetVolume(file_interact))

# Fig 1C: three-panel overlay
png(here("figures", "fig1C_voxelwise_perinat_adult.png"),
    units = "cm", width = 36, height = 20, res = 150)

sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_common,   low = 2.85, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 5%") |>
  addtitle("Exposure to EE: common effect (P43 + P96)") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_age,      low = 2.35, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 5%") |>
  addtitle("Effect of age at perfusion") |>
  sliceSeries(nrow = 9, ncol = 2, begin = 45, end = 235) |>
  anatomy(anat_pa, low = 400, high = 3500) |>
  overlay(stat_interact, low = 2.60, high = 6, symmetric = TRUE) |>
  legend("t-stats FDR 20%") |>
  addtitle("Perinatal vs adulthood EE (interaction)") |>
  draw()

dev.off()

# ── 3. Atlas volumes for Fig 1D ───────────────────────────────────────────────

make_hvols <- function(data, defs, allen) {
  vols <- anatGetAll(data$Atlases, method = "labels", defs = defs, side = "both")
  addVolumesToHierarchy(makeMICeDefsHierachical(defs, allen), vols)
}

hvols_p43_males    <- make_hvols(data_p43_males,    DSURQE_MAP, ALLEN_JSON)
hvols_p43_sh_males <- make_hvols(data_p43_sh_males, DSURQE_MAP, ALLEN_JSON)
hvols_all_males    <- make_hvols(data_all_males,    DSURQE_MAP, ALLEN_JSON)
hvols_p96          <- make_hvols(data_p96,          DSURQE_MAP, ALLEN_JSON)
hvols_p96_sh       <- make_hvols(data_p96_sh,       DSURQE_MAP, ALLEN_JSON)

# ── 4. Z-score helper ─────────────────────────────────────────────────────────
# roi_node: FindNode result for the target region in the relevant hierarchy.
# ref_node: FindNode result in the SH-only hierarchy (reference mean).
# sd_node:  FindNode result in the P43+P96 combined hierarchy (shared SD).

zscore_roi <- function(data, df, roi_node, ref_node, sd_node, var_name) {
  brain_node <- FindNode(roi_node$root, "root2")
  ref_brain  <- FindNode(ref_node$root, "root2")
  sd_brain   <- FindNode(sd_node$root, "root2")

  df |>
    mutate(
      !!var_name := (roi_node$volumes / brain_node$volumes -
                       mean(ref_node$volumes / ref_brain$volumes)) /
                      sd(sd_node$volumes  / sd_brain$volumes)
    )
}

# ── 5. Fig 1D: ROI panels ─────────────────────────────────────────────────────

cbP <- c("#ede3cb", "#bf9000")   # perinatal palette
cbA <- c("#dee9d4", "#70ad47")   # adult palette

# Shared point-range geom
roi_geom <- function(fill_pal) {
  list(
    geom_point(aes(fill = factor(Housing)),
               position = position_jitterdodge(),
               alpha = 0.7, color = "black", shape = 21, size = 3,
               show.legend = FALSE),
    stat_summary(aes(fill = factor(Housing)),
                 fun.data = mean_cl_boot, geom = "pointrange",
                 position = position_jitterdodge(),
                 color = "black", shape = 21, size = 1,
                 show.legend = FALSE),
    scale_fill_manual(values = fill_pal),
    ylim(-3, 3),
    theme_classic(base_size = 14)
  )
}

make_roi_plot <- function(df, hvols, hvols_sh, hvols_sd, roi_name, y_label, pal) {
  rn <- FindNode(hvols,    roi_name)
  sn <- FindNode(hvols_sh, roi_name)
  dn <- FindNode(hvols_sd, roi_name)
  br <- FindNode(hvols,    "root2")
  bs <- FindNode(hvols_sh, "root2")
  bd <- FindNode(hvols_sd, "root2")

  df |>
    mutate(z = (rn$volumes / br$volumes - mean(sn$volumes / bs$volumes)) /
                 sd(dn$volumes / bd$volumes)) |>
    ggplot(aes(x = Housing, y = z)) +
    roi_geom(pal) +
    labs(x = "Housing", y = paste0(y_label, " (z-score)"))
}

save_roi_fig <- function(p_perinat, p_adult, filename) {
  fig <- plot_grid(p_perinat, p_adult, labels = "AUTO")
  ggsave(here("figures", filename), fig,
         width = 14, height = 10, units = "cm", dpi = 150)
}

rois <- list(
  list(name = "Hippocampal region",  file = "fig1D_roi_zscores_hippocampus.png"),
  list(name = "Visual areas",        file = "fig1D_roi_zscores_visual_areas.png"),
  list(name = "Somatomotor areas",   file = "fig1D_roi_zscores_somatomotor.png"),
  list(name = "Orbital area",        file = "fig1D_roi_zscores_orbital.png"),
  list(name = "Auditory areas",      file = "fig1D_roi_zscores_auditory.png")
)

for (roi in rois) {
  p_p43 <- make_roi_plot(data_p43_males, hvols_p43_males, hvols_p43_sh_males,
                         hvols_all_males, roi$name, roi$name, cbP)
  p_p96 <- make_roi_plot(data_p96,       hvols_p96,       hvols_p96_sh,
                         hvols_all_males, roi$name, roi$name, cbA)
  save_roi_fig(p_p43, p_p96, roi$file)
  message("Saved: ", roi$file)
}

message("fig1_clean.R complete — fig1C and fig1D panels saved to figures/")

sessionInfo()
