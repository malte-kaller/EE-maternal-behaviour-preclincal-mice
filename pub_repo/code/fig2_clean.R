# fig2_clean.R
# Figure 2: Voxelwise EE effect at P7 (Fig 2A) and brain/ROI volumes (Fig 2C)
#
# Fig 2C brain volume — LOCAL: reproduced from MaternalCareAndVolume_1_.xlsx.
# Fig 2A and Fig 2C ROI z-scores — CLUSTER ONLY (RMINC + .mnc files, Oxford BMRC).
#
# Dataset:
#   N (neonatal, P7): n = 114 total; n = 88 good-quality images used in models.
#   Atlas subset: males only for ROI z-scores (matching published analysis).
#
# Model (Fig 2A):
#   HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex
#   FDR 1% (t ≈ 3.23 for Housing)
#
# Fig 2C ROI z-score formula:
#   (roi/brain - mean_SH_males) / SD_all_males
#   Six atlas regions: Striatum dorsal region, Hindbrain, Visual areas,
#   Cerebellum, Hippocampal region, Olfactory areas.
#
# Outputs:
#   figures/fig2A_voxelwise_P7_housing.png        [cluster]
#   figures/fig2C_brain_volume.png                [local]
#   figures/fig2C_roi_striatum.png                [cluster]
#   figures/fig2C_roi_hindbrain.png               [cluster]
#   figures/fig2C_roi_visual_areas.png            [cluster]
#   figures/fig2C_roi_cerebellum.png              [cluster]
#   figures/fig2C_roi_hippocampus.png             [cluster]
#   figures/fig2C_roi_olfactory.png               [cluster]

library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(effsize)

# ─────────────────────────────────────────────────────────────────────────────
# Fig 2C — BRAIN VOLUME (local, no RMINC needed)
# ─────────────────────────────────────────────────────────────────────────────

mcv <- read_excel(here("data", "derived", "MaternalCareAndVolume_1_.xlsx")) |>
  filter(image_quality == "Good") |>
  mutate(Housing = factor(housing, levels = c("SH", "EE")))

d_brain <- cohen.d(
  mcv$volumes[mcv$Housing == "EE"],
  mcv$volumes[mcv$Housing == "SH"]
)

fig2C_brainvol <- ggplot(mcv, aes(x = Housing, y = volumes)) +
  geom_jitter(aes(shape = Housing), position = position_jitter(0.2), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_shape_manual(values = c(SH = 16, EE = 17)) +
  annotate("text", x = 1.5, y = max(mcv$volumes, na.rm = TRUE) * 1.03,
           label = paste0("d = ", round(d_brain$estimate, 2)), size = 3.5) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 2, legend.position = "none") +
  labs(x = "Housing", y = expression("Total brain volume ("*mu*"l)"))

ggsave(here("figures", "fig2C_brain_volume.png"),
       fig2C_brainvol, width = 5, height = 8, units = "cm", dpi = 300)

message("Fig 2C brain volume saved — local section complete")

# ─────────────────────────────────────────────────────────────────────────────
# CLUSTER ONLY below this point
# ─────────────────────────────────────────────────────────────────────────────

library(RMINC)
library(MRIcrotome)
library(data.tree)
library(Hmisc)
library(cowplot)

# ── Cluster paths ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_P7    <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
ANAT_P7    <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3.mnc")
P7_DEFS    <- file.path(DATA_DIR, "P7_ECox_mapping_of_labels_ABI.csv")
ALLEN_JSON <- "/well/lerch/shared/tools/atlases/Allen_Brain/Allen_hierarchy_definitions.json"

# ── 1. Load and filter P7 data ────────────────────────────────────────────────

p7 <- read.csv(file.path(DATA_DIR, "Project_3_Data_Classification_DETBMRC.csv")) |>
  filter(Image_Quality == "Good") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "SH"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Litter  = relevel(factor(Litter),  ref = "1")
  )

atlas <- read.csv(file.path(DATA_DIR, "Project_3_Voted_Atlases_ClassificationBMRC.csv")) |>
  filter(Image_Quality == "Good") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "SH"),
    Sex     = relevel(factor(Sex),     ref = "M"),
    Litter  = relevel(factor(Litter),  ref = "1")
  )

# Subsets for atlas z-scores (males only, SH reference)
atlas_males    <- filter(atlas, Sex == "M")
atlas_sh_males <- filter(atlas, Sex == "M", Housing == "SH")

# ── 2. Fig 2A: voxelwise model ────────────────────────────────────────────────

lm_p7  <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex,
                  p7, mask = MASK_P7)
fdr_p7 <- mincFDR(lm_p7, mask = MASK_P7, method = "FDR")
print(fdr_p7)  # t ≈ 3.23 at FDR 1%

file_p7_tstat <- file.path(tempdir(), "fig2A_P7_housing_tstat.mnc")
mincWriteVolume(lm_p7, file_p7_tstat, "tvalue-HousingEE")

anat_p7 <- mincArray(mincGetVolume(ANAT_P7))
stat_p7 <- mincArray(mincGetVolume(file_p7_tstat))

png(here("figures", "fig2A_voxelwise_P7_housing.png"),
    units = "cm", width = 18, height = 20, res = 150)

sliceSeries(nrow = 9, ncol = 2, begin = 55, end = 195) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_p7, low = 3.23, high = 6, symmetric = TRUE) |>
  legend("t-statistic: Housing (EE vs SH), FDR 1%") |>
  draw()

dev.off()

# ── 3. Atlas hierarchy for ROI z-scores ───────────────────────────────────────

make_hvols_p7 <- function(data, defs, allen) {
  vols <- anatGetAll(data$Atlases, method = "labels", defs = defs, side = "both")
  addVolumesToHierarchy(makeMICeDefsHierachical(defs, allen), vols)
}

hvols_males    <- make_hvols_p7(atlas_males,    P7_DEFS, ALLEN_JSON)
hvols_sh_males <- make_hvols_p7(atlas_sh_males, P7_DEFS, ALLEN_JSON)

# hanatLm on P7 atlas (all good-quality subjects, controlling for volumes)
hvolsP7_all <- make_hvols_p7(atlas, P7_DEFS, ALLEN_JSON)
hlm_p7  <- hanatLm(~ Housing + Litter + Size_Litter + Sex + Volumes,
                    data = atlas, anatTree = hvolsP7_all)
qP7 <- hanatFDR(hlm_p7)
thresholds(qP7)

# ── 4. Fig 2C ROI z-score helper ──────────────────────────────────────────────
# Formula: (roi/brain - mean_SH_males) / SD_all_males

make_roi_p7 <- function(df, hvols, hvols_sh, roi_name) {
  rn <- FindNode(hvols,    roi_name)
  sn <- FindNode(hvols_sh, roi_name)
  br <- FindNode(hvols,    "root2")
  bs <- FindNode(hvols_sh, "root2")

  df |>
    mutate(z = (rn$volumes / br$volumes - mean(sn$volumes / bs$volumes)) /
                 sd(rn$volumes / br$volumes)) |>
    ggplot(aes(x = Housing, y = z)) +
    geom_point(aes(fill = factor(Housing)),
               position = position_jitterdodge(),
               alpha = 0.7, color = "black", shape = 21, size = 4,
               show.legend = FALSE) +
    stat_summary(aes(fill = factor(Housing)),
                 fun.data = mean_cl_boot, geom = "pointrange",
                 position = position_jitterdodge(),
                 color = "black", shape = 21, size = 1,
                 show.legend = FALSE) +
    scale_fill_manual(values = c("#EDD7C9", "#B76029")) +
    ylim(-3, 3) +
    theme_classic(base_size = 14) +
    labs(x = "Housing", y = paste0(roi_name, " (z-score)"))
}

# ── 5. Fig 2C ROI panels ──────────────────────────────────────────────────────

rois <- list(
  list(name = "Striatum dorsal region", file = "fig2C_roi_striatum.png"),
  list(name = "Hindbrain",              file = "fig2C_roi_hindbrain.png"),
  list(name = "Visual areas",           file = "fig2C_roi_visual_areas.png"),
  list(name = "Cerebellum",             file = "fig2C_roi_cerebellum.png"),
  list(name = "Hippocampal region",     file = "fig2C_roi_hippocampus.png"),
  list(name = "Olfactory areas",        file = "fig2C_roi_olfactory.png")
)

for (roi in rois) {
  p <- make_roi_p7(atlas_males, hvols_males, hvols_sh_males, roi$name)
  ggsave(here("figures", roi$file), p,
         width = 7, height = 10, units = "cm", dpi = 150)
  message("Saved: ", roi$file)
}

message("fig2_clean.R complete — fig2A and fig2C panels saved to figures/")

sessionInfo()
