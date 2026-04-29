# fig3_clean.R
# Figure 3: Maternal behaviour (Fig 3B, 3C), voxelwise overlays (Fig 3D, 3E),
#           and striatum-contact scatter (Fig 3F)
#
# Fig 3B, 3C, 3F — LOCAL: reproduced from CSV/xlsx data.
# Fig 3D, 3E    — CLUSTER ONLY: requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
#
# Outputs:
#   figures/fig3B_total_contact_condition.png
#   figures/fig3B_active_nursing_condition.png
#   figures/fig3C_crossover_dams.png
#   figures/fig3D_housing_contact_overlay.png   [cluster]
#   figures/fig3E_voxel_tstat_scatter.png       [cluster]
#   figures/fig3F_striatum_contact_scatter.png

library(here)
library(ggplot2)
library(dplyr)
library(readxl)
library(effsize)

# ─────────────────────────────────────────────────────────────────────────────
# Fig 3B + 3C — LOCAL: maternal behaviour by condition
# Data: MatBehav_per_Mother_Malte.csv (one row per dam × litter, n = 10)
# Note: Total_Contact is stored in seconds; convert to minutes for plotting.
# ─────────────────────────────────────────────────────────────────────────────

df_f <- read.csv(here("data", "derived", "MatBehav_per_Mother_Malte.csv")) |>
  mutate(
    Condition = relevel(factor(Condition), ref = "SH"),
    Litter    = relevel(factor(Litter),    ref = "Litter 1"),
    Total_Contact_min        = Total_Contact / 60,
    ContactTime_Active_Nurse_min = ContactTime_Active_Nurse / 60
  )

df_f_crossover <- filter(df_f, Cage %in% c("CageB", "CageC"))

# Fig 3B left panel: total contact time by condition
fig3B_tc <- ggplot(df_f, aes(x = Condition, y = Total_Contact_min, group = 1)) +
  geom_jitter(aes(shape = Condition), position = position_jitter(0.2), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_shape_manual(values = c(SH = 16, EE = 17)) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 2, legend.position = "none") +
  labs(x = "Condition", y = "Total contact time (min)")

# Fig 3B right panel: active nursing time by condition
fig3B_an <- ggplot(df_f, aes(x = Condition, y = ContactTime_Active_Nurse_min, group = 1)) +
  geom_jitter(aes(shape = Condition), position = position_jitter(0.2), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_shape_manual(values = c(SH = 16, EE = 17)) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 2, legend.position = "none") +
  labs(x = "Condition", y = "Active nursing time (min)")

# Fig 3C: crossover dams only — all black lines and points
fig3C <- ggplot(df_f_crossover, aes(x = Litter, y = Total_Contact_min, group = Cage)) +
  geom_line(color = "black", linetype = "longdash", linewidth = 0.5) +
  geom_point(aes(shape = Condition), color = "black", size = 4) +
  scale_shape_manual(values = c(SH = 16, EE = 17)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(x = "Litter", y = "Total contact time (min)")

ggsave(here("figures", "fig3B_total_contact_condition.png"),
       fig3B_tc, height = 8, width = 5, units = "cm", dpi = 300)
ggsave(here("figures", "fig3B_active_nursing_condition.png"),
       fig3B_an, height = 8, width = 5, units = "cm", dpi = 300)
ggsave(here("figures", "fig3C_crossover_dams.png"),
       fig3C, height = 8, width = 10, units = "cm", dpi = 300)

message("Fig 3B and 3C saved (local)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 3F — LOCAL: striatum z-score vs total contact scatter
# Data: MaternalCareAndVolume_1_.xlsx
# Note: Total_Contact stored in seconds; x-axis label reflects this.
# ─────────────────────────────────────────────────────────────────────────────

CB_PALETTE <- c("#EDD7C9", "#B76029")

zscore_ref <- function(x, housing, ref = "SH") {
  ref_mean <- mean(x[housing == ref], na.rm = TRUE)
  all_sd   <- sd(x, na.rm = TRUE)
  (x - ref_mean) / all_sd
}

mcv <- read_excel(here("data", "derived", "MaternalCareAndVolume_1_.xlsx"))
names(mcv) <- tolower(gsub("\\.", "_", names(mcv)))

roi_df <- mcv |>
  filter(image_quality == "Good") |>
  mutate(
    Housing       = factor(housing, levels = c("SH", "EE")),
    striatum_norm = (left_striatum + right_striatum) / 2 / volumes,
    striatum_z    = zscore_ref(striatum_norm, housing)
  )

r_val <- cor(roi_df$total_contact, roi_df$striatum_z, use = "complete.obs")
cat("Striatum z-score vs total contact, r =", round(r_val, 2), "\n")

fig3F <- ggplot(roi_df, aes(x = total_contact, y = striatum_z)) +
  geom_point(aes(fill = Housing), alpha = 0.7, color = "black", shape = 21, size = 4,
             show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  scale_fill_manual(values = CB_PALETTE) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(r_val, 2)),
           hjust = 1.2, vjust = 1.5, size = 4) +
  ylim(-3, 3) +
  theme_classic(base_size = 14) +
  labs(x = "Total contact time (s)",
       y = "Striatum volume (z-score)")

ggsave(here("figures", "fig3F_striatum_contact_scatter.png"),
       fig3F, width = 8, height = 8, units = "cm", dpi = 300)

message("Fig 3F saved (local)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 3D + 3E — CLUSTER ONLY below this line
# ─────────────────────────────────────────────────────────────────────────────

library(RMINC)
library(MRIcrotome)

DATA_DIR <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_P7  <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
ANAT_P7  <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3.mnc")
OUT_VOLS <- here("data", "derived", "stat_volumes")

# ── Load P7 data merged with maternal care ────────────────────────────────────

data_contact <- read.csv(here("data", "derived",
  "Project_3_Data_with_MaternalCare_IDsorted_VoxWiseBMRC.csv")) |>
  filter(Image_Quality == "Good") |>
  mutate(
    Housing = relevel(factor(Housing), ref = "SH"),
    Litter  = relevel(factor(Litter),  ref = "Litter 1")
  )

# ── Voxelwise models ──────────────────────────────────────────────────────────

lm_housing <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter,
                     data_contact, mask = MASK_P7)

lm_contact <- mincLm(HighB.Distortion.Corrected ~ Housing + Total_Contact + Size_Litter,
                     data_contact, mask = MASK_P7)

file_housing <- file.path(OUT_VOLS, "lmP7_tstat_Housing.mnc")
file_contact <- file.path(OUT_VOLS, "lmContact_tstat_TotalContact.mnc")

mincWriteVolume(lm_housing, file_housing, "tvalue-HousingEE")
mincWriteVolume(lm_contact, file_contact, "tvalue-Total_Contact")

# ── Load arrays ───────────────────────────────────────────────────────────────

anat_p7      <- mincArray(mincGetVolume(ANAT_P7))
stat_housing <- mincArray(mincGetVolume(file_housing))
stat_contact <- mincArray(mincGetVolume(file_contact))

# Fig 3D: housing and contact t-stat overlays side by side
png(here("figures", "fig3D_housing_contact_overlay.png"),
    units = "cm", width = 36, height = 20, res = 150)
sliceSeries(nrow = 7, ncol = 2, begin = 50, end = 200) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_housing, low = 2.48, high = 6, symmetric = TRUE) |>
  legend("Housing effect at P7 (t-stat)") |>
  sliceSeries(nrow = 7, ncol = 2, begin = 50, end = 200) |>
  anatomy(anat_p7, low = 400, high = 3500) |>
  overlay(stat_contact, low = 2.48, high = 6, symmetric = TRUE) |>
  legend("Maternal contact effect at P7 (t-stat)") |>
  draw()
dev.off()

# ── Fig 3E: voxel-level t-stat scatter ───────────────────────────────────────

mask_vals <- mincGetVolume(MASK_P7)
in_mask   <- which(mask_vals > 0.5)

scatter_df <- data.frame(
  t_housing = as.numeric(mincGetVolume(file_housing))[in_mask],
  t_contact = as.numeric(mincGetVolume(file_contact))[in_mask]
) |> filter(!is.na(t_housing), !is.na(t_contact))

r_vox <- cor(scatter_df$t_housing, scatter_df$t_contact)
cat("Voxel-wise correlation (housing vs contact t-stats):", round(r_vox, 3), "\n")

fig3E <- ggplot(scatter_df, aes(x = t_housing, y = t_contact)) +
  geom_point(alpha = 0.15, size = 0.4, color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(r_vox, 2)),
           hjust = 1.2, vjust = 1.5, size = 4) +
  theme_classic(base_size = 12) +
  labs(x = "t-statistic: Housing (EE vs SH)",
       y = "t-statistic: Maternal contact time")

ggsave(here("figures", "fig3E_voxel_tstat_scatter.png"),
       fig3E, width = 8, height = 8, units = "cm", dpi = 300)

message("fig3_clean.R complete — fig3B, 3C, 3D, 3E, 3F saved to figures/")

sessionInfo()
