# fig3_original.R
# Figure 3: Maternal behaviour (Fig 3B, 3C), voxelwise overlays (Fig 3D, 3E),
#           and striatum-contact scatter (Fig 3F)
#
# Source: Paper_PlottingMaternalCare.Rmd  (Fig 3B, 3C)
#         06_figures.R / 2502_PerinatalEnrichment_Publication(2).Rmd  (Fig 3D, 3E)
#         03b_P7_roi_local.R  (Fig 3F)
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
# Source: Paper_PlottingMaternalCare.Rmd
# Data:   MatBehav_per_Mother_Malte.csv (one row per dam × litter, n = 10)
# ─────────────────────────────────────────────────────────────────────────────

df_f <- read.csv(here("data", "derived", "MatBehav_per_Mother_Malte.csv"))

# Setting reference conditions (as in original Rmd)
df_f$Condition <- factor(df_f$Condition)
df_f$Condition <- relevel(df_f$Condition, ref = "SH")
df_f$Litter <- factor(df_f$Litter)
df_f$Litter <- relevel(df_f$Litter, ref = "Litter 1")

# Defining cross over and non cross over dams
df_f_bold  <- subset(df_f, Cage %in% c("CageB", "CageC"))
df_f_other <- subset(df_f, Cage %in% c("CageA", "CageD", "CageE"))

# Fig 3B left panel: total contact time by condition
plot_TC1 = ggplot(df_f, aes(x = Condition, y = Total_Contact, group = 1)) +
  geom_jitter(aes(shape = Condition), position = position_jitter(0.2), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_color_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(aspect.ratio = 2, legend.position='none', text = element_text(size=12)) +
  labs(x = "Condition", y = "Total Contact time")

# Fig 3B right panel: active nursing time by condition
plot_AN1 = ggplot(df_f, aes(x = Condition, y = ContactTime_Active_Nurse, group = 1)) +
  geom_jitter(aes(shape = Condition), position = position_jitter(0.2), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_color_manual(values = c("black", "black")) +
  scale_shape_manual(values = c(16, 17)) +
  theme_classic() +
  theme(aspect.ratio = 2, legend.position='none', text = element_text(size=12)) +
  labs(x = "Condition", y = "Active nursing time")

# Fig 3C: crossover dams only, all black lines
plot_TC3 = ggplot(df_f_bold, aes(x = Litter, y = Total_Contact, group = Cage)) +
  geom_line(data = df_f_bold, aes(linetype = "longdash", color = "black"), size = 0.5) +
  geom_point(aes(color = "black", shape=Condition), size=4) +
  scale_color_manual(values = c("black", "black")) +
  theme(legend.position = "none", aspect.ratio = 1, text = element_text(size=12)) +
  theme_classic() +
  labs(x = "Litter", y = "Total Contact time") +
  guides(linetype = "none")

ggsave(here("figures", "fig3B_total_contact_condition.png"),
       plot_TC1, height = 8, width = 5, units = "cm", dpi = 300)
ggsave(here("figures", "fig3B_active_nursing_condition.png"),
       plot_AN1, height = 8, width = 5, units = "cm", dpi = 300)
ggsave(here("figures", "fig3C_crossover_dams.png"),
       plot_TC3, height = 8, width = 10, units = "cm", dpi = 300)

message("Fig 3B and 3C saved (local)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 3F — LOCAL: striatum z-score vs total contact scatter
# Source: 03b_P7_roi_local.R
# Data:   MaternalCareAndVolume_1_.xlsx
# ─────────────────────────────────────────────────────────────────────────────

dataMCV <- read_excel(here("data", "derived", "MaternalCareAndVolume_1_.xlsx"))

# Clean column names (xlsx uses dots; match cleaned_data.RData convention)
names(dataMCV) <- tolower(gsub("\\.", "_", names(dataMCV)))

roi_df <- dataMCV |>
  filter(image_quality == "Good") |>
  mutate(
    brain_vol      = volumes,
    striatum       = (left_striatum + right_striatum) / 2,
    striatum_norm  = striatum / brain_vol,
    Housing        = factor(housing, levels = c("SH", "EE"))
  )

# Z-score: (roi/brain - SH_mean) / SD over all subjects
zscore_n <- function(x, housing, ref = "SH") {
  ref_mean <- mean(x[housing == ref], na.rm = TRUE)
  all_sd   <- sd(x, na.rm = TRUE)
  (x - ref_mean) / all_sd
}

roi_df <- roi_df |>
  mutate(striatum_z = zscore_n(striatum_norm, Housing))

cbPalette <- c("#EDD7C9", "#B76029")

# Pearson correlation
r_val <- cor(roi_df$total_contact, roi_df$striatum_z, use = "complete.obs")
cat("Striatum z-score vs total contact, r =", round(r_val, 2), "\n")

fig3F <- ggplot(roi_df, aes(x = total_contact, y = striatum_z)) +
  geom_point(aes(fill = Housing), alpha = 0.7, color = "black", shape = 21, size = 4,
             show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  scale_fill_manual(values = cbPalette) +
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
library(cowplot)

DATA_DIR  <- "/well/lerch/users/wwk430/P7_EE_MaternalBehaviour"
MASK_P7   <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
ANAT_P7   <- file.path(DATA_DIR, "Pups_DWFSE_P7_2-nlin-3.mnc")
out_vols  <- here("data", "derived", "stat_volumes")

# ── Load P7 data merged with maternal care ────────────────────────────────────
# Project_3_Data_with_MaternalCare_IDsorted_VoxWiseBMRC.csv merges pup-level
# imaging data with dam-level contact time (one Total_Contact per litter).

data_contact <- read.csv(here("data", "derived",
  "Project_3_Data_with_MaternalCare_IDsorted_VoxWiseBMRC.csv")) |>
  filter(Image_Quality == "Good") |>
  mutate(
    Housing    = relevel(factor(Housing), ref = "SH"),
    Litter     = relevel(factor(Litter),  ref = "Litter 1")
  )

# ── Voxelwise housing model (Fig 3D left panel) ───────────────────────────────

lmAllP7 <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter,
                   data_contact, mask = MASK_P7)

file_housing <- file.path(out_vols, "lmP7_tstat_Housing.mnc")
mincWriteVolume(lmAllP7, file_housing, "tvalue-HousingEE")

# ── Voxelwise contact model (Fig 3D right panel) ──────────────────────────────
# Outcome: Housing + Total_Contact + Size_Litter
# This gives the t-stat for maternal contact controlling for housing.

lmContact <- mincLm(HighB.Distortion.Corrected ~ Housing + Total_Contact + Size_Litter,
                     data_contact, mask = MASK_P7)

file_contact <- file.path(out_vols, "lmContact_tstat_TotalContact.mnc")
mincWriteVolume(lmContact, file_contact, "tvalue-Total_Contact")

# ── Load arrays for overlay ───────────────────────────────────────────────────

anatP7      <- mincArray(mincGetVolume(ANAT_P7))
statHousing <- mincArray(mincGetVolume(file_housing))
statContact <- mincArray(mincGetVolume(file_contact))

# Fig 3D: housing and contact t-stat overlays side by side
png(here("figures", "fig3D_housing_contact_overlay.png"),
    units = "cm", width = 36, height = 20, res = 150)
sliceSeries(nrow=7, ncol=2, begin=50, end=200) %>%
  anatomy(anatP7, low=400, high=3500) %>%
  overlay(statHousing, low=2.48, high=6, symmetric = T) %>%
  legend("Housing effect at P7 (t-stat)") %>%
  sliceSeries(nrow=7, ncol=2, begin=50, end=200) %>%
  anatomy(anatP7, low=400, high=3500) %>%
  overlay(statContact, low=2.48, high=6, symmetric = T) %>%
  legend("Maternal contact effect at P7 (t-stat)") %>%
  draw()
dev.off()

# ── Fig 3E: voxel-level t-stat scatter ───────────────────────────────────────

mask_vals   <- mincGetVolume(MASK_P7)
in_mask     <- which(mask_vals > 0.5)

tstat_house   <- as.numeric(mincGetVolume(file_housing))[in_mask]
tstat_contact_vals <- as.numeric(mincGetVolume(file_contact))[in_mask]

scatter_df <- data.frame(
  t_housing = tstat_house,
  t_contact = tstat_contact_vals
) |> filter(!is.na(t_housing), !is.na(t_contact))

r_val_vox <- cor(scatter_df$t_housing, scatter_df$t_contact)
cat("Voxel-wise correlation (EE vs contact t-stats):", round(r_val_vox, 3), "\n")

fig3E <- ggplot(scatter_df, aes(x = t_housing, y = t_contact)) +
  geom_point(alpha = 0.15, size = 0.4, color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(r_val_vox, 2)),
           hjust = 1.2, vjust = 1.5, size = 4) +
  theme_classic(base_size = 12) +
  labs(x = "t-statistic: Housing (EE vs SH)",
       y = "t-statistic: Maternal contact time")

ggsave(here("figures", "fig3E_voxel_tstat_scatter.png"),
       fig3E, width = 8, height = 8, units = "cm", dpi = 300)

message("fig3_original.R complete — fig3B, 3C, 3D, 3E, 3F saved to figures/")

sessionInfo()
