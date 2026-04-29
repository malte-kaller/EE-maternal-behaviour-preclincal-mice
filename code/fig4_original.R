# fig4_original.R
# Figure 4: Voxelwise mediation of housing effect via maternal contact at P7
#           (Fig 4A, 4B), and ROI-level ACME/ADE bar charts (Fig 4C, 4D)
#
# Source: MediationBrainFigure.Rmd
#
# Fig 4C, 4D — LOCAL: hardcoded mediation effect sizes (no data files needed).
# Fig 4A, 4B — CLUSTER ONLY: requires RMINC, MRIcrotome, and .mnc files
#              at /gpfs3/well/lerch/users/flg293/MatBehavProj/EE_MatCare_Mediation_2/
#
# Outputs:
#   figures/fig4A_housing_effect_P7.png         [cluster]
#   figures/fig4B_mediation_maps.png            [cluster]
#   figures/fig4C_striatum_ACME_ADE.png
#   figures/fig4D_brainstem_ACME_ADE.png

library(here)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────
# Fig 4C — LOCAL: striatum ACME and ADE bar chart
# Source: MediationBrainFigure.Rmd
# Values: mean z-scores ± SD from voxel-wise mediation over striatum mask
# ─────────────────────────────────────────────────────────────────────────────

cbPaletteP  <- c("#ede3cb", "#bf9000")
cbPaletteA  <- c("#dee9d4", "#70ad47")
cbPaletteAP <- c("#bf9000", "#70ad47")

striatum <- data.frame(
  Effect = c("ACME", "ADE"),
  Mean   = c(-0.302900, -0.373731),
  SD     = c( 0.193988,  0.230320)
)

fig4C <- ggplot(striatum, aes(x = Effect, y = Mean, fill = Effect)) +
  geom_bar(stat = "identity",
           width = 0.45,
           colour = "black",
           linewidth = 0.6) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.12,
                linewidth = 0.8) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.6) +
  scale_fill_manual(values = c("ADE" = cbPaletteA[2],
                               "ACME" = cbPaletteP[2])) +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 14),
    axis.text     = element_text(size = 13),
    legend.position = "none"
  ) +
  xlab("Effect") +
  ylab("Z-score") +
  ggtitle("Striatum")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 4D — LOCAL: brainstem ACME and ADE bar chart
# Source: MediationBrainFigure.Rmd
# ─────────────────────────────────────────────────────────────────────────────

brainstem <- data.frame(
  Effect = c("ACME", "ADE"),
  Mean   = c(-0.087640,  0.982687),
  SD     = c( 0.163396,  0.221164)
)

fig4D <- ggplot(brainstem, aes(x = Effect, y = Mean, fill = Effect)) +
  geom_bar(stat = "identity",
           width = 0.45,
           colour = "black",
           linewidth = 0.6) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.12,
                linewidth = 0.8) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.6) +
  scale_fill_manual(values = c("ADE" = cbPaletteA[2],
                               "ACME" = cbPaletteP[2])) +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    axis.title    = element_text(size = 14),
    axis.text     = element_text(size = 13),
    legend.position = "none"
  ) +
  xlab("Effect") +
  ylab("Z-score") +
  ggtitle("Brainstem")

ggsave(here("figures", "fig4C_striatum_ACME_ADE.png"),
       fig4C, width = 8, height = 10, units = "cm", dpi = 300)
ggsave(here("figures", "fig4D_brainstem_ACME_ADE.png"),
       fig4D, width = 8, height = 10, units = "cm", dpi = 300)

message("Fig 4C and 4D saved (local)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 4A + 4B — CLUSTER ONLY below this line
# ─────────────────────────────────────────────────────────────────────────────

library(RMINC)
library(MRIcrotome)

MED_DIR  <- "/gpfs3/well/lerch/users/flg293/MatBehavProj/EE_MatCare_Mediation_2"
THR_DIR  <- file.path(MED_DIR, "data", "thr_data")

ANAT_P7  <- file.path(MED_DIR, "Pups_DWFSE_P7_2-nlin-3.mnc")

# ── Fig 4A: housing effect at P7 ─────────────────────────────────────────────

anatVol      <- mincArray(mincGetVolume(ANAT_P7))
statHousing  <- mincArray(mincGetVolume(file.path(MED_DIR, "notebook_lmAllH.mnc")))

png(here("figures", "fig4A_housing_effect_P7.png"),
    units = "cm", width = 18, height = 20, res = 150)
sliceSeries(nrow = 9, ncol = 2, begin = 50, end = 200) %>%
  anatomy(anatVol, low = 400, high = 3500) %>%
  overlay(statHousing, low = 2.48, high = 6, symmetric = T) %>%
  legend("t-stats Housing effect at P7 (FDR 5%)") %>%
  draw()
dev.off()

# ── Fig 4B: four-panel mediation z-score overlays ────────────────────────────

statTotalH  <- mincArray(mincGetVolume(file.path(THR_DIR, "z_TH_masked.mnc")))
statPropMed <- mincArray(mincGetVolume(file.path(THR_DIR, "prop_med_map_full_Housingmasked.mnc")))
statADE     <- mincArray(mincGetVolume(file.path(THR_DIR, "z_ade_masked.mnc")))
statACME    <- mincArray(mincGetVolume(file.path(THR_DIR, "z_acme_masked.mnc")))

png(here("figures", "fig4B_mediation_maps.png"),
    units = "cm", width = 24, height = 20, res = 150)
sliceSeries(nrow = 6, ncol = 1, begin = 60, end = 200) %>%
  anatomy(anatVol, low = 400, high = 3500) %>%
  overlay(statTotalH, low = 0, high = 1.5, symmetric = T) %>%
  legend("Total effect of housing (z-score)") %>%
  sliceSeries(nrow = 6, ncol = 1, begin = 60, end = 200) %>%
  anatomy(anatVol, low = 400, high = 3500) %>%
  overlay(statPropMed, low = 0, high = 1, symmetric = T) %>%
  legend("Proportional Mediation (ratio)") %>%
  sliceSeries(nrow = 6, ncol = 1, begin = 50, end = 200) %>%
  anatomy(anatVol, low = 400, high = 3500) %>%
  overlay(statADE, low = 0, high = 2, symmetric = T) %>%
  legend("Direct effect (ADE)") %>%
  sliceSeries(nrow = 6, ncol = 1, begin = 50, end = 200) %>%
  anatomy(anatVol, low = 400, high = 3500) %>%
  overlay(statACME, low = 0, high = 2, symmetric = T) %>%
  legend("Indirect effect (ACME)") %>%
  draw()
dev.off()

message("fig4_original.R complete — fig4A, 4B, 4C, 4D saved to figures/")

sessionInfo()
