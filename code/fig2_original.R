# fig2_original.R
# Figure 2: Voxelwise EE effect at P7 (Fig 2A) and brain/ROI volumes (Fig 2C)
#
# Source: 2502_PerinatalEnrichment_Publication(2).Rmd
# Code copied directly from original analysis notebook.
#
# Fig 2A  — CLUSTER ONLY: requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
# Fig 2C brain volume — LOCAL: reproduced from MaternalCareAndVolume_1_.xlsx.
# Fig 2C ROI z-scores — CLUSTER ONLY: requires hanatLm / atlas .mnc files.
#
# Cluster data directory: /well/lerch/users/wwk430/P7_EE_MaternalBehaviour/
#
# Outputs:
#   figures/fig2A_voxelwise_P7_housing.png
#   figures/fig2C_brain_volume.png
#   figures/fig2C_roi_striatum.png
#   figures/fig2C_roi_hindbrain.png
#   figures/fig2C_roi_visual_areas.png
#   figures/fig2C_roi_cerebellum.png
#   figures/fig2C_roi_hippocampus.png
#   figures/fig2C_roi_olfactory.png

library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(effsize)

# ─────────────────────────────────────────────────────────────────────────────
# Fig 2C (brain volume) — LOCAL: runs from MaternalCareAndVolume_1_.xlsx
# ─────────────────────────────────────────────────────────────────────────────

dataMCV <- read_excel(here("data", "derived", "MaternalCareAndVolume_1_.xlsx"))

roi_df <- dataMCV |>
  filter(image_quality == "Good") |>
  mutate(
    brain_vol = volumes,
    Housing   = factor(housing, levels = c("SH", "EE"))
  )

d_brain <- cohen.d(
  roi_df$brain_vol[roi_df$Housing == "EE"],
  roi_df$brain_vol[roi_df$Housing == "SH"]
)

fig2C_brainvol <- ggplot(roi_df, aes(x = Housing, y = brain_vol)) +
  geom_jitter(aes(shape = Housing), position = position_jitter(0.2), size = 4) +
  stat_summary(fun = mean, geom = "crossbar", linewidth = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  scale_shape_manual(values = c(SH = 16, EE = 17)) +
  annotate("text", x = 1.5, y = max(roi_df$brain_vol, na.rm = TRUE) * 1.03,
           label = paste0("d = ", round(d_brain$estimate, 2)), size = 3.5) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 2, legend.position = "none") +
  labs(x = "Housing", y = expression("Total brain volume ("*mu*"l)"))

ggsave(here("figures", "fig2C_brain_volume.png"),
       fig2C_brainvol, width = 5, height = 8, units = "cm", dpi = 300)

message("fig2C brain volume saved (local)")

# ─────────────────────────────────────────────────────────────────────────────
# Fig 2A + Fig 2C ROI z-scores — CLUSTER ONLY below this line
# ─────────────────────────────────────────────────────────────────────────────

library(RMINC)
library(MRIcrotome)
library(data.tree)
library(Hmisc)
library(cowplot)
library(tidyverse)

dataDir <- '/well/lerch/users/wwk430/P7_EE_MaternalBehaviour/'

# ── Data loading ──────────────────────────────────────────────────────────────

inputFileP7    <- file.path(dataDir, 'Project_3_Data_Classification_DETBMRC.csv')
inputFileAtlas <- file.path(dataDir, 'Project_3_Voted_Atlases_ClassificationBMRC.csv')
defs           <- "P7_ECox_mapping_of_labels_ABI.csv"

mydataTableP7 <- read.csv(inputFileP7)
mydataTableP7 <- mydataTableP7 %>% filter(mydataTableP7$Image_Quality == "Good")

indivAtlasesP7 <- read.csv(inputFileAtlas)
indivAtlasesP7 <- indivAtlasesP7 %>% filter(indivAtlasesP7$Image_Quality == "Good")

# ── Factor releveling ─────────────────────────────────────────────────────────

mydataTableP7$Housing <- factor(mydataTableP7$Housing)
mydataTableP7$Housing <- relevel(mydataTableP7$Housing, ref = "SH")
mydataTableP7$Sex <- factor(mydataTableP7$Sex)
mydataTableP7$Sex <- relevel(mydataTableP7$Sex, ref = "M")
mydataTableP7$Litter <- factor(mydataTableP7$Litter)
mydataTableP7$Litter <- relevel(mydataTableP7$Litter, ref = "1")

indivAtlasesP7$Housing <- factor(indivAtlasesP7$Housing)
indivAtlasesP7$Housing <- relevel(indivAtlasesP7$Housing, ref = "SH")
indivAtlasesP7$Sex <- factor(indivAtlasesP7$Sex)
indivAtlasesP7$Sex <- relevel(indivAtlasesP7$Sex, ref = "M")
indivAtlasesP7$Litter <- factor(indivAtlasesP7$Litter)
indivAtlasesP7$Litter <- relevel(indivAtlasesP7$Litter, ref = "1")

# Males-only atlas subset (used for Fig 2C ROI z-scores)
indivAtlasesP7Standard <- indivAtlasesP7 %>% filter(indivAtlasesP7$Housing == "SH")
indivAtlasesP7StandardMales <- indivAtlasesP7 %>% filter(indivAtlasesP7$Housing == "SH" & indivAtlasesP7$Sex == "M")
indivAtlasesP7Males <- indivAtlasesP7 %>% filter(indivAtlasesP7$Sex == "M")

# ── Voxelwise model (Fig 2A) ──────────────────────────────────────────────────

lmAllP7 <- mincLm(HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex,
                   mydataTableP7, mask = "Pups_DWFSE_P7_2-nlin-3_mask.mnc")
print(lmAllP7)

qAllP7 <- mincFDR(lmAllP7, mask = "Pups_DWFSE_P7_2-nlin-3_mask.mnc", method = "FDR")
print(qAllP7)

fileAllP7 <- file.path(tempdir(), "lmAllP7Housing1.mnc")
mincWriteVolume(lmAllP7, fileAllP7, "tvalue-HousingEE")

# Fig 2A: P7 housing effect overlay (FDR 1%, t = 3.23)
anatVolP7 <- mincArray(mincGetVolume("Pups_DWFSE_P7_2-nlin-3.mnc"))
statVolP7  <- mincArray(mincGetVolume(fileAllP7))

png(here("figures", "fig2A_voxelwise_P7_housing.png"),
    units = "cm", width = 18, height = 20, res = 150)
sliceSeries(nrow=9, ncol=2, begin=55, end=195) %>%
  anatomy(anatVolP7, low=400, high=3500) %>%
  overlay(statVolP7, low=3.23, high=6, symmetric = T) %>%
  legend("t-stats Housing effect at P7 FDR 1%") %>%
  draw()
dev.off()

# ── Atlas volumes for Fig 2C ROI z-scores ────────────────────────────────────

allen_json  <- "/well/lerch/shared/tools/atlases/Allen_Brain/Allen_hierarchy_definitions.json"
p7_mapping  <- file.path(dataDir, "P7_ECox_mapping_of_labels_ABI.csv")

volsP7 <- anatGetAll(indivAtlasesP7$Atlases, method="labels", defs=defs, side="both")
volsP7Standard <- anatGetAll(indivAtlasesP7Standard$Atlases, method="labels", defs=defs, side="both")
volsP7StandardMales <- anatGetAll(indivAtlasesP7StandardMales$Atlases, method="labels", defs=defs, side="both")
volsP7Males <- anatGetAll(indivAtlasesP7Males$Atlases, method="labels", defs=defs, side="both")

hdefs <- makeMICeDefsHierachical(p7_mapping, allen_json)
hvolsP7 <- addVolumesToHierarchy(hdefs, volsP7)

hdefs <- makeMICeDefsHierachical(p7_mapping, allen_json)
hvolsP7Standard <- addVolumesToHierarchy(hdefs, volsP7Standard)

hdefs <- makeMICeDefsHierachical(p7_mapping, allen_json)
hvolsP7StandardMales <- addVolumesToHierarchy(hdefs, volsP7StandardMales)

hdefs <- makeMICeDefsHierachical(p7_mapping, allen_json)
hvolsP7Males <- addVolumesToHierarchy(hdefs, volsP7Males)

# hanatLm on P7 atlas (males only, controlling for volumes)
hlmP7 <- hanatLm(~ Housing + Litter + Size_Litter + Sex + Volumes,
                   data = indivAtlasesP7,
                   anatTree = hvolsP7)
qP7 <- hanatFDR(hlmP7)
thresholds(qP7)

# ── Colour palette ────────────────────────────────────────────────────────────

cbPalette <- c("#EDD7C9", "#B76029")

# ── Fig 2C ROI z-scores (males only, mean-centred on SH males, SD of all males)
# Published version: indivAtlasesP7Males / hvolsP7Males / hvolsP7StandardMales

Striatumdorsalregion <- indivAtlasesP7Males %>%
mutate(Striatumdorsalregion = (FindNode(hvolsP7Males, "Striatum dorsal region")$volumes/FindNode(hvolsP7Males, "root2")$volumes-mean(FindNode(hvolsP7StandardMales, "Striatum dorsal region")$volumes/FindNode(hvolsP7StandardMales, "root2")$volumes))/sd(FindNode(hvolsP7Males, "Striatum dorsal region")$volumes/FindNode(hvolsP7Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = Striatumdorsalregion, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=4, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPalette) +
ylim(-3,3)+
theme_classic(base_size=14)

Hindbrain <- indivAtlasesP7Males %>%
mutate(Hindbrain = (FindNode(hvolsP7Males, "Hindbrain")$volumes/FindNode(hvolsP7Males, "root2")$volumes-mean(FindNode(hvolsP7StandardMales, "Hindbrain")$volumes/FindNode(hvolsP7StandardMales, "root2")$volumes))/sd(FindNode(hvolsP7Males, "Hindbrain")$volumes/FindNode(hvolsP7Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = Hindbrain, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=4, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPalette) +
ylim(-3,3)+
theme_classic(base_size=14)

VisualAreas <- indivAtlasesP7Males %>%
mutate(VisualAreas = (FindNode(hvolsP7Males, "Visual areas")$volumes/FindNode(hvolsP7Males, "root2")$volumes-mean(FindNode(hvolsP7StandardMales, "Visual areas")$volumes/FindNode(hvolsP7StandardMales, "root2")$volumes))/sd(FindNode(hvolsP7Males, "Visual areas")$volumes/FindNode(hvolsP7Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = VisualAreas, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=4, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPalette) +
ylim(-3,3)+
theme_classic(base_size=14)

Cerebellum <- indivAtlasesP7Males %>%
mutate(Cerebellum = (FindNode(hvolsP7Males, "Cerebellum")$volumes/FindNode(hvolsP7Males, "root2")$volumes-mean(FindNode(hvolsP7StandardMales, "Cerebellum")$volumes/FindNode(hvolsP7StandardMales, "root2")$volumes))/sd(FindNode(hvolsP7Males, "Cerebellum")$volumes/FindNode(hvolsP7Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = Cerebellum, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=4, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPalette) +
ylim(-3,3)+
theme_classic(base_size=14)

Hippocampus <- indivAtlasesP7Males %>%
mutate(Hippocampus = (FindNode(hvolsP7Males, "Hippocampal region")$volumes/FindNode(hvolsP7Males, "root2")$volumes-mean(FindNode(hvolsP7StandardMales, "Hippocampal region")$volumes/FindNode(hvolsP7StandardMales, "root2")$volumes))/sd(FindNode(hvolsP7Males, "Hippocampal region")$volumes/FindNode(hvolsP7Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = Hippocampus, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=4, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPalette) +
ylim(-3,3)+
theme_classic(base_size=14)

OlfaAreasM <- indivAtlasesP7Males %>%
mutate(OlfaAreasM = (FindNode(hvolsP7Males, "Olfactory areas")$volumes/FindNode(hvolsP7Males, "root2")$volumes-mean(FindNode(hvolsP7StandardMales, "Olfactory areas")$volumes/FindNode(hvolsP7StandardMales, "root2")$volumes))/sd(FindNode(hvolsP7Males, "Olfactory areas")$volumes/FindNode(hvolsP7Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = OlfaAreasM, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=4, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPalette) +
ylim(-3,3)+
theme_classic(base_size=14)

# Save ROI panels individually (matches original Rmd output)
png(here("figures", "fig2C_roi_striatum.png"),   units="cm", width=7, height=10, res=150)
plot_grid(Striatumdorsalregion, labels = "AUTO"); dev.off()

png(here("figures", "fig2C_roi_hindbrain.png"),  units="cm", width=7, height=10, res=150)
plot_grid(Hindbrain, labels = "AUTO"); dev.off()

png(here("figures", "fig2C_roi_visual_areas.png"), units="cm", width=7, height=10, res=150)
plot_grid(VisualAreas, labels = "AUTO"); dev.off()

png(here("figures", "fig2C_roi_cerebellum.png"), units="cm", width=7, height=10, res=150)
plot_grid(Cerebellum, labels = "AUTO"); dev.off()

png(here("figures", "fig2C_roi_hippocampus.png"), units="cm", width=7, height=10, res=150)
plot_grid(Hippocampus, labels = "AUTO"); dev.off()

png(here("figures", "fig2C_roi_olfactory.png"),  units="cm", width=7, height=10, res=150)
plot_grid(OlfaAreasM, labels = "AUTO"); dev.off()

message("fig2_original.R complete — fig2A and fig2C panels saved to figures/")

sessionInfo()
