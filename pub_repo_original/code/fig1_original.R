# fig1_original.R
# Figure 1: Voxelwise EE effects at P43/P96 (Fig 1C) and ROI z-scores (Fig 1D)
#
# Source: 2502_PerinatalEnrichment_Publication(2).Rmd
# Code copied directly from original analysis notebook.
#
# CLUSTER ONLY — requires RMINC, MRIcrotome, and .mnc files on Oxford BMRC.
# Cluster data directory: /well/lerch/users/wwk430/P7_EE_MaternalBehaviour/
#
# Outputs:
#   figures/fig1C_voxelwise_perinat_adult.png
#   figures/fig1D_roi_zscores_hippocampus.png
#   figures/fig1D_roi_zscores_visual_areas.png
#   figures/fig1D_roi_zscores_somatomotor.png
#   figures/fig1D_roi_zscores_orbital.png
#   figures/fig1D_roi_zscores_auditory.png

library(RMINC)
library(data.tree)
library(rjson)
library(MRIcrotome)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(cowplot)
library(here)

# ── Data loading ──────────────────────────────────────────────────────────────

dataDir <- '/well/lerch/users/wwk430/P7_EE_MaternalBehaviour/'

inputFile <- file.path(dataDir, 'AdultPerinatLSQ6/DataBMRC.csv')
mydataTable2 <- read.csv(inputFile)
print(mydataTable2)

mydataTable2 <- mydataTable2 %>% filter(mydataTable2$Removed_from_EE == "no")
dataAdult <- mydataTable2 %>% filter(mydataTable2$Group == "adult" & mydataTable2$Removed_from_EE == "no")
dataPerinatal <- mydataTable2 %>% filter(mydataTable2$Group == "perinatal" & mydataTable2$Removed_from_EE == "no")

# ── Factor releveling ─────────────────────────────────────────────────────────

mydataTable2$Housing <- factor(mydataTable2$Housing)
mydataTable2$Housing <- relevel(mydataTable2$Housing, ref = "shoebox")
mydataTable2$Sex <- factor(mydataTable2$Sex)
mydataTable2$Sex <- relevel(mydataTable2$Sex, ref = "M")
mydataTable2$Group <- factor(mydataTable2$Group)
mydataTable2$Group <- relevel(mydataTable2$Group, ref = "adult")
mydataTable2$Cage <- factor(mydataTable2$Cage)

dataAdult$Housing <- factor(dataAdult$Housing)
dataAdult$Housing <- relevel(dataAdult$Housing, ref = "shoebox")
dataAdult$Sex <- factor(dataAdult$Sex)
dataAdult$Sex <- relevel(dataAdult$Sex, ref = "M")
dataAdult$Cage <- factor(dataAdult$Cage)

dataPerinatal$Housing <- factor(dataPerinatal$Housing)
dataPerinatal$Housing <- relevel(dataPerinatal$Housing, ref = "shoebox")
dataPerinatal$Sex <- factor(dataPerinatal$Sex)
dataPerinatal$Sex <- relevel(dataPerinatal$Sex, ref = "M")
dataPerinatal$Cage <- factor(dataPerinatal$Cage)

# ── Males-only subsets ────────────────────────────────────────────────────────

mydataTable2PeriMales <- mydataTable2 %>% filter(mydataTable2$Group == "perinatal" & mydataTable2$Sex == "M")
mydataTable2SHPeriMales <- mydataTable2 %>% filter(mydataTable2$Housing == "shoebox" & mydataTable2$Group == "perinatal" & mydataTable2$Sex == "M")
mydataTable2SHMales <- mydataTable2 %>% filter(mydataTable2$Housing == "shoebox" & mydataTable2$Sex == "M")
mydataTable2Males <- mydataTable2 %>% filter(mydataTable2$Sex == "M")
mydataTable2SHAdult <- mydataTable2 %>% filter(mydataTable2$Housing == "shoebox" & mydataTable2$Group == "adult")

# ── Voxelwise models (Fig 1C) ─────────────────────────────────────────────────

# Common EE effect across P43 and P96 (Housing + Sex + Group)
lmPerinatAdultHousing <- mincLm(Relative_Jacobians ~ Housing + Sex + Group, mydataTable2,
                                 mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
print(lmPerinatAdultHousing)

qvPerinatAdultHousing <- mincFDR(lmPerinatAdultHousing,
                                  mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc",
                                  method = "FDR")
print(qvPerinatAdultHousing)

# Perinatal-specific effect (Housing * Group interaction)
lmPerinatAdultHousingInt <- mincLm(Relative_Jacobians ~ Housing*Group + Sex, mydataTable2,
                                    mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc")
print(lmPerinatAdultHousingInt)

qvPerinatAdultHousingInt <- mincFDR(lmPerinatAdultHousingInt,
                                     mask = "AdultPerinatLSQ6/MagetAdultPerinat-nlin-3_mask.mnc",
                                     method = "FDR")
print(qvPerinatAdultHousingInt)

# Write t-stat volumes to temporary files
fileAdultPerinatH <- file.path(tempdir(), "lmAdultPerinatHousing.mnc")
mincWriteVolume(lmPerinatAdultHousing, fileAdultPerinatH, "tvalue-Housingenriched")
fileAdultPerinatAge <- file.path(tempdir(), "lmAdultPerinatAge.mnc")
mincWriteVolume(lmPerinatAdultHousing, fileAdultPerinatAge, "tvalue-Groupperinatal")
fileAdultPerinatInt <- file.path(tempdir(), "lmAdultPerinatInt.mnc")
mincWriteVolume(lmPerinatAdultHousingInt, fileAdultPerinatInt, "tvalue-Housingenriched:Groupperinatal")

# Load anatomy and stat volumes
anatVol <- mincArray(mincGetVolume("AdultPerinatLSQ6/MagetAdultPerinat-nlin-3.mnc"))
statVolPerinatAdultH <- mincArray(mincGetVolume(fileAdultPerinatH))
statVolPerinatAdultAge <- mincArray(mincGetVolume(fileAdultPerinatAge))
statVolPerinatAdultInt <- mincArray(mincGetVolume(fileAdultPerinatInt))

# Fig 1C: Common housing effect and perinatal-specific interaction overlay
png(here("figures", "fig1C_voxelwise_perinat_adult.png"),
    units = "cm", width = 36, height = 20, res = 150)
sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultH, low=2.85, high=6, symmetric = T) %>%
  legend("t-stats FDR 5%") %>%
  addtitle('Exposure to EE common effects') %>%
  sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultAge, low=2.35, high=6, symmetric = T) %>%
  legend("t-stats FDR 5%") %>%
  addtitle('Effect of age difference at perfusion') %>%
  sliceSeries(nrow=9, ncol=2, begin=45, end=235) %>%
  anatomy(anatVol, low=400, high=3500) %>%
  overlay(statVolPerinatAdultInt, low=2.60, high=6, symmetric = T) %>%
  legend("t-stats FDR 20%") %>%
  addtitle('Perinatal vs Adulthood exposure to EE') %>%
  draw()
dev.off()

# ── Atlas data loading (Fig 1D) ───────────────────────────────────────────────

defsAdultPerinat <- "AdultPerinatLSQ6/DSURQE_40micron_R_mapping.csv"
allen_json <- "/well/lerch/shared/tools/atlases/Allen_Brain/Allen_hierarchy_definitions.json"
dsurqe_mapping <- "/well/lerch/shared/tools/atlases/Dorr_2008_Steadman_2013_Ullmann_2013_Richards_2011_Qiu_2016_Egan_2015_40micron/mappings/DSURQE_40micron_R_mapping.csv"

# All P43+P96 males (used for shared SD in z-score denominator)
volsP43P96Males <- anatGetAll(mydataTable2Males$Atlases, method="labels",
                               defs=defsAdultPerinat, side="both")

# SH-only P43 males (used for z-score reference mean)
volsP43StandardMales <- anatGetAll(mydataTable2SHPeriMales$Atlases, method="labels",
                                    defs=defsAdultPerinat, side="both")

# All P43 males
volsP43Males <- anatGetAll(mydataTable2PeriMales$Atlases, method="labels",
                            defs=defsAdultPerinat, side="both")

# SH-only P96 (used for z-score reference mean for adult dataset)
volsP96Standard <- anatGetAll(mydataTable2SHAdult$Atlases, method="labels",
                               defs=defsAdultPerinat, side="both")

# All P96
volsP96 <- anatGetAll(dataAdult$Atlases, method="labels",
                       defs=defsAdultPerinat, side="both")

# Build hierarchy trees
hvolsP43Males <- addVolumesToHierarchy(
  makeMICeDefsHierachical(dsurqe_mapping, allen_json), volsP43Males)

hvolsP43StandardMales <- addVolumesToHierarchy(
  makeMICeDefsHierachical(dsurqe_mapping, allen_json), volsP43StandardMales)

hvolsP43P96Males <- addVolumesToHierarchy(
  makeMICeDefsHierachical(dsurqe_mapping, allen_json), volsP43P96Males)

hvolsP96Standard <- addVolumesToHierarchy(
  makeMICeDefsHierachical(dsurqe_mapping, allen_json), volsP96Standard)

hvolsAdult <- addVolumesToHierarchy(
  makeMICeDefsHierachical(dsurqe_mapping, allen_json), volsP96)

# ── Colour palettes ───────────────────────────────────────────────────────────

cbPaletteP <- c("#ede3cb", "#bf9000")   # perinatal: SH (light), EE (dark gold)
cbPaletteA <- c("#dee9d4", "#70ad47")   # adult: SH (light), EE (dark green)

# ── Fig 1D: ROI z-score plots ─────────────────────────────────────────────────
# Z-score formula: (roi/brain_vol - mean_SH_P43) / sd_all_P43P96
# Males only for P43; all animals for P96.

# Hippocampal region
HippocampalRegionPerinat <- mydataTable2PeriMales %>%
mutate(HippocampalRegion = (FindNode(hvolsP43Males, "Hippocampal region")$volumes/FindNode(hvolsP43Males, "root2")$volumes-mean(FindNode(hvolsP43StandardMales, "Hippocampal region")$volumes/FindNode(hvolsP43StandardMales, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Hippocampal region")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = HippocampalRegion, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteP) +
ylim(-3,3)+
theme_classic(base_size=14)

HippocampalRegionAdult <- dataAdult %>%
mutate(HippocampalRegionAd = (FindNode(hvolsAdult, "Hippocampal region")$volumes/FindNode(hvolsAdult, "root2")$volumes-mean(FindNode(hvolsP96Standard, "Hippocampal region")$volumes/FindNode(hvolsP96Standard, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Hippocampal region")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = HippocampalRegionAd, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteA) +
ylim(-3,3)+
theme_classic(base_size=14)

png(here("figures", "fig1D_roi_zscores_hippocampus.png"),
    units = "cm", width = 14, height = 10, res = 150)
plot_grid(HippocampalRegionPerinat, HippocampalRegionAdult, labels = "AUTO")
dev.off()

# Lateral septal complex
LateralSeptalPerinat <- mydataTable2PeriMales %>%
mutate(LateralSeptal = (FindNode(hvolsP43Males, "Lateral septal complex")$volumes/FindNode(hvolsP43Males, "root2")$volumes-mean(FindNode(hvolsP43StandardMales, "Lateral septal complex")$volumes/FindNode(hvolsP43StandardMales, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Lateral septal complex")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = LateralSeptal, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteP) +
ylim(-3,3)+
theme_classic(base_size=14)

LateralSeptalAdult <- dataAdult %>%
mutate(LateralSeptalAd = (FindNode(hvolsAdult, "Lateral septal complex")$volumes/FindNode(hvolsAdult, "root2")$volumes-mean(FindNode(hvolsP96Standard, "Lateral septal complex")$volumes/FindNode(hvolsP96Standard, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Lateral septal complex")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = LateralSeptalAd, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteA) +
ylim(-3,3)+
theme_classic(base_size=14)

# Visual areas
VisualAreasPerinat <- mydataTable2PeriMales %>%
mutate(VisualAreasPerinat = (FindNode(hvolsP43Males, "Visual areas")$volumes/FindNode(hvolsP43Males, "root2")$volumes-mean(FindNode(hvolsP43StandardMales, "Visual areas")$volumes/FindNode(hvolsP43StandardMales, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Visual areas")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = VisualAreasPerinat, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteP) +
ylim(-3,3)+
theme_classic(base_size=14)

VisualAreasAdult <- dataAdult %>%
mutate(VisualAreasAdult = (FindNode(hvolsAdult, "Visual areas")$volumes/FindNode(hvolsAdult, "root2")$volumes-mean(FindNode(hvolsP96Standard, "Visual areas")$volumes/FindNode(hvolsP96Standard, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Visual areas")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = VisualAreasAdult, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteA) +
ylim(-3,3) +
theme_classic(base_size=14)

png(here("figures", "fig1D_roi_zscores_visual_areas.png"),
    units = "cm", width = 14, height = 10, res = 150)
plot_grid(VisualAreasPerinat, VisualAreasAdult, labels = "AUTO")
dev.off()

# Somatomotor areas
SomatoMotorPerinat <- mydataTable2PeriMales %>%
mutate(SomatoMotorPerinat = (FindNode(hvolsP43Males, "Somatomotor areas")$volumes/FindNode(hvolsP43Males, "root2")$volumes-mean(FindNode(hvolsP43StandardMales, "Somatomotor areas")$volumes/FindNode(hvolsP43StandardMales, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Somatomotor areas")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = SomatoMotorPerinat, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteP) +
ylim(-3,3)+
theme_classic(base_size=14)

SomatoMotorAdult <- dataAdult %>%
mutate(SomatoMotorAdult = (FindNode(hvolsAdult, "Somatomotor areas")$volumes/FindNode(hvolsAdult, "root2")$volumes-mean(FindNode(hvolsP96Standard, "Somatomotor areas")$volumes/FindNode(hvolsP96Standard, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Somatomotor areas")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = SomatoMotorAdult, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteA) +
ylim(-3,3) +
theme_classic(base_size=14)

png(here("figures", "fig1D_roi_zscores_somatomotor.png"),
    units = "cm", width = 14, height = 10, res = 150)
plot_grid(SomatoMotorPerinat, SomatoMotorAdult, labels = "AUTO")
dev.off()

# Orbital area
OrbitalPerinat <- mydataTable2PeriMales %>%
mutate(OrbitalPerinat = (FindNode(hvolsP43Males, "Orbital area")$volumes/FindNode(hvolsP43Males, "root2")$volumes-mean(FindNode(hvolsP43StandardMales, "Orbital area")$volumes/FindNode(hvolsP43StandardMales, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Orbital area")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = OrbitalPerinat, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteP) +
ylim(-3,3)+
theme_classic(base_size=14)

OrbitalAdult <- dataAdult %>%
mutate(OrbitalAdult = (FindNode(hvolsAdult, "Orbital area")$volumes/FindNode(hvolsAdult, "root2")$volumes-mean(FindNode(hvolsP96Standard, "Orbital area")$volumes/FindNode(hvolsP96Standard, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Orbital area")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = OrbitalAdult, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteA) +
ylim(-3,3) +
theme_classic(base_size=14)

png(here("figures", "fig1D_roi_zscores_orbital.png"),
    units = "cm", width = 14, height = 10, res = 150)
plot_grid(OrbitalPerinat, OrbitalAdult, labels = "AUTO")
dev.off()

# Auditory areas
AuditPerinat <- mydataTable2PeriMales %>%
mutate(AuditPerinat = (FindNode(hvolsP43Males, "Auditory areas")$volumes/FindNode(hvolsP43Males, "root2")$volumes-mean(FindNode(hvolsP43StandardMales, "Auditory areas")$volumes/FindNode(hvolsP43StandardMales, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Auditory areas")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = AuditPerinat, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteP) +
ylim(-3,3)+
theme_classic(base_size=14)

AuditAdult <- dataAdult %>%
mutate(AuditAdult = (FindNode(hvolsAdult, "Auditory areas")$volumes/FindNode(hvolsAdult, "root2")$volumes-mean(FindNode(hvolsP96Standard, "Auditory areas")$volumes/FindNode(hvolsP96Standard, "root2")$volumes))/sd(FindNode(hvolsP43P96Males, "Auditory areas")$volumes/FindNode(hvolsP43P96Males, "root2")$volumes)) %>%
ggplot() +
aes(x = Housing, y = AuditAdult, colour = Housing) +
geom_point(position = position_jitterdodge(),
alpha=0.7, color="black", shape=21, size=3, aes(fill=factor(Housing)), show.legend = FALSE) +
stat_summary(fun.data = mean_cl_boot,
geom = "pointrange",
position = position_jitterdodge(), color="black", shape=21, size=1, aes(fill=factor(Housing)), show.legend = FALSE) +
scale_fill_manual(values=cbPaletteA) +
ylim(-3,3) +
theme_classic(base_size=14)

png(here("figures", "fig1D_roi_zscores_auditory.png"),
    units = "cm", width = 14, height = 10, res = 150)
plot_grid(AuditPerinat, AuditAdult, labels = "AUTO")
dev.off()

message("fig1_original.R complete — fig1C and fig1D panels saved to figures/")

sessionInfo()
