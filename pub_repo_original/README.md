# Analysis code — Kaller & Ligneul (2025)

*Environmental enrichment during gestation alters neonatal brain structure via maternal care*

---

## Overview

This repository contains the original R code used to produce all figures in the paper. The study examines how perinatal environmental enrichment (EE) affects neonatal (P7) brain structure and whether this effect is mediated by changes in maternal care behaviour. Three imaging datasets are used:

| Label | Age at perfusion | n | Description |
|-------|-----------------|---|-------------|
| N | P7 (neonatal) | 114 | Pups from EE or standard-housing (SH) dams |
| P | P43 (perinatal) | 49 | Animals raised in EE or SH |
| A | P96 (adult) | 25 | Animals switched to EE or SH in adulthood |

---

## Repository structure

```
code/
├── fig1_original.R          Figure 1 — voxelwise EE effects P43/P96, ROI z-scores
├── fig2_original.R          Figure 2 — voxelwise EE effect P7, brain and ROI volumes
├── fig3_original.R          Figure 3 — maternal behaviour, voxelwise contact effect, striatum scatter
├── fig4_original.R          Figure 4 — voxelwise mediation, ACME/ADE bar charts
├── suppfig1_original.R      Suppl Fig 1 — P43 housing effect: all subjects vs males only
├── suppfig2_original.R      Suppl Fig 2 — Effect of age difference at perfusion
├── suppfig3_original.R      Suppl Fig 3 — P7 model components (housing, litter size, litter order)
├── suppfig4_original.R      Suppl Fig 4 — Combined EE models, 4-panel overlay
└── suppfig5_original.R      Suppl Fig 5 — Maternal behaviour panels B–F
```

---

## Two types of code

Each script is self-contained and clearly divided into two sections:

**Local section** — runs on any machine from the CSV/xlsx data files listed below. Produces the data-based figure panels.

**Cluster-only section** — requires [RMINC](https://github.com/Mouse-Imaging-Centre/RMINC), [MRIcrotome](https://github.com/Mouse-Imaging-Centre/MRIcrotome), and access to the `.mnc` neuroimaging files stored on the Oxford BMRC cluster. These sections are included for full transparency of the methods but cannot be re-run without cluster access. The cluster path is defined as `DATA_DIR` at the top of each cluster section.

| Script | Local panels | Cluster-only panels |
|--------|-------------|---------------------|
| `fig1_original.R` | — | Fig 1C (voxelwise), Fig 1D (ROI z-scores) |
| `fig2_original.R` | Fig 2C (brain volume) | Fig 2A (voxelwise), Fig 2C (atlas ROIs) |
| `fig3_original.R` | Fig 3B, 3C (behaviour), Fig 3F (scatter) | Fig 3D (overlays), Fig 3E (t-stat scatter) |
| `fig4_original.R` | Fig 4C, 4D (bar charts) | Fig 4A, 4B (mediation overlays) |
| `suppfig1_original.R` | — | 2-panel P43 housing overlay |
| `suppfig2_original.R` | — | Single age-difference overlay |
| `suppfig3_original.R` | — | 3-panel P7 model components |
| `suppfig4_original.R` | — | 4-panel combined EE models |
| `suppfig5_original.R` | Panels B–F | — |

---

## Data

The input data files are not included in this repository. The local sections of the scripts require the following files in `data/derived/`:

| File | Used by |
|------|---------|
| `MatBehav_per_Mother_Malte.csv` | fig3, suppfig5 |
| `MaternalCareAndVolume_1_.xlsx` | fig2, fig3, suppfig5 |
| `Project_3_Data_with_MaternalCare_IDsorted_VoxWiseBMRC.csv` | fig3 (cluster) |

The voxelwise neuroimaging files (`.mnc` Jacobian maps, atlas files, and masks) are stored on the Oxford BMRC cluster and referenced by the absolute paths defined in each script's cluster section.

---

## Key models

**P7 voxelwise (dataset N, Fig 2A, Suppl Fig 3):**
```
HighB.Distortion.Corrected ~ Housing + Litter + Size_Litter + Sex
```
FDR 1%, t ≈ 3.23 for Housing.

**P43 and P96 voxelwise (datasets P and A, Fig 1C, Suppl Figs 1–4):**
```
Relative_Jacobians ~ Housing + Sex + Group          # common EE effect
Relative_Jacobians ~ Housing * Group + Sex          # perinatal-specific interaction
```

**Voxelwise mediation (dataset N, Fig 4):**
```
Mediator : Total_Contact ~ Housing + Size_Litter
Outcome  : voxel ~ Housing + Total_Contact + Size_Litter
```
5 000 quasi-Bayesian simulations via `mediation::mediate()`, parallelised with `RMINC::pMincApply()` (500 batches, seed 12345).

**ROI z-score formula (Figs 1D, 2C):**
```
z = (roi/brain − mean_SH_males) / SD_all_males
```
Males only; SD shared across P43 and P96 datasets to allow cross-age comparison.

---

## Software

**Local analysis:**
```r
install.packages(c("here", "readxl", "dplyr", "ggplot2",
                   "cowplot", "effsize", "Hmisc", "mediation"))
```

**Cluster analysis (additional):**
- [RMINC](https://github.com/Mouse-Imaging-Centre/RMINC) ≥ 1.5
- [MRIcrotome](https://github.com/Mouse-Imaging-Centre/MRIcrotome)
- `data.tree`

R version: 4.4.x

All scripts use `here::here()` for paths. Open the `.Rproj` file in RStudio before running, or set the working directory to the project root manually.

---

## Citation

Kaller M & Ligneul R (2025). *Environmental enrichment during gestation alters neonatal brain structure via maternal care.* [journal TBC]
