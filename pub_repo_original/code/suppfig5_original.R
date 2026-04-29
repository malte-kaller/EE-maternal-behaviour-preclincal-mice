# suppfig5_original.R
# Supplementary Figure 5: Maternal behaviour panels B–F
#
# Panel A (study design schematic) cannot be reproduced from code.
#
# Source: Paper_PlottingMaternalCare.Rmd (panels B, C)
#         MatBehav_per_Mother_Malte.csv (panels D, E)
#         MaternalCareAndVolume_1_.xlsx (panel F)
#
# LOCAL — all panels reproduced from CSV/xlsx data.
#
# Outputs:
#   figures/suppfig5B_total_contact_all_cages.png
#   figures/suppfig5C_active_nursing_all_cages.png
#   figures/suppfig5D_behavioral_categories.png
#   figures/suppfig5E_contact_vs_littersize.png
#   figures/suppfig5F_brainstem_contact_scatter.png

library(here)
library(ggplot2)
library(dplyr)
library(readxl)

# ─────────────────────────────────────────────────────────────────────────────
# Load per-dam summary data
# ─────────────────────────────────────────────────────────────────────────────

df_f <- read.csv(here("data", "derived", "MatBehav_per_Mother_Malte.csv"))

df_f$Condition <- factor(df_f$Condition)
df_f$Condition <- relevel(df_f$Condition, ref = "SH")
df_f$Litter <- factor(df_f$Litter)
df_f$Litter <- relevel(df_f$Litter, ref = "Litter 1")
df_f$Cage <- factor(df_f$Cage)

df_f_bold  <- subset(df_f, Cage %in% c("CageB", "CageC"))
df_f_other <- subset(df_f, Cage %in% c("CageA", "CageD", "CageE"))

# ─────────────────────────────────────────────────────────────────────────────
# Panel B: Total contact time — all dams, colored by cage
# Source: Paper_PlottingMaternalCare.Rmd (plot_TC2)
# ─────────────────────────────────────────────────────────────────────────────

plot_TC2 = ggplot(df_f, aes(x = Litter, y = Total_Contact, group = Cage)) +
  geom_line(data = df_f_other, aes(linetype = "solid", color = Cage)) +
  geom_line(data = df_f_bold,  aes(linetype = "longdash", color = Cage), size = 1) +
  geom_point(aes(color = Cage, shape = Condition), size = 4) +
  theme(legend.position = "none", aspect.ratio = 1, text = element_text(size = 12)) +
  theme_classic() +
  labs(x = "Litter", y = "Total Contact time") +
  guides(linetype = "none")

ggsave(here("figures", "suppfig5B_total_contact_all_cages.png"),
       plot_TC2, height = 8, width = 10, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel C: Active nursing time — all dams, colored by cage
# Source: Paper_PlottingMaternalCare.Rmd (plot_AN2)
# ─────────────────────────────────────────────────────────────────────────────

plot_AN2 = ggplot(df_f, aes(x = Litter, y = ContactTime_Active_Nurse, group = Cage)) +
  geom_line(data = df_f_other, aes(linetype = "solid", color = Cage)) +
  geom_line(data = df_f_bold,  aes(linetype = "longdash", color = Cage), size = 1.5) +
  geom_point(aes(color = Cage, shape = Condition), size = 4) +
  theme(legend.position = "none", aspect.ratio = 1, text = element_text(size = 12)) +
  theme_classic() +
  labs(x = "Litter", y = "Active nursing time") +
  guides(linetype = "none")

ggsave(here("figures", "suppfig5C_active_nursing_all_cages.png"),
       plot_AN2, height = 8, width = 10, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel D: Behavioral sub-categories by condition
# Columns from MatBehav_per_Mother_Malte.csv
# ─────────────────────────────────────────────────────────────────────────────

behav_cols <- c(
  "ContactTime_Inactive Nurse",
  "ContactTime_Passive contact",
  "ContactTime_Nurse/Groom",
  "ContactTime_Building"
)

behav_labels <- c("Inactive nursing time (sec)", "Passive Contact",
                   "Nurse and Groom", "Building")

make_behav_plot <- function(col, label) {
  ggplot(df_f, aes_string(x = "Condition", y = paste0("`", col, "`"), group = 1)) +
    geom_jitter(aes(shape = Condition), position = position_jitter(0.2), size = 4) +
    stat_summary(fun = mean, geom = "crossbar", size = 0.3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
    scale_shape_manual(values = c(16, 17)) +
    theme_classic() +
    theme(aspect.ratio = 2, legend.position = "none", text = element_text(size = 12)) +
    labs(x = "Condition", y = label)
}

library(cowplot)
behav_plots <- mapply(make_behav_plot, behav_cols, behav_labels, SIMPLIFY = FALSE)
suppfig5D <- plot_grid(plotlist = behav_plots, nrow = 1)

ggsave(here("figures", "suppfig5D_behavioral_categories.png"),
       suppfig5D, height = 8, width = 20, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel E: Maternal contact time vs litter size, colored by cage
# ─────────────────────────────────────────────────────────────────────────────

plot_contact_littersize <- ggplot(df_f, aes(x = LitterSize, y = Total_Contact, color = Cage)) +
  geom_point(aes(shape = Condition), size = 4) +
  geom_smooth(aes(group = Litter), method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  theme_classic() +
  theme(text = element_text(size = 12)) +
  labs(x = "Litter Size", y = "Maternal contact time")

ggsave(here("figures", "suppfig5E_contact_vs_littersize.png"),
       plot_contact_littersize, height = 8, width = 10, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel F: Brainstem z-score vs total contact time
# Data: MaternalCareAndVolume_1_.xlsx
# ─────────────────────────────────────────────────────────────────────────────

dataMCV <- read_excel(here("data", "derived", "MaternalCareAndVolume_1_.xlsx"))
names(dataMCV) <- tolower(gsub("\\.", "_", names(dataMCV)))

roi_df <- dataMCV |>
  filter(image_quality == "Good") |>
  mutate(
    Housing       = factor(housing, levels = c("SH", "EE")),
    brainstem     = medulla + pons + midbrain,
    brainstem_norm = brainstem / volumes
  )

zscore_n <- function(x, housing, ref = "SH") {
  ref_mean <- mean(x[housing == ref], na.rm = TRUE)
  all_sd   <- sd(x, na.rm = TRUE)
  (x - ref_mean) / all_sd
}

roi_df <- roi_df |>
  mutate(brainstem_z = zscore_n(brainstem_norm, housing))

cbPalette <- c("#EDD7C9", "#B76029")

r_val <- cor(roi_df$total_contact, roi_df$brainstem_z, use = "complete.obs")
cat("Brainstem z vs total contact, r =", round(r_val, 2), "\n")

suppfig5F <- ggplot(roi_df, aes(x = total_contact, y = brainstem_z)) +
  geom_point(aes(fill = Housing), alpha = 0.7, color = "black", shape = 21, size = 4,
             show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  scale_fill_manual(values = cbPalette) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(r_val, 2)),
           hjust = 1.2, vjust = 1.5, size = 4) +
  ylim(-3, 3) +
  theme_classic(base_size = 14) +
  labs(x = "Maternal contact time (s)",
       y = "Brainstem volume (z-score)")

ggsave(here("figures", "suppfig5F_brainstem_contact_scatter.png"),
       suppfig5F, width = 8, height = 8, units = "cm", dpi = 300)

message("suppfig5_original.R complete")

sessionInfo()
