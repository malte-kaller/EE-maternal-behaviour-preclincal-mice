# suppfig5_clean.R
# Supplementary Figure 5: Maternal behaviour panels B–F
#
# Panel A (study design schematic) cannot be reproduced from code.
#
# LOCAL — all panels reproduced from CSV/xlsx data.
#
# Data:
#   MatBehav_per_Mother_Malte.csv   — per-dam summary (panels B, C, D, E)
#   MaternalCareAndVolume_1_.xlsx   — P7 pup ROI volumes + contact time (panel F)
#
# Note: Total_Contact is stored in seconds; panels B, C, E convert to minutes.
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
library(cowplot)

CB_PALETTE <- c("#EDD7C9", "#B76029")

# ─────────────────────────────────────────────────────────────────────────────
# Load per-dam summary data
# ─────────────────────────────────────────────────────────────────────────────

df_f <- read.csv(here("data", "derived", "MatBehav_per_Mother_Malte.csv")) |>
  mutate(
    Condition = relevel(factor(Condition), ref = "SH"),
    Litter    = relevel(factor(Litter),    ref = "Litter 1"),
    Cage      = factor(Cage),
    Total_Contact_min        = Total_Contact / 60,
    ContactTime_Active_Nurse_min = ContactTime_Active_Nurse / 60
  )

df_f_bold  <- filter(df_f, Cage %in% c("CageB", "CageC"))
df_f_other <- filter(df_f, Cage %in% c("CageA", "CageD", "CageE"))

# ─────────────────────────────────────────────────────────────────────────────
# Panel B: Total contact time — all dams, colored by cage
# Crossover dams (CageB, CageC) drawn with long-dash; others solid
# ─────────────────────────────────────────────────────────────────────────────

suppfig5B <- ggplot(df_f, aes(x = Litter, y = Total_Contact_min, group = Cage)) +
  geom_line(data = df_f_other, aes(color = Cage), linetype = "solid",   linewidth = 0.7) +
  geom_line(data = df_f_bold,  aes(color = Cage), linetype = "longdash", linewidth = 1.0) +
  geom_point(aes(color = Cage, shape = Condition), size = 4) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1) +
  labs(x = "Litter", y = "Maternal contact time (min)")

ggsave(here("figures", "suppfig5B_total_contact_all_cages.png"),
       suppfig5B, height = 8, width = 10, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel C: Active nursing time — all dams, colored by cage
# ─────────────────────────────────────────────────────────────────────────────

suppfig5C <- ggplot(df_f, aes(x = Litter, y = ContactTime_Active_Nurse_min, group = Cage)) +
  geom_line(data = df_f_other, aes(color = Cage), linetype = "solid",   linewidth = 0.7) +
  geom_line(data = df_f_bold,  aes(color = Cage), linetype = "longdash", linewidth = 1.0) +
  geom_point(aes(color = Cage, shape = Condition), size = 4) +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1) +
  labs(x = "Litter", y = "Active nursing time (min)")

ggsave(here("figures", "suppfig5C_active_nursing_all_cages.png"),
       suppfig5C, height = 8, width = 10, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel D: Behavioral sub-categories by condition (jitter + mean ± SE)
# ─────────────────────────────────────────────────────────────────────────────

behav_vars <- list(
  list(col = "ContactTime_Inactive Nurse", label = "Inactive nursing time (sec)"),
  list(col = "ContactTime_Passive contact", label = "Passive contact (sec)"),
  list(col = "ContactTime_Nurse/Groom",    label = "Nurse and Groom (sec)"),
  list(col = "ContactTime_Building",        label = "Building (sec)")
)

make_behav_plot <- function(col, label) {
  ggplot(df_f, aes(x = Condition, y = .data[[col]], group = 1)) +
    geom_jitter(aes(shape = Condition), position = position_jitter(0.2), size = 4) +
    stat_summary(fun = mean, geom = "crossbar", linewidth = 0.3) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
    scale_shape_manual(values = c(SH = 16, EE = 17)) +
    theme_classic(base_size = 12) +
    theme(aspect.ratio = 2, legend.position = "none") +
    labs(x = "Condition", y = label)
}

behav_plots <- lapply(behav_vars, function(v) make_behav_plot(v$col, v$label))
suppfig5D   <- plot_grid(plotlist = behav_plots, nrow = 1)

ggsave(here("figures", "suppfig5D_behavioral_categories.png"),
       suppfig5D, height = 8, width = 20, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel E: Maternal contact time vs litter size, colored by cage
# Separate regression line per litter
# ─────────────────────────────────────────────────────────────────────────────

suppfig5E <- ggplot(df_f, aes(x = LitterSize, y = Total_Contact_min, color = Cage)) +
  geom_point(aes(shape = Condition), size = 4) +
  geom_smooth(aes(group = Litter), method = "lm", se = TRUE,
              color = "black", linewidth = 0.7) +
  theme_classic(base_size = 12) +
  labs(x = "Litter size", y = "Maternal contact time (min)")

ggsave(here("figures", "suppfig5E_contact_vs_littersize.png"),
       suppfig5E, height = 8, width = 10, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────────────────────────
# Panel F: Brainstem z-score vs total contact time (scatter only)
# Data: MaternalCareAndVolume_1_.xlsx
# Note: contact time converted to minutes for x-axis
# ─────────────────────────────────────────────────────────────────────────────

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
    Housing        = factor(housing, levels = c("SH", "EE")),
    brainstem_norm = (medulla + pons + midbrain) / volumes,
    brainstem_z    = zscore_ref(brainstem_norm, housing),
    total_contact_min = total_contact / 60
  )

r_val <- cor(roi_df$total_contact_min, roi_df$brainstem_z, use = "complete.obs")
cat("Brainstem z vs total contact, r =", round(r_val, 2), "\n")

suppfig5F <- ggplot(roi_df, aes(x = total_contact_min, y = brainstem_z)) +
  geom_point(aes(fill = Housing), alpha = 0.7, color = "black", shape = 21, size = 4,
             show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.7) +
  scale_fill_manual(values = CB_PALETTE) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("r = ", round(r_val, 2)),
           hjust = 1.2, vjust = 1.5, size = 4) +
  ylim(-3, 3) +
  theme_classic(base_size = 14) +
  labs(x = "Maternal contact time (min)",
       y = "Brainstem volume (z-score)")

ggsave(here("figures", "suppfig5F_brainstem_contact_scatter.png"),
       suppfig5F, width = 8, height = 8, units = "cm", dpi = 300)

message("suppfig5_clean.R complete — panels B–F saved to figures/")

sessionInfo()
