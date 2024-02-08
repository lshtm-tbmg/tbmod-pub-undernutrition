# --------------------------------
# Plots for Ukraine / Undernutrition paper
# Rebecca Clark
# Updated 4 February 2024
# -------------------------------- 

rm(list=ls())

suppressPackageStartupMessages({
  library(rlang)
  library(fs)
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(qs)
  library(uuid)
  library(gridExtra)
  library(ggpubr)
  library(digest)
  library(patchwork)
  library(dplyr, warn.conflicts = FALSE)
  options(dplyr.summarise.inform = FALSE)
  
  theme_set(theme_minimal_grid() + panel_border(color = "black"))
  
})

# ------------------------------
# FIGURE 1 FOR MAIN TEXT
# ------------------------------
total_rel <- fread("./epi_output_UKR/rate_differences_UKR_BMI.csv")

total_rel$Scenario <- factor(total_rel$Scenario,
                              levels = c("Export Restriction Shock", "Energy Price Shock",
                                         "Export and Energy Shocks"),
                              labels = c("Export Restriction Shock", "Energy Price Shock",
                                         "Export and Energy Shocks"))

total_rel$BMI <- factor(total_rel$BMI,
                         levels = c("Moderate to severe thinness", "Mild thinness",
                                    "Normal BMI", "Overweight BMI"),
                         labels = c("Moderate to severe thinness", "Mild thinness",
                                    "Normal BMI", "Overweight BMI"))


ggplot(data = total_rel, aes(x = Scenario, y = medval*-100, fill = BMI)) +
  geom_bar(stat="identity", width=.5, position = "dodge") +
  geom_errorbar(aes(x = Scenario, ymin = highval*-100, ymax = lowval*-100),
                position = position_dodge(width = .5), lwd = 0.5, width = 0.25) +
  scale_colour_viridis_d(direction = -1, option = "viridis") +
  scale_fill_viridis_d(direction = -1, option = "viridis") +
  ylab("Relative difference (%)") +
  facet_grid(Indicator ~ .) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.justification = "right",
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 1, size = 14),
    axis.text = element_text(size = 14)
  ) 

# ------------------------------


### Appendix Figures

# ------------------------------
## APPENDIX Figure S5.1
# ------------------------------
rate_plot <- fread("./epi_output_UKR/raw_output_undernutrition_UKR.csv")

rate_plot$Scenario <- factor(rate_plot$Scenario,
                              levels = c("No Shocks", "Export Restriction Shock",
                                         "Energy Price Shock", "Export and Energy Shocks"),
                              labels = c("No Shocks", "Export Restriction Shock",
                                         "Energy Price Shock", "Export and Energy Shocks"))

ggplot(data = rate_plot) +
  geom_ribbon(aes(x = Year, ymin = lowval, ymax = highval, fill = Scenario), alpha = 0.25) +
  geom_line(aes(x = Year, y = medval, col = Scenario)) +
  scale_colour_viridis_d(direction = -1, option = "viridis") +
  scale_fill_viridis_d(direction = -1, option = "viridis") +
  ylab("Rate (per 100,000)") + ylim(c(0,NA)) +
  facet_wrap(~ Indicator, scales = "free_y") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.justification = "right",
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 1, size = 14),
    axis.text = element_text(size = 14)
  ) 

# ---------------



# ---------------------------------------
# APPENDIX: FIGURE S5.2 and S5.3 by BMI
# ---------------------------------------

total_plot <- fread("./epi_output_UKR/raw_output_undernutrition_UKR_BMI.csv")

total_plot$BMI <- factor(total_plot$BMI,
                         levels = c("Moderate to severe thinness", "Mild thinness",
                                    "Normal BMI", "Overweight BMI"),
                         labels = c("Moderate to severe thinness", "Mild thinness",
                                    "Normal BMI", "Overweight BMI"))

total_plot$Scenario <- factor(total_plot$Scenario,
                              levels = c("No Shocks", "Export Restriction Shock",
                                         "Energy Price Shock", "Export and Energy Shocks"),
                              labels = c("No Shocks", "Export Restriction Shock",
                                         "Energy Price Shock", "Export and Energy Shocks"))

ggplot(data = total_plot[(Indicator == "TB incidence rate")]) +
  geom_ribbon(aes(x = Year, ymin = lowval, ymax = highval, fill = Scenario), alpha = 0.25) +
  geom_line(aes(x = Year, y = medval, col = Scenario)) +
  scale_colour_viridis_d(direction = -1, option = "viridis") +
  scale_fill_viridis_d(direction = -1, option = "viridis") +
  ylab("Rate (per 100,000)") + ylim(c(0,NA)) +
  facet_wrap(BMI ~ Indicator, scales = "free_y") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.justification = "right",
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 1, size = 14),
    axis.text = element_text(size = 14)
  ) 

ggplot(data = total_plot[(Indicator == "Number of TB cases")]) +
  geom_ribbon(aes(x = Year, ymin = lowval, ymax = highval, fill = Scenario), alpha = 0.25) +
  geom_line(aes(x = Year, y = medval, col = Scenario)) +
  scale_colour_viridis_d(direction = -1, option = "viridis") +
  scale_fill_viridis_d(direction = -1, option = "viridis") +
  ylab("Number (1000s)") + ylim(c(0,NA)) +
  facet_wrap(BMI ~ Indicator, scales = "free_y") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.justification = "right",
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_text(hjust = 1, size = 14),
    axis.text = element_text(size = 14)
  ) 


# ------ end


