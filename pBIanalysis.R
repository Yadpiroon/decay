# pBI143 and temperature 4, 15, 25 c#

##### Load necessary libraries####
library(readxl)
library(tidyverse)
library(gapminder)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(corrplot)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)
library(ggsignif)
library(readxl)
library(readxl)
library(patchwork)
library(cowplot)
library(lattice)
library(PASWR)
library(ggpattern)
library(rlang)

####Fig 1 bar####
library(readxl)
tcecoc <- read_excel("tcecoc.xlsx", col_types = c("text", "numeric", "numeric", "numeric"))
head(tcecoc)

df <- tcecoc %>%
  rename(
    Location = site,
    `Total Coliforms` = total,
    `E. coli` = ecoli,
    `Other Coliforms` = other
  ) %>%
  pivot_longer(cols = -Location, names_to = "Bacteria", values_to = "Concentration")

df$Location <- factor(df$Location, levels = c("Upstream", "Downstream", "Wastewater"))
comparisons <- list(c("Upstream", "Wastewater"), c("Downstream", "Wastewater"))

df_summary <- df %>%
  group_by(Location, Bacteria) %>%
  summarise(
    mean_conc = mean(Concentration, na.rm = TRUE),
    sd_conc = sd(Concentration, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  )

df_summary$Bacteria <- factor(df_summary$Bacteria,
                              levels = c("Total Coliforms", "E. coli", "Other Coliforms"))

comparisons <- list(c("Upstream", "Wastewater"), c("Downstream", "Wastewater"))

bar <- ggplot(df_summary, aes(x = Location, y = mean_conc, fill = Bacteria)) +
    geom_bar_pattern(
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    alpha = 0.5,
    pattern_density = 0.3,
    pattern_spacing = 0.02,
    pattern_key_scale_factor = 0.6,
    aes(pattern = Bacteria, fill = Bacteria)
  ) +
  scale_pattern_manual(values = c(
    "stripe",     # Total Coliforms
    "crosshatch", # E. coli
    "circle"      # Other Coliforms
  )) +
  scale_fill_manual(
    values = c("lightpink", "skyblue", "purple"),
    labels = c("Total Coliforms", "E. coli", "Other Coliforms")
  ) +
  geom_errorbar(aes(ymin = mean_conc - sd_conc, ymax = mean_conc + sd_conc),
                position = position_dodge(width = 0.8), width = 0.3) +
  geom_text(aes(label = paste0("n=", n), y = 0),
            position = position_dodge(width = 0.8),
            vjust = 1.5,
            size = 3.5, color = "black") +
  labs(
    y = expression(Concentration~(log[10]~CFU/mL)),
    x = NULL
  ) +

  theme_minimal() +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_signif(
    comparisons = comparisons,
    test = "t.test",
    map_signif_level = function(p) {
      if (p < 0.05) {
        return("*")
      } else {
        return(sprintf("p = %.3f", p))
      }
    },
    textsize = 4,
    y_position = c(6, 7)
  )

bar
ggsave(file="barfig1.jpeg", bar, width= 150, height = 120, units = "mm", dpi=800)

### figure 1 bar and mean####
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

tcecoc <- read_excel("tcecoc.xlsx", col_types = c("text",
                                                  "numeric", "numeric", "numeric", "numeric"))
head(tcecoc)

df2 <- tcecoc %>%
  rename(
    Location = site,
    `Total Coliforms` = total,
    `E. coli` = ecoli
  ) %>%
  pivot_longer(cols = c(`Total Coliforms`, `E. coli`, sfmD), names_to = "Bacteria", values_to = "Concentration")

df2$Location <- factor(df2$Location, levels = c("Upstream", "Downstream", "Wastewater"))

comparisons <- list(c("Upstream", "Wastewater"), c("Downstream", "Wastewater", c("Upstream", "Downstream")))

Meanbar <- ggplot(df2, aes(x = Location, y = Concentration, color = Bacteria, shape = Bacteria)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = 0.5)) +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    map_signif_level = function(p) if (p < 0.05) "*" else sprintf("p = %.3f", p),
    y_position = c(7, 8),
    textsize = 4,
    aes(group = Bacteria),
    color = "black"
  ) +
  scale_color_manual(values = c("deeppink", "steelblue", "blue")) +
  labs(
    y = expression(Concentration~(log[10]~CFU/mL)~or~(log[10]~copies/L)),
    x = NULL
  ) +
  theme_minimal() +
  theme(
    axis.ticks.length = unit(0.20, "cm"),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )
Meanbar

ggsave(file="Meanbarfig1.jpeg", Meanbar, width= 150, height = 120, units = "mm", dpi=800)

#### fig 2 tc ec ####
library(readxl)
lnww <- read_excel("lnww.xlsx", col_types = c("text",
                                              "text", "numeric", "numeric"))
head(lnww)

library(readxl)
fig2newgrtcec <- read_excel("fig2newgrtcec.xlsx",
                            col_types = c("text", "text", "numeric",
                                          "numeric", "numeric", "numeric",
                                          "numeric"))
head(fig2newgrtcec)

wastewater <- fig2newgrtcec[46:129, ]
dw <- fig2newgrtcec %>% filter(Type == "Downstream")
up <- fig2newgrtcec %>% filter(Type == "Upstream")

head(wastewater)
head(dw)
head(up)

wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
TCww <- ggplot(wastewater, aes(x = Days, y = TC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Total coliforms in wastewater",
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

TCww
ggsave(file="TCww.jpeg", TCww, width= 150, height = 150, units = "mm", dpi=800)

ECww <- ggplot(wastewater, aes(x = Days, y = EC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = expression(italic("E. coli")~"in wastewater"),
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

ECww
ggsave(file="ECww.jpeg", ECww, width= 150, height = 150, units = "mm", dpi=800)

sfmDww <- ggplot(wastewater, aes(x = Days, y = sfmD, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = expression(italic("sfmD")~"in wastewater"),
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~copies/L))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

sfmDww
ggsave(file="sfmDww.jpeg", sfmDww, width= 150, height = 150, units = "mm", dpi=800)


OCww <- ggplot(wastewater, aes(x = Days, y = OC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Other coliforms",
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(1, 10, by = 1)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

OCww
ggsave(file="OCww.jpeg", OCww, width= 150, height = 150, units = "mm", dpi=800)

### Com
Fig2ww <- plot_grid(TCww, ECww, OCww, ncol = 1,
                    labels = c("a", "b", "c"), label_size = 12)

ggsave(file="Fig2ww.jpeg", Fig2ww, width= 100, height = 240, units = "mm", dpi=900)

# new combine total other and e. coli

Fig2ww2 <- plot_grid(TCww, OCww, ECww, ncol = 1,
                    labels = c("a", "b", "c"), label_size = 12)

ggsave(file="Fig2ww2.jpeg", Fig2ww2, width= 100, height = 240, units = "mm", dpi=900)

#### UP DW

library(readxl)
lnww <- read_excel("lnww.xlsx", col_types = c("text",
                                              "text", "numeric", "numeric"))
head(lnww)

library(readxl)
fig2newgrtcec <- read_excel("fig2newgrtcec.xlsx",
                            col_types = c("text", "text", "numeric",
                                          "numeric", "numeric", "numeric"))
head(fig2newgrtcec)

wastewater <- fig2newgrtcec[46:129, ]
dw <- fig2newgrtcec %>% filter(Type == "Downstream")
up <- fig2newgrtcec %>% filter(Type == "Upstream")

head(wastewater)
head(dw)
head(up)

wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
dw$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))

TCup <- ggplot(up, aes(x = Days, y = TC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Total coliforms in upstream river",
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

TCup
ggsave(file="TCup.jpeg", TCup, width= 150, height = 150, units = "mm", dpi=800)

ECup <- ggplot(up, aes(x = Days, y = EC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = expression(italic("E. coli")~"in upstream river"),
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

ECup
ggsave(file="ECup.jpeg", ECup, width= 150, height = 150, units = "mm", dpi=800)


sfmDup <- ggplot(up, aes(x = Days, y = sfmD, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = expression(italic("sfmD")~"in upstream river"),
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~copies/L))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

sfmDup
ggsave(file="sfmDup.jpeg", sfmDup, width= 150, height = 150, units = "mm", dpi=800)


OCup <- ggplot(up, aes(x = Days, y = OC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Other coliforms in upstream river",
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 0.5)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

OCup
ggsave(file="OCup.jpeg", OCup, width= 150, height = 150, units = "mm", dpi=800)

TCdw <- ggplot(dw, aes(x = Days, y = TC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Total coliforms in downstream river",
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

TCdw
ggsave(file="TCdw.jpeg", TCdw, width= 150, height = 150, units = "mm", dpi=800)

ECdw <- ggplot(dw, aes(x = Days, y = EC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = expression(italic("E. coli")~"in downstream river"),
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

ECdw
ggsave(file="ECdw.jpeg", ECdw, width= 150, height = 150, units = "mm", dpi=800)


sfmDdw <- ggplot(dw, aes(x = Days, y = sfmD, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = expression(italic("sfmD")~"in downstream river"),
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~copies/L))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 6)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

sfmDdw
ggsave(file="sfmDdw.jpeg", sfmDdw, width= 150, height = 150, units = "mm", dpi=800)


OCdw <- ggplot(dw, aes(x = Days, y = OC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Other coliforms in downstream river",
    x = "Time (Days)",
    y = expression(Concentration~(log[10]~CFU/mL))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(0, 10, by = 0.5)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

OCdw
ggsave(file="OCdw.jpeg", OCdw, width= 150, height = 150, units = "mm", dpi=800)


######Combine p5-p10

Fig3river <- plot_grid(TCup, ECup, OCup, TCdw, ECdw, OCdw, ncol = 2,
                       labels = c("a", "b", "c", "d", "e", "f"), label_size = 12)


ggsave(file="Fig3river.jpeg", Fig3river, width= 200, height = 240, units = "mm", dpi=900)

# new combine tc up dw, oc up dw, ec up dw

Fig3river2 <- plot_grid(TCup, TCdw, OCup, OCdw, ECup, ECdw, ncol = 2,
                       labels = c("a", "b", "c", "d", "e", "f"), label_size = 12)


ggsave(file="Fig3river2.jpeg", Fig3river2, width= 200, height = 240, units = "mm", dpi=900)


#### fig 4 boxplot methods####

library(readxl)
boxplot <- read_excel("boxplot.xlsx", col_types = c("text",
                                                    "numeric"))
head(boxplot)

library(ggplot2)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)

EMF <- boxplot %>% filter(Method == "EMF")
PEG <- boxplot %>% filter(Method == "PEG")

summary(EMF$Concentration)
sd(EMF$Concentration)
summary(PEG$Concentration)
sd(PEG$Concentration)

boxplot_n <- boxplot %>%
  group_by(Method) %>%
  summarise(n = n(), .groups = "drop")

boxplot_n <- boxplot %>%
  group_by(Method) %>%
  summarise(n = n(), .groups = "drop")

bx <- ggplot(boxplot, aes(x = Method, y = Concentration, color = Concentration)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 17,
               position = position_dodge(width = 0.5)) +
  geom_text(data = boxplot_n, aes(x = Method, y = 7.2, label = paste0("n=", n)),
            inherit.aes = FALSE,
            vjust = 0,
            size = 4) +

  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # show dots for each sample

  labs(
    title = NULL,
    y = expression(Concentration~(log[10]~copies/L))
  ) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), limits = c(1, 8)) +
  scale_fill_manual(values = c("skyblue", "lightpink")) +

  theme_minimal() +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5)
  )

bx
ggsave(file="bx.jpeg", bx, width= 80, height = 80, units = "mm", dpi=800)

#### peg and emf new to mean####

library(readxl)
boxplot <- read_excel("boxplot.xlsx", col_types = c("text",
                                                    "numeric"))
head(boxplot)

library(ggplot2)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)

EMF <- boxplot %>% filter(Method == "EMF")
PEG <- boxplot %>% filter(Method == "PEG")

summary(EMF$Concentration)
sd(EMF$Concentration)
summary(PEG$Concentration)
sd(PEG$Concentration)



ggsave(file="Meanbarfig1.jpeg", Meanbar, width= 150, height = 120, units = "mm", dpi=800)


#####presentation####
library(readxl)
graphset1Ln <- read_excel("graphset1Ln.xlsx",
                          col_types = c("text", "text", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric"))
head(graphset1Ln)

wastewater <- graphset1Ln[1:84, ]
river <- graphset1Ln[85:174, ]
dw <- graphset1Ln %>% filter(Type == "Downstream")
up <- graphset1Ln %>% filter(Type == "Upstream")

#####wwpbi####

wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww <- ggplot(wastewater, aes(x = Days, y = pBI, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww
ggsave(file="pBIww.jpeg", pBIww, width= 150, height = 150, units = "mm", dpi=900)

#####uppbi####

up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup <- ggplot(up, aes(x = Days, y = pBI, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in upstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup
ggsave(file="pBIriverup.jpeg", pBIriverup, width= 150, height = 150, units = "mm", dpi=800)

#####dwpbi####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw <- ggplot(dw, aes(x = Days, y = pBI, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in downstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw
ggsave(file="pBIriverdw.jpeg", pBIriverdw, width= 150, height = 150, units = "mm", dpi=800)

#####wwcrass####

wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww1 <- ggplot(wastewater, aes(x = Days, y = crass, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "crAssphage in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww1
ggsave(file="pBIww.jpeg", pBIww1, width= 150, height = 150, units = "mm", dpi=800)
#####upcrass####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup1 <- ggplot(up, aes(x = Days, y = crass, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "crAssphage in upstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup1
ggsave(file="pBIriverup1.jpeg", pBIriverup1, width= 150, height = 150, units = "mm", dpi=800)

#####dwcrass####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw1 <- ggplot(dw, aes(x = Days, y = crass, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "crAssphage in downstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw1
ggsave(file="pBIriverdw1.jpeg", pBIriverdw1, width= 150, height = 150, units = "mm", dpi=800)
#####wwgyrb####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww11 <- ggplot(wastewater, aes(x = Days, y = gyrB, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "gyrB in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww11
ggsave(file="pBIww11.jpeg", pBIww11, width= 150, height = 150, units = "mm", dpi=800)
#####upgyrb####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup11 <- ggplot(up, aes(x = Days, y = gyrB, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "gyrB in upstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup11
ggsave(file="pBIriverup11.jpeg", pBIriverup11, width= 150, height = 150, units = "mm", dpi=800)

#####dwgyrb####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw11 <- ggplot(dw, aes(x = Days, y = gyrB, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "gyrB in downstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw11
ggsave(file="pBIriverdw11.jpeg", pBIriverdw11, width= 150, height = 150, units = "mm", dpi=800)

#####wwhf####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111 <- ggplot(wastewater, aes(x = Days, y = HE183, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "HF183 in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111
ggsave(file="pBIww111.jpeg", pBIww111, width= 150, height = 150, units = "mm", dpi=800)
#####uphf####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111 <- ggplot(up, aes(x = Days, y = HE183, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "HF183 in upstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111
ggsave(file="pBIriverup111.jpeg", pBIriverup111, width= 150, height = 150, units = "mm", dpi=800)

#####dwhf####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111 <- ggplot(dw, aes(x = Days, y = HE183, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "HF183 in downstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111
ggsave(file="pBIriverdw111.jpeg", pBIriverdw111, width= 150, height = 150, units = "mm", dpi=800)

#####wwintl####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111intl <- ggplot(wastewater, aes(x = Days, y = intl, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "intl1 in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111intl
ggsave(file="pBIww111intl.jpeg", pBIww111intl, width= 150, height = 150, units = "mm", dpi=800)
#####upintl####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111intl <- ggplot(up, aes(x = Days, y = intl, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "intl in upstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111intl
ggsave(file="pBIriverup111intl.jpeg", pBIriverup111intl, width= 150, height = 150, units = "mm", dpi=800)

#####dwintl####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111intl <- ggplot(dw, aes(x = Days, y = intl, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "intl in downstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111intl
ggsave(file="pBIriverdw111intl.jpeg", pBIriverdw111intl, width= 150, height = 150, units = "mm", dpi=800)

#####wwsul####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111sul <- ggplot(wastewater, aes(x = Days, y = sul, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "sul1 in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111sul
ggsave(file="pBIww111sul.jpeg", pBIww111sul, width= 150, height = 150, units = "mm", dpi=800)
#####upsul####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111sul <- ggplot(up, aes(x = Days, y = sul, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "sul1 in upstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111sul
ggsave(file="pBIriverup111sul.jpeg", pBIriverup111sul, width= 150, height = 150, units = "mm", dpi=800)

#####dwsul####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111sul <- ggplot(dw, aes(x = Days, y = sul, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "sul1 in downstream river",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111sul
ggsave(file="pBIriverdw111sul.jpeg", pBIriverdw111sul, width= 150, height = 150, units = "mm", dpi=800)

#####wwinvA####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111invA <- ggplot(wastewater, aes(x = Days, y = invA, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "invA",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111invA
ggsave(file="pBIww111invA.jpeg", pBIww111invA, width= 150, height = 150, units = "mm", dpi=800)
#####upinvA####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111invA <- ggplot(up, aes(x = Days, y = invA, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "invA",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111invA
ggsave(file="pBIriverup111invA.jpeg", pBIriverup111invA, width= 150, height = 150, units = "mm", dpi=800)

#####dwinvA####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111invA <- ggplot(dw, aes(x = Days, y = invA, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "invA",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111invA
ggsave(file="pBIriverdw111invA.jpeg", pBIriverdw111invA, width= 150, height = 150, units = "mm", dpi=800)

#####wwsfmD####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111sfmD <- ggplot(wastewater, aes(x = Days, y = sfmD, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "sfmD",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111sfmD
ggsave(file="pBIww111sfmD.jpeg", pBIww111sfmD, width= 150, height = 150, units = "mm", dpi=800)
#####upsfmD####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111sfmD <- ggplot(up, aes(x = Days, y = sfmD, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "sfmD",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111sfmD
ggsave(file="pBIriverup111sfmD.jpeg", pBIriverup111sfmD, width= 150, height = 150, units = "mm", dpi=800)

#####dwsfmD####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111sfmD <- ggplot(dw, aes(x = Days, y = sfmD, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "sfmD",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111sfmD
ggsave(file="pBIriverdw111sfmD.jpeg", pBIriverdw111sfmD, width= 150, height = 150, units = "mm", dpi=800)

#####wwHAdVs####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111HAdVs <- ggplot(wastewater, aes(x = Days, y = HAdVs, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "HAdVs",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111HAdVs
ggsave(file="pBIww111HAdVs.jpeg", pBIww111HAdVs, width= 150, height = 150, units = "mm", dpi=800)
#####upHAdVs####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111HAdVs <- ggplot(up, aes(x = Days, y = HAdVs, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "HAdVs",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111HAdVs
ggsave(file="pBIriverup111HAdVs.jpeg", pBIriverup111HAdVs, width= 150, height = 150, units = "mm", dpi=800)

#####dwHAdVs####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111HAdVs <- ggplot(dw, aes(x = Days, y = HAdVs, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "HAdVs",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111HAdVs
ggsave(file="pBIriverdw111HAdVs.jpeg", pBIriverdw111HAdVs, width= 150, height = 150, units = "mm", dpi=800)


#####wwNoVGII####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111NoVGII <- ggplot(wastewater, aes(x = Days, y = NoVGII, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "NoVGII",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111NoVGII
ggsave(file="pBIww111NoVGII.jpeg", pBIww111NoVGII, width= 150, height = 150, units = "mm", dpi=800)
#####upNoVGII####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111NoVGII <- ggplot(up, aes(x = Days, y = NoVGII, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "NoVGII",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111NoVGII
ggsave(file="pBIriverup111NoVGII.jpeg", pBIriverup111NoVGII, width= 150, height = 150, units = "mm", dpi=800)

#####dwNoVGII####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111NoVGII <- ggplot(dw, aes(x = Days, y = NoVGII, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "NoVGII",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111NoVGII
ggsave(file="pBIriverdw111NoVGII.jpeg", pBIriverdw111NoVGII, width= 150, height = 150, units = "mm", dpi=800)

#####wwTC####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111TC <- ggplot(wastewater, aes(x = Days, y = TC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "TC",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111TC
ggsave(file="pBIww111TC.jpeg", pBIww111TC, width= 150, height = 150, units = "mm", dpi=800)
#####upTC####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111TC <- ggplot(up, aes(x = Days, y = TC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "TC",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111TC
ggsave(file="pBIriverup111TC.jpeg", pBIriverup111TC, width= 150, height = 150, units = "mm", dpi=800)

#####dwTC####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111TC <- ggplot(dw, aes(x = Days, y = TC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "TC",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111TC
ggsave(file="pBIriverdw111TC.jpeg", pBIriverdw111TC, width= 150, height = 150, units = "mm", dpi=800)


#####wwEC####
wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww111EC <- ggplot(wastewater, aes(x = Days, y = EC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "EC",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-20, 0, by = 2), limits = c(-20, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww111EC
ggsave(file="pBIww111EC.jpeg", pBIww111EC, width= 150, height = 150, units = "mm", dpi=800)
#####upEC####
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup111EC <- ggplot(up, aes(x = Days, y = EC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "EC",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup111EC
ggsave(file="pBIriverup111EC.jpeg", pBIriverup111EC, width= 150, height = 150, units = "mm", dpi=800)

#####dwEC####
dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw111EC <- ggplot(dw, aes(x = Days, y = EC, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "EC",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-16, 0, by = 2), limits = c(-12, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "black", size = 14, face = "bold"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    strip.text.y = element_text(color = "black", size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw111EC
ggsave(file="pBIriverdw111EC.jpeg", pBIriverdw111EC, width= 150, height = 150, units = "mm", dpi=800)


######Combine ww up dw####

Fig9w <- plot_grid(pBIww, pBIww1, pBIww11, pBIww111, ncol = 2,
                  labels = c("a", "b", "c", "d"), label_size = 14)

ggsave(file="Fig9w.jpeg", Fig9w, width= 170, height = 120, units = "mm", dpi=900)

Fig9u <- plot_grid(pBIriverup, pBIriverup1, pBIriverup11, pBIriverup111, ncol = 2,
                   labels = c("a", "b", "c", "d"), label_size = 14)

ggsave(file="Fig9u.jpeg", Fig9u, width= 170, height = 120, units = "mm", dpi=900)

Fig9d <- plot_grid(pBIriverdw, pBIriverdw1, pBIriverdw11, pBIriverdw111, ncol = 2,
                   labels = c("a", "b", "c", "d"), label_size = 14)

ggsave(file="Fig9d.jpeg", Fig9d, width= 170, height = 120, units = "mm", dpi=900)


#### Export fig 2 TC/EC/sfmD Temp WW UP DW ####
Fig2TCECSfmD <- plot_grid(TCww, TCdw, TCup,
                          ECww, ECdw, ECup,
                          sfmDww, sfmDdw, sfmDup, ncol = 3,
                   labels = c("a", "b", "c",
                              "d", "e", "f",
                              "g", "h", "i"), label_size = 14)

ggsave(file="Fig2TCECSfmD.jpeg", Fig2TCECSfmD, width= 250, height = 250, units = "mm", dpi=900)

#### Export fig 4 Ln pBI crAss gyrB HF183 Temp WW UP DW ####
Fig4decay<- plot_grid(pBIww, pBIriverdw, pBIriverup,
                      pBIww1, pBIriverdw1, pBIriverup1,
                      pBIww11, pBIriverdw11, pBIriverup11,
                      pBIww111, pBIriverdw111, pBIriverup111,ncol = 3,
                      labels = c("a", "b", "c",
                                 "d", "e", "f",
                                 "g", "h", "i",
                                 "j", "k", "l"), label_size = 14)

ggsave(file="Fig4decay.jpeg", Fig4decay, width= 260, height = 320, units = "mm", dpi=900)

Fig4decay3<- plot_grid(pBIww, pBIriverdw, pBIriverup,
                      pBIww111, pBIriverdw111, pBIriverup111,
                      pBIww11, pBIriverdw11, pBIriverup11,
                      pBIww1, pBIriverdw1, pBIriverup1,ncol = 3,
                      labels = c("a", "b", "c",
                                 "d", "e", "f",
                                 "g", "h", "i",
                                 "j", "k", "l"), label_size = 14)

ggsave(file="Fig4decay3.jpeg", Fig4decay3, width= 260, height = 320, units = "mm", dpi=900)



#### Export fig 5 Ln intl1 and sul Temp WW DW UP ####
Fig5decayARG <- plot_grid(pBIww111intl, pBIriverdw111intl,
                      pBIriverup111intl, pBIww111sul,
                      pBIriverdw111sul, pBIriverup111sul, ncol = 3,
                      labels = c("a", "b",
                                 "c", "d",
                                 "e", "f"), label_size = 14)

ggsave(file="Fig5decayARG.jpeg", Fig5decayARG, width= 260, height = 160, units = "mm", dpi=900)


####anova####

library(readxl)
graphset1Ln <- read_excel("graphset1Ln.xlsx",
                          col_types = c("text", "text", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric"))
head(graphset1Ln)

wastewater <- graphset1Ln[1:84, ]
river <- graphset1Ln[85:174, ]
dw <- graphset1Ln %>% filter(Type == "Downstream")
up <- graphset1Ln %>% filter(Type == "Upstream")

parameters <- c("pBI", "crass", "gyrB", "HE183", "intl", "sul")

run_anova <- function(data, param) {
  formula <- as.formula(paste(param, "~ Temp"))
  model <- aov(formula, data = data)
  summary_res <- summary(model)
  tukey_res <- TukeyHSD(model)
  list(summary = summary_res, tukey = tukey_res)
}

datasets <- list(wastewater = wastewater, river = river, dw = dw, up = up)

results <- list()

for (dataset_name in names(datasets)) {
  cat("\n### Dataset:", dataset_name, "###\n")
  data <- datasets[[dataset_name]]

  for (param in parameters) {
    cat("\n--- Parameter:", param, "---\n")
    res <- run_anova(data, param)
    print(res$summary)
    print(res$tukey)
  }
}

library(readxl)
library(dplyr)
library(writexl)
library(dplyr)

#### export table 7 ####
library(readxl)
library(dplyr)
library(writexl)

# Load your data
graphset1Ln <- read_excel("graphset1Ln.xlsx",
                          col_types = c("text", "text", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric", "numeric", "numeric",
                                        "numeric"))

# Create three datasets
wastewater <- graphset1Ln[1:84, ] %>%
  mutate(Type = "Wastewater")
dw <- graphset1Ln %>% filter(Type == "Downstream")
up <- graphset1Ln %>% filter(Type == "Upstream")

# Combine them into one dataset
combined_data <- bind_rows(wastewater, dw, up)

# Define target parameters
parameters <- c("pBI", "crass", "gyrB", "HE183", "intl", "sul")
temperatures <- unique(combined_data$Temp)

# Store results
all_results <- list()

# Loop through each temperature
for (temp in temperatures) {
  temp_data <- combined_data %>% filter(Temp == temp)

  tukey_results <- list()

  for (param in parameters) {
    # Run ANOVA
    formula <- as.formula(paste(param, "~ Type"))
    model <- aov(formula, data = temp_data)

    # Run Tukey HSD
    tukey <- TukeyHSD(model)
    tukey_df <- as.data.frame(tukey$Type)
    tukey_df$Comparison <- rownames(tukey_df)

    # Keep only comparison and p-value
    tukey_results[[param]] <- tukey_df[, c("Comparison", "p adj")]
    colnames(tukey_results[[param]]) <- c("Comparison", param)
  }

  # Merge all parameter p-values by comparison
  merged_p <- Reduce(function(x, y) full_join(x, y, by = "Comparison"), tukey_results)

  # Add temperature label
  merged_p <- merged_p %>%
    mutate(Temperature = paste0(temp, "C")) %>%
    dplyr::select(Temperature, dplyr::everything())

  # Store in list
  all_results[[as.character(temp)]] <- merged_p
}

# Combine results across all temperatures
final_pvalues <- bind_rows(all_results)

# Export to Excel
write_xlsx(final_pvalues, path = "Tukey_pvalues_WW_DW_UP.xlsx")
cat("Exported: Tukey_pvalues_WW_DW_UP.xlsx\n")


####Fig 5 pbi Ln####
library(readxl)
lnww <- read_excel("lnww.xlsx", col_types = c("text",
                                              "text", "numeric", "numeric"))
View(lnww)

wastewater <- lnww[1:84, ]
river <- lnww[85:174, ]
dw <- lnww %>% filter(Type == "Downstream")
up <- lnww %>% filter(Type == "Upstream")

head(wastewater)
head(river)
head(dw)
head(up)

river$Temp <- factor(river$Temp, levels = c("4C", "15C", "25C"))
pBIriver <- ggplot(river, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 1, aes(group = Temp), linetype = "dashed", alpha = 0.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp)) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 12, face = "bold"),
    legend.position = "top",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriver
ggsave(file="pBIriver.jpeg", pBIriver, width= 150, height = 150, units = "mm", dpi=800)


wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
pBIww <- ggplot(wastewater, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIww
ggsave(file="pBIww.jpeg", pBIww, width= 150, height = 150, units = "mm", dpi=800)

#### generate the tC and EC for Ln####
#TC
library(readxl)
TCLn <- read_excel("TCLn.xls", col_types = c("text", 
                                             "text", "text", "numeric"))
head(TCLn)

wastewater <- TCLn[1:84, ]
river <- TCLn[85:174, ]
dw <- TCLn %>% filter(Type == "Downstream")
up <- TCLn %>% filter(Type == "Upstream")

head(wastewater)
head(river)
head(dw)
head(up)

wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
wastewater$Days <- as.numeric(wastewater$Days)

TCww <- ggplot(wastewater, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Total coliform in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

TCww
ggsave(file="TCww.jpeg", TCww, width= 150, height = 150, units = "mm", dpi=800)


up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
up$Days <- as.numeric(up$Days)
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
TCriverup <- ggplot(up, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Total coliform in upstream river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

TCriverup
ggsave(file="TCriverup.jpeg", TCriverup, width= 150, height = 150, units = "mm", dpi=800)


dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
dw$Days <- as.numeric(dw$Days)

dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
TCriverdw <- ggplot(dw, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "Total coliform in downstream river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

TCriverdw

ggsave(file="TCriverdw.jpeg", TCriverdw, width= 150, height = 150, units = "mm", dpi=800)


####Ecoli Ln####
library(readxl)
ECLn <- read_excel("ECLn.xlsx", col_types = c("text", 
                                              "text", "text", "numeric"))
head(ECLn)

wastewater <- ECLn[1:84, ]
river <- ECLn[85:174, ]
dw <- ECLn %>% filter(Type == "Downstream")
up <- ECLn %>% filter(Type == "Upstream")

head(wastewater)
head(river)
head(dw)
head(up)

wastewater$Temp <- factor(wastewater$Temp, levels = c("4C", "15C", "25C"))
wastewater$Days <- as.numeric(wastewater$Days)

ECww <- ggplot(wastewater, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "E. coli in wastewater",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

ECww
ggsave(file="ECww.jpeg", ECww, width= 150, height = 150, units = "mm", dpi=800)


up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
up$Days <- as.numeric(up$Days)
up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
ECriverup <- ggplot(up, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "E. coli in upstream river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

ECriverup
ggsave(file="ECriverup.jpeg", ECriverup, width= 150, height = 150, units = "mm", dpi=800)


dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
dw$Days <- as.numeric(dw$Days)

dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
ECriverdw <- ggplot(dw, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "E. coli in downstream river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

ECriverdw

ggsave(file="ECriverdw.jpeg", ECriverdw, width= 150, height = 150, units = "mm", dpi=800)

####combine tc and ec####

TCww
TCriverup
TCriverdw

ECww
ECriverup
ECriverdw

Figreview <- plot_grid(TCww, TCriverup, TCriverdw, 
                       ECww, ECriverup, ECriverdw, ncol = 3,
                   labels = c("a", "b", "c", 
                              "d", "e", "f"), label_size = 14)

ggsave(file="Figreview.jpeg", Figreview, width= 260, height = 160, units = "mm", dpi=900)

  ##########

up$Temp <- factor(up$Temp, levels = c("4C", "15C", "25C"))
pBIriverup <- ggplot(up, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in upstream river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverup
ggsave(file="pBIriverup.jpeg", pBIriverup, width= 150, height = 150, units = "mm", dpi=800)


dw$Temp <- factor(dw$Temp, levels = c("4C", "15C", "25C"))
pBIriverdw <- ggplot(dw, aes(x = Days, y = Ln, color = Temp, fill = Temp, shape = Temp, group = Temp)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = Temp, shape = Temp)) +
  stat_summary(fun = mean, geom = "line", size = 0.5, aes(group = Temp), linetype = "dashed", alpha = 0.3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.3, aes(group = Temp), alpha = 1) +
  scale_color_manual(
    values = c("#000000", "#000000", "#000000"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_fill_manual(
    values = c("#2d376c", "#f0f1f1", "#FCB97D"),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    labels = c("4°C", "15°C", "25°C"),
    name = "Temperature"
  ) +
  labs(
    title = "pBI143 in downstream river water",
    x = "Time (Days)",
    y = expression(Ln(C[t]/C[0]))
  ) +
  scale_x_continuous(breaks = seq(0, 14, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 0, by = 2), limits = c(-10, 0)) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )

pBIriverdw

ggsave(file="pBIriverdw.jpeg", pBIriverdw, width= 150, height = 150, units = "mm", dpi=800)

######Combine p5-p10
Fig5 <- plot_grid(pBIww, pBIriverup, pBIriverdw, ncol = 1,
                   labels = c("a", "b", "c"), label_size = 12)

ggsave(file="Fig5.jpeg", Fig5, width= 100, height = 240, units = "mm", dpi=900)

##### Fig 6 corr ####
#####COR RIVER#### #CANT COMBINE DUE TO UNEQUAL N

library(readr)
corr_cpx <- read_csv("corr.csv", col_types = cols(pbi = col_number(),
                                                  tc = col_number(), ec = col_number()))

head(corr_cpx)

conc <- corr_cpx[, 6:8]
head(conc)

summary(conc$pBI143)
sd(conc$pBI143)
shapiro.test(conc$pBI143) #not normal
summary(conc$TC)
sd(conc$TC)
shapiro.test(conc$TC) #not normal
summary(conc$EC)
sd(conc$EC)
shapiro.test(conc$EC) #not normal

corr <- round(cor(conc, method = "spearman"), 2)

p.df <- as.data.frame(ggcorrplot::cor_pmat(conc, method = "spearman"))

labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 ~ "*")
}

p.labs = p.df %>%
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#2d376c", "#f0f1f1", "#FCB97D")) +
  theme(legend.text = element_text(color = "black", size = 10), #detail site label
        legend.title = element_text(color = "black", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", size = 2),
        panel.background = element_rect(fill = "#F1F1F1", colour = NA),
        panel.grid.minor = element_line(colour = "white", size = 0.2),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)


cor.plot.labs

ggsave(file="cor.plot.labsa.jpeg", cor.plot.labs,
       width= 80, height = 60, units = "mm", dpi=600)

library(readr)
corsite <- read_csv("corsite.csv", col_types = cols(pBI143_up = col_number(),
                                                    Total_coliform_up = col_number(), E.coli_up = col_number(),
                                                    pBI143_dw = col_number(), Total_coliform_dw = col_number(),
                                                    E.coli_dw = col_number()))
head(corsite)

corr <- round(cor(corsite, use = "pairwise.complete.obs"), 2)
p_mat <- ggcorrplot::cor_pmat(corsite)
p_mat

labs.function <- function(x) {
  ifelse(x < 0.05, "*", "")
}

p.df <- as.data.frame(p_mat)
p.df$Var1 <- rownames(p.df)
p.df_long <- melt(p.df, id.vars = "Var1", variable.name = "Var2", value.name = "pval")

p.df_long <- p.df_long %>%
  mutate(lab = labs.function(pval)) %>%
  select(Var1, Var2, lab)
cor_plot <- ggcorrplot(corr, hc.order = FALSE, type = "lower", lab = TRUE,
                       ggtheme = theme_gray(),
                       colors = c("#f1c78a", "white", "#9fc995")) +
  theme(
    legend.text = element_text(color = "black", family = "Arial", size = 10),
    legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"),
    axis.text = element_text(color = "black", family = "Arial", size = 10),
    panel.background = element_rect(fill = "#F1F1F1", colour = NA),
    panel.grid.minor = element_line(colour = "white", size = 0.2),
    panel.grid.major = element_line(colour = "white", size = 0.2)
  ) +
  labs(title = "River water")
cor_plot

cor_data <- cor_plot$data
p.labs_final <- inner_join(cor_data, p.df_long, by = c("Var1", "Var2"))

cor.plot.labsa <- cor_plot +
  geom_text(data = p.labs_final,
            aes(x = Var1, y = Var2, label = lab),
            nudge_y = 0.25,
            size = 5)

cor.plot.labsa

ggsave(file="cor.plot.labsa.jpeg", cor.plot.labsa,
       width= 120, height = 100, units = "mm", dpi=600)

# Fig 6 here
getwd()
library(tidyverse)
library(gapminder)
library(ggplot2)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
library(dplyr)
library(survival)
library(corrplot)
library(NADA)
library(NADA2)
library(EnvStats)
library(stats)
library(base)
library(ggsignif)
library(readxl)
library(readxl)
library(patchwork)
library(cowplot)
library(lattice)
library(PASWR)
library(factoextra)
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)

# Load libraries
library(readr)
library(ggplot2)
library(ggcorrplot)
library(dplyr)
library(reshape2)

# Load the data
CORWW <- read_csv("CORWW.csv", col_types = cols(
  pBI143_ww = col_number(),
  Total_coliform_ww = col_number(),
  E.coli_ww = col_number()))
head(CORWW)

# Correlation matrix
corr <- round(cor(CORWW, use = "pairwise.complete.obs"), 2)

# p-value matrix
p_mat <- ggcorrplot::cor_pmat(CORWW)
p_mat
# Function to label significance
labs.function <- function(x) {
  ifelse(x < 0.05, "*", "")
}

# Prepare long-format p-value label table
p.df <- as.data.frame(p_mat)
p.df$Var1 <- rownames(p.df)
p.df_long <- melt(p.df, id.vars = "Var1", variable.name = "Var2", value.name = "pval")

# Add significance labels
p.df_long <- p.df_long %>%
  mutate(lab = labs.function(pval)) %>%
  select(Var1, Var2, lab)

# Create base correlation plot
cor_plot <- ggcorrplot(corr, hc.order = FALSE, type = "lower", lab = TRUE,
                       ggtheme = theme_gray(),
                       colors = c("#3EDBF0", "white", "#FF75A0")) +
  theme(
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 10, face = "bold"),
    axis.text = element_text(color = "black", size = 10),
    panel.background = element_rect(fill = "#F1F1F1", colour = NA),
    panel.grid.minor = element_line(colour = "white", size = 0.2),
    panel.grid.major = element_line(colour = "white", size = 0.2)
  ) +
  labs(title = "Wastewater")
cor_plot
cor_data <- cor_plot$data
p.labs_final <- inner_join(cor_data, p.df_long, by = c("Var1", "Var2"))

cor.plot.labs3 <- cor_plot +
  geom_text(data = p.labs_final,
            aes(x = Var1, y = Var2, label = lab),
            nudge_y = 0.25,
            size = 5)

cor.plot.labs3


ggsave(file="cor.plot.labsa3.jpeg", cor.plot.labs3,
       width= 80, height = 60, units = "mm", dpi=600)

Fig6 <- plot_grid(cor.plot.labsa, cor.plot.labs3, ncol = 1,
                  labels = c("a", "b"), label_size = 12)

ggsave(file="Fig6.jpeg", Fig6, width= 120, height = 190, units = "mm", dpi=900)

####other coliforms


library(readr)
corpc <- read_csv("corpc.csv", col_types = cols(pBI143_up = col_number(),
                                                Total_coliform_up = col_number(), E.coli_up = col_number(),
                                                Other_coliform_up = col_number(), pBI143_dw = col_number(),
                                                Total_coliform_dw = col_number(), E.coli_dw = col_number(),
                                                Other_coliform_dw = col_number()))
head(corpc)

corr <- round(cor(corpc, use = "pairwise.complete.obs"), 2)
p_mat <- ggcorrplot::cor_pmat(corpc)

labs.function <- function(x) {
  ifelse(x < 0.05, "*", "")
}

p.df <- as.data.frame(p_mat)
p.df$Var1 <- rownames(p.df)
p.df_long <- melt(p.df, id.vars = "Var1", variable.name = "Var2", value.name = "pval")

p.df_long <- p.df_long %>%
  mutate(lab = labs.function(pval)) %>%
  select(Var1, Var2, lab)
cor_plot <- ggcorrplot(corr, hc.order = FALSE, type = "lower", lab = TRUE,
                       ggtheme = theme_gray(),
                       colors = c("#f1c78a", "white", "#9fc995")) +
  theme(
    legend.text = element_text(color = "black", family = "Arial", size = 10),
    legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"),
    axis.text = element_text(color = "black", family = "Arial", size = 10),
    panel.background = element_rect(fill = "#F1F1F1", colour = NA),
    panel.grid.minor = element_line(colour = "white", size = 0.2),
    panel.grid.major = element_line(colour = "white", size = 0.2)
  ) +
  labs(title = "River water")

cor_data <- cor_plot$data
p.labs_final <- inner_join(cor_data, p.df_long, by = c("Var1", "Var2"))

cor.plot.labsa1 <- cor_plot +
  geom_text(data = p.labs_final,
            aes(x = Var1, y = Var2, label = lab),
            nudge_y = 0.25,
            size = 5)

cor.plot.labsa1

ggsave(file="cor.plot.labsa1.jpeg", cor.plot.labsa1,
       width= 120, height = 100, units = "mm", dpi=600)


##### cor ww
# Load libraries
library(readr)
library(ggplot2)
library(ggcorrplot)
library(dplyr)
library(reshape2)


library(readr)
CORocWW <- read_csv("CORocWW.csv", col_types = cols(pBI143_ww = col_number(),
                                                    Total_colifrom_ww = col_number(), E.coli_ww = col_number(),
                                                    Other_coliform_ww = col_number()))
head(CORocWW)

# Correlation matrix
corr <- round(cor(CORocWW, use = "pairwise.complete.obs"), 2)

# p-value matrix
p_mat <- ggcorrplot::cor_pmat(CORocWW)
p_mat
# Function to label significance
labs.function <- function(x) {
  ifelse(x < 0.05, "*", "")
}

# Prepare long-format p-value label table
p.df <- as.data.frame(p_mat)
p.df
p.df$Var1 <- rownames(p.df)
p.df$Var1
p.df_long <- melt(p.df, id.vars = "Var1", variable.name = "Var2", value.name = "pval")
p.df_long
# Add significance labels
p.df_long <- p.df_long %>%
  mutate(lab = labs.function(pval)) %>%
  select(Var1, Var2, lab)
p.df_long
# Create base correlation plot
cor_plot <- ggcorrplot(corr, hc.order = FALSE, type = "lower", lab = TRUE,
                       ggtheme = theme_gray(),
                       colors = c("#3EDBF0", "white", "#FF75A0")) +
  theme(
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 10, face = "bold"),
    axis.text = element_text(color = "black", size = 10),
    panel.background = element_rect(fill = "#F1F1F1", colour = NA),
    panel.grid.minor = element_line(colour = "white", size = 0.2),
    panel.grid.major = element_line(colour = "white", size = 0.2)
  ) +
  labs(title = "Wastewater")
cor_plot

cor_data <- cor_plot$data
p.labs_final <- inner_join(cor_data, p.df_long, by = c("Var1", "Var2"))
p.labs_final 

cor.plot.labs4 <- cor_plot +
  geom_text(data = p.labs_final,
            aes(x = Var1, y = Var2, label = lab),
            nudge_y = 0.25,
            size = 5)

cor.plot.labs4


ggsave(file="cor.plot.labsa4.jpeg", cor.plot.labs4,
       width= 80, height = 60, units = "mm", dpi=600)




Fig6new <- plot_grid(
  cor.plot.labsa1,   # first plot (a)
  cor.plot.labs4,    # second plot (b)
  ncol = 1,          # vertical layout
  labels = c("a", "b"),
  label_size = 12,
  rel_heights = c(1, 0.6)
)

ggsave(file="Fig6new.jpeg", Fig6new, width= 140, height = 190, units = "mm", dpi=900)

#### sfmD ecoli correlation p<0.001####
library(ggcorrplot)
library(reshape2)
library(ggcorrplot)
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(reshape2)
library(tidyverse)
library(vegan)
library(grid)
library(gridExtra)
library(readr)
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)

library(readr)
sfmdcorr <- read_csv("sfmdcorr.csv", col_types = cols(sfmD = col_number(), 
                                                      EC = col_number()))
head(sfmdcorr)

attach(sfmdcorr)
p.mat <- cor_pmat(sfmdcorr) #correlation matrix with p-values
head(p.mat)
p.mat
m <- cor(sfmdcorr) 
m

corr <- round(cor(sfmdcorr), 2)
corr

p.df <- as.data.frame(ggcorrplot::cor_pmat(sfmdcorr))
p.df

labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")
}

p.labs = p.df %>%
  mutate_all(labs.function)
p.labs
p.labs$Var1 = as.factor(rownames(p.labs))

p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")
p.labs
cor_plot = ggcorrplot(corr, hc.order = F, type = "lower",
                      lab = T, ggtheme = ggplot2::theme_gray, colors = c("#F4EEEE", "#FFE6E6", "#E1AFD1")) +
  theme(legend.text = element_text(color = "black", family = "Arial", size = 10), #detail site label
        legend.title = element_text(color = "black", family = "Arial", size = 10, face = "bold"), #site detail label
        axis.text = element_text(color = "black", family = "Arial", size = 2),
        panel.background = element_rect(fill = "grey95", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.2))

p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                  paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot.labs = cor_plot +
  geom_text(aes(x = p.labs$Var1,
                y = p.labs$Var2),
            label = p.labs$lab,
            nudge_y = 0.25,
            size = 5)

cor.plot.labs

ggsave(file="cor.plot.labs.jpeg", cor.plot.labs, width= 250, height = 250, units = "mm", dpi=600)



#####Model tc ONLY####
library(stats)
library(car)
library(MASS)


library(readr)
corr_cpx <- read_csv("corr.csv", col_types = cols(pbi = col_number(),
                                                  tc = col_number(), ec = col_number()))


head(corr_cpx)
LM <- corr_cpx[, 6:7]
head(LM)

model <- lm(pBI143 ~ TC, data = LM)
summary(model)
r_squared <- summary(model)$r.squared
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope <- coef(model)[2]
cat("pBI143 = ", round(intercept, 4), " + (", round(slope, 4), " * TC)\n")
cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.765
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

####tc EC temp #VIF > 5 # 84.4%#### #OC = 86.7 %
head(corr_cpx)
LM2 <- corr_cpx[, 5:8]
head(LM2)
model <- lm(pBI143 ~ TC + EC + Temp, data = LM2)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic

f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")


# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = LM2)
vif_values <- vif(model)
print(vif_values)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

# remove EC
# TC        EC      Temp
# 11.034418 10.772935  1.082981
# new oc 10.591900 10.309234  1.093569

####ALL tc+temp ONLY VIF < 5 #84.4%#### oc 86.7 %

model_tc_temp <- lm(pBI143 ~ TC + Temp, data = LM2)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

# pBI143 =  2.9864  + ( 0.9609  * TC) + ( -0.0638  * Temperature)
# new pBI143 =  3.0349  + ( 0.9876  * TC) + ( -0.061  * Temperature)


#### KM KF UP DW TC + EC TC EC####
head(corr_cpx)

df_Upstream <- corr_cpx %>% filter(Site == "Upstream")
print(df_Upstream)
df_Downstream <- corr_cpx %>% filter(Site == "Downstream")
print(df_Downstream)
df_Kamanashi <- corr_cpx %>% filter(Site == "Kamanashi")
print(df_Kamanashi)
df_Kofu <- corr_cpx %>% filter(Site == "Kofu")
print(df_Kofu)
df_Kofu_Kamanashi <- corr_cpx %>% filter(Site %in% c("Kofu", "Kamanashi"))
print(df_Kofu_Kamanashi)
df_Upstream_Downstream <- corr_cpx %>% filter(Site %in% c("Upstream", "Downstream"))
print(df_Upstream_Downstream)

#####km+KF TC + EC #VIF>5 #72.1%#### oc 75.3 %
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Kofu_Kamanashi)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Kofu_Kamanashi)
vif_values <- vif(model)
print(vif_values)

cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#pBI143 =  3.532  + ( 0.5559  * TC) + ( 0.4039  * EC) + ( -0.0742  * Temperature)
# oc pBI143 =  3.532  + ( 0.5559  * TC) + ( 0.4039  * EC) + ( -0.0742  * Temperature)

#####km+KF TC + temp #VIF<5 #72.1%####
model_tc_temp <- lm(pBI143 ~ TC + Temp, data = df_Kofu_Kamanashi)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####KM TC + EC #VIF > 5 #72.4%####
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Kamanashi)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")


# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Kamanashi)
vif_values <- vif(model)
print(vif_values)

cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")


#####KM TC + temp #VIF < 5 #72.4%####
model_tc_temp <- lm(pBI143 ~ TC + Temp, data = df_Kamanashi)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")




#####KF TC + EC #VIF > 5 #77.4%####
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Kofu)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Kofu)
vif_values <- vif(model)
print(vif_values)

cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####KF TC + TEMP #VIF < 5 #77.4%####
model_tc_temp <- lm(pBI143 ~ TC + Temp, data = df_Kofu)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####UP TC + EC #VIF < 5 #38.9%####
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Upstream)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Upstream)
vif_values <- vif(model)
print(vif_values)

cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####UP TC + TEMP #VIF < 5 #38.9%####
model_tc_temp <- lm(pBI143 ~ TC + Temp, data = df_Upstream)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####DW TC + EC #VIF < 5 #56.1%####
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Downstream)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Downstream)
vif_values <- vif(model)
print(vif_values)

cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####DW TC + TEMP #VIF < 5 #56.1%####
model_tc_temp <- lm(pBI143 ~ TC + Temp, data = df_Downstream)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####UPDW TC + EC #VIF < 5 #47.8%####
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Upstream_Downstream)
summary(model)
r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model)[1]
slope_TC <- coef(model)[2]
slope_EC <- coef(model)[3]
slope_Temp <- coef(model)[4]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_EC, 4), " * EC) + (", round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

# Check Variance Inflation Factor (VIF) for Model 2 to confirm multicollinearity
model <- lm(pBI143 ~ TC + EC + Temp, data = df_Upstream_Downstream)
vif_values <- vif(model)
print(vif_values)

cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#####UPDW TC + TEMP #VIF < 5 #47.8%####
model_tc_temp <- lm(pBI143 ~ TC + Temp, data = df_Upstream_Downstream)
summary(model_tc_temp)

r_squared <- summary(model)$r.squared
f_statistic <- summary(model)$fstatistic
f_value <- f_statistic[1]  # Extract F-value
p_value <- pf(f_value, f_statistic[2], f_statistic[3], lower.tail = FALSE)  # Compute p-value
p_value
aic_value <- AIC(model)
dw_test <- durbinWatsonTest(model)
dw_test
condition_number <- kappa(model.matrix(model), exact = TRUE)

intercept <- coef(model_tc_temp)[1]
slope_TC <- coef(model_tc_temp)[2]
slope_Temp <- coef(model_tc_temp)[3]

cat("pBI143 = ", round(intercept, 4), " + (", round(slope_TC, 4), " * TC) + (",
    round(slope_Temp, 4), " * Temperature)\n")

cat("Model Performance:\n")
cat("R² =", round(r_squared, 3), "→ The model explains", round(r_squared * 100, 1), "% of the variance in pBI143.\n") #R² = 0.844
cat("F-statistic =", round(f_value, 3), ", p =", round(p_value, 5), "→ Overall model significance.\n")
cat("AIC =", round(aic_value, 2), "→ Lower AIC indicates better model fit.\n")
cat("Condition Number =", round(condition_number, 1), "→ Multicollinearity concern if > 30.\n")

vif_values_new <- vif(model_tc_temp)
print(vif_values_new)

# Interpretation
cat("\nInterpretation of VIF:\n")
cat(" - VIF < 5: No multicollinearity concern.\n")
cat(" - VIF 5-10: Moderate multicollinearity.\n")
cat(" - VIF > 10: High multicollinearity (Consider removing variables).\n")

#### mean sd calculation####

library(readr)
library(dplyr)
library(tidyr)
library(writexl)

genedata <- read_csv("genedata.csv", col_types = cols(
  pBI = col_number(), crass = col_number(), gyrB = col_number(), HE183 = col_number(),
  intl = col_number(), sul = col_number(), invA = col_number(), sfmD = col_number(),
  HAdVs = col_number(), NoVGII = col_number(), NoVGI = col_number(),
  TC = col_number(), EC = col_number()
))
head(genedata)

markers <- c("pBI", "crass", "gyrB", "HE183", "intl", "sul", "invA", "sfmD", "HAdVs", "NoVGII", "NoVGI")
days_of_interest <- c(0, 1, 3, 5, 7, 10, 14)
sites <- c("wastewater", "downstream", "upstream")
temps <- c("4C", "15C", "25C")

wastewater <- genedata[1:76, ]
downstream <- genedata %>% filter(Type == "Downstream")
upstream <- genedata %>% filter(Type == "Upstream")

site_data_list <- list(
  wastewater = wastewater,
  downstream = downstream,
  upstream = upstream
)

output_df <- expand.grid(marker = markers, Days = days_of_interest) %>%
  arrange(marker, Days)

for (temp in temps) {
  for (site in names(site_data_list)) {
    data_site <- site_data_list[[site]]

    day0_data <- data_site %>%
      filter(Days == 0) %>%
      pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
      group_by(marker, Days) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = n(),
        sdm = sd / sqrt(n)
      ) %>%
      mutate(mean_sd = sprintf("%.2f ± %.2f", mean, sdm)) %>%
      ungroup()

    other_days_data <- data_site %>%
      filter(Days != 0, Temp == temp, Days %in% days_of_interest) %>%
      pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
      group_by(marker, Days) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = n(),
        sdm = sd / sqrt(n)
      ) %>%
      mutate(mean_sd = sprintf("%.1f ± %.1f", mean, sdm)) %>%
      ungroup()

    subdata <- bind_rows(day0_data, other_days_data)

    colname <- paste(temp, site, sep = "_")
    output_df[[colname]] <- NA

    for (i in 1:nrow(subdata)) {
      idx <- which(output_df$marker == subdata$marker[i] & output_df$Days == subdata$Days[i])
      output_df[[colname]][idx] <- subdata$mean_sd[i]
    }
  }
}

write_xlsx(list("Marker Data Summary" = output_df), "genedata_summary_sites_day0.xlsx")


#### mean sd significant test####
# Load libraries
library(tidyverse)
library(rstatix)
library(ggpubr)

# Check column names (just to be sure)
print(colnames(genedata))

# 1. Filter Day 0 rows
day0_data <- genedata %>%
  filter(Days == 0)

# 2. Choose only relevant columns — use exact names from colnames() above
# Include only DNA/RNA marker columns you care about
marker_columns <- c("pBI", "crass", "gyrB", "HE183", "intl", "sul", "HAdVs")

# 3. Add an ID for pairing in statistical tests
day0_long <- day0_data %>%
  mutate(SampleID = row_number()) %>%
  pivot_longer(cols = all_of(marker_columns), names_to = "Marker", values_to = "Concentration")

# 4. Summary: mean per marker
summary_stats <- day0_long %>%
  group_by(Marker) %>%
  summarise(mean_conc = mean(Concentration, na.rm = TRUE)) %>%
  arrange(desc(mean_conc))

print(summary_stats)

# 5. Plot the marker concentrations
ggplot(day0_long, aes(x = Marker, y = Concentration)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  labs(title = "Marker Concentrations on Day 0", y = "log10 Concentration", x = "Marker")

# 6. Friedman test (non-parametric repeated measures)
friedman_result <- day0_long %>%
  friedman_test(Concentration ~ Marker | SampleID)

print(friedman_result)

# 7. Post-hoc test if significant
if (friedman_result$p < 0.05) {
  pairwise_result <- day0_long %>%
    pairwise_wilcox_test(Concentration ~ Marker, paired = TRUE, p.adjust.method = "BH")
  print(pairwise_result)
}

# Load libraries
library(tidyverse)
library(rstatix)
library(ggpubr)

# Filter only Day 0
day0_data <- genedata %>%
  filter(Days == 0)

# Define markers to include
marker_columns <- c("pBI", "crass", "gyrB", "HE183", "intl", "sul", "HAdVs")

# Reshape to long format with site info
day0_long <- day0_data %>%
  mutate(SampleID = row_number()) %>%
  pivot_longer(cols = all_of(marker_columns), names_to = "Marker", values_to = "Concentration")

# 1. Reshape data and keep Type before grouping
day0_long <- genedata %>%
  filter(Days == 0) %>%
  mutate(SampleID = row_number()) %>%
  pivot_longer(cols = c(pBI, crass, gyrB, HE183, intl, sul, HAdVs),
               names_to = "Marker",
               values_to = "Concentration")

# 2. Run Friedman test per site
friedman_by_site <- day0_long %>%
  group_by(Type) %>%
  group_modify(~ {
    friedman_result <- .x %>%
      friedman_test(Concentration ~ Marker | SampleID)

    friedman_result <- friedman_result %>%
      mutate(Site = .y$Type) %>%
      dplyr::select(Site, everything())

    return(friedman_result)
  })

print(friedman_by_site)

pairwise_by_site <- day0_long %>%
  group_by(Type) %>%
  group_modify(~ {
    test <- .x %>%
      pairwise_wilcox_test(Concentration ~ Marker, paired = TRUE, p.adjust.method = "BH") %>%
      mutate(Site = .y$Type) %>%
      dplyr::select(Site, everything())
    return(test)
  })

signif_pairwise_by_site <- pairwise_by_site %>%
  filter(p.adj.signif != "ns")  # "ns" = not significant

print(signif_pairwise_by_site)

print(pairwise_by_site)

min(crass)
