# pBI143 and temperature 4, 15, 25 c#

##### Load necessary libraries####
getwd()
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
library(segmented)
library(openxlsx)
library(readr)

#### k and T90 value but unknown breakpoint day####
lnraw <- read_csv("lnraw.csv")
head(lnraw)

df <- lnraw %>%
  pivot_longer(
    cols = pBI143:EC,   # pivot all gene columns (4th to 16th) into long format
    names_to = "gene",
    values_to = "ln_Ct_C0"
  )
head(df)

results <- list()
df <- df %>%
  mutate(site = case_when(
    site %in% c("Kamanashi", "Kofu") ~ "wastewater",
    TRUE ~ site
  ))

results <- list()

df %>%
  group_by(gene, temperature, site) %>%
  group_split() %>%
  walk(function(subdf) {
    info <- subdf[1, c("gene", "temperature", "site")]  # Metadata

    # Fit basic linear model
    mod <- lm(ln_Ct_C0 ~ time, data = subdf)

    # Try segmented model
    has_break <- FALSE
    try({
      seg_mod <- segmented(mod, seg.Z = ~time, psi = list(time = 5))
      pval <- davies.test(mod, seg.Z = ~ time)$p.value
      if (pval < 0.05) has_break <- TRUE
    }, silent = TRUE)

    if (has_break) {
      s <- slope(seg_mod)$time[, "Est."]
      ks <- -s
      t90s <- 2.303 / ks
      results <<- append(results, list(data.frame(
        gene = info$gene,
        temperature = info$temperature,
        site = info$site,
        k_type = c("k1", "k2"),
        k = ks,
        T90 = t90s
      )))
    } else {
      k <- -coef(mod)["time"]
      T90 <- 2.303 / k
      results <<- append(results, list(data.frame(
        gene = info$gene,
        temperature = info$temperature,
        site = info$site,
        k_type = "k",
        k = k,
        T90 = T90
      )))
    }
  })

final_results <- bind_rows(results)
head(final_results)
write.xlsx(final_results, "decay_rates_output2.xlsx")


#####Find 95% CI but k type only not for k1 and k2 so can not use this####

library(tidyverse)
library(segmented)
library(openxlsx)

# Load and clean data
lnraw <- read_csv("lnraw.csv")
head(lnraw)

df <- lnraw %>%
  pivot_longer(
    cols = pBI143:EC,   # pivot all gene columns (4th to 16th) into long format
    names_to = "gene",
    values_to = "ln_Ct_C0"
  )
head(df)

df <- df %>%
  mutate(
    site = case_when(site %in% c("Kamanashi", "Kofu") ~ "wastewater", TRUE ~ site),
    temperature = as.character(temperature),
    ln_Ct_C0 = as.numeric(ln_Ct_C0)
  )

# Prepare container
results <- list()

# Group and analyze
df %>%
  group_by(gene, temperature, site) %>%
  group_split() %>%
  walk(function(subdf) {
    info <- subdf[1, c("gene", "temperature", "site")]

    if (nrow(subdf) < 3 || all(is.na(subdf$ln_Ct_C0))) return(NULL)

    mod <- lm(ln_Ct_C0 ~ time, data = subdf)
    has_break <- FALSE
    seg_mod <- NULL

    # Attempt segmented model
    try({
      seg_mod <- segmented(mod, seg.Z = ~time, psi = list(time = 5))
      pval <- davies.test(mod, seg.Z = ~ time)$p.value
      if (!is.na(pval) && pval < 0.05) has_break <- TRUE
    }, silent = TRUE)

    if (has_break && !is.null(seg_mod)) {
      # Try to extract slope info safely
      slopes_df <- tryCatch(slope(seg_mod)$time, error = function(e) NULL)

      if (!is.null(slopes_df) &&
          all(c("Est.", "St.Err") %in% colnames(slopes_df)) &&
          nrow(slopes_df) == 2) {

        s <- slopes_df[, "Est."]
        se <- slopes_df[, "St.Err"]
        ks <- -s
        k_lower <- -s + 1.96 * se
        k_upper <- -s - 1.96 * se
        t90s <- 2.303 / ks
        t90_lower <- 2.303 / k_upper
        t90_upper <- 2.303 / k_lower

        results <<- append(results, list(data.frame(
          gene = info$gene,
          temperature = info$temperature,
          site = info$site,
          k_type = c("k1", "k2"),
          k = ks,
          k_lower_95 = k_lower,
          k_upper_95 = k_upper,
          T90 = t90s,
          T90_lower_95 = t90_lower,
          T90_upper_95 = t90_upper
        )))
        return()
      }
    }

    # Fallback: simple model
    k <- -coef(mod)["time"]
    se <- summary(mod)$coefficients["time", "Std. Error"]
    k_lower <- k - 1.96 * se
    k_upper <- k + 1.96 * se
    T90 <- 2.303 / k
    T90_lower <- 2.303 / k_upper
    T90_upper <- 2.303 / k_lower

    results <<- append(results, list(data.frame(
      gene = info$gene,
      temperature = info$temperature,
      site = info$site,
      k_type = "k",
      k = k,
      k_lower_95 = k_lower,
      k_upper_95 = k_upper,
      T90 = T90,
      T90_lower_95 = T90_lower,
      T90_upper_95 = T90_upper
    )))
  })

# Export results
final_results <- bind_rows(results)
head(final_results)
write.xlsx(final_results, "decay_rates_output_with_CI.xlsx")

######## decay k_type to find break-point days####
library(tidyverse)
library(segmented)

# Prepare a results list
k_results <- list()

# Group and iterate by gene, site, and temperature
df %>%
  group_by(gene, site, temperature) %>%
  group_split() %>%
  walk(function(subdf) {
    if (nrow(subdf) < 2 || all(is.na(subdf$ln_Ct_C0))) return(NULL)

    # Fit linear model
    base_mod <- lm(ln_Ct_C0 ~ time, data = subdf)

    # Default values
    k_type <- "k"
    k1 <- k2 <- NA
    bp <- NA

    # Try segmented
    try({
      seg_mod <- segmented(base_mod, seg.Z = ~ time, psi = list(time = 5))
      pval <- davies.test(base_mod, seg.Z = ~ time)$p.value
      if (!is.na(pval) && pval < 0.05) {
        k_type <- "k1_k2"
        bp <- as.numeric(seg_mod$psi[2])
        slopes <- slope(seg_mod)$time[, "Est."]
        k1 <- slopes[1]
        k2 <- slopes[2]
      } else {
        # Simple decay slope
        k1 <- coef(base_mod)["time"]
      }
    }, silent = TRUE)

    # Save results
    k_results <<- append(k_results, list(tibble(
      gene = subdf$gene[1],
      site = subdf$site[1],
      temperature = subdf$temperature[1],
      k_type = k_type,
      slope1 = k1,
      slope2 = k2,
      breakpoint = bp
    )))
  })

# Combine results into a data frame
k_summary <- bind_rows(k_results)
head(k_summary)

# View the summary
print(k_summary)
write.xlsx(k_summary, "k_summary.xlsx")

##### Summarize mean ± SD for k NOT WORK because no SD ####
library(dplyr)
library(tidyr)
library(writexl)

# Example: Load your data (final_results)
# final_results <- read_csv("decay_rates_output_with_CI.csv")

# Prepare mean ± SD summary for k and T90
summary_stats <- final_results %>%
  group_by(gene, temperature, site) %>%
  summarise(
    mean_k = mean(k, na.rm = TRUE),
    sd_k = sd(k, na.rm = TRUE),
    mean_T90 = mean(T90, na.rm = TRUE),
    sd_T90 = sd(T90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    k_mean_sd = sprintf("%.4f ± %.4f", mean_k, sd_k),
    T90_mean_sd = sprintf("%.2f ± %.2f", mean_T90, sd_T90)
  )

# Merge k and T90 mean ± SD into one long format table
summary_long <- summary_stats %>%
  pivot_longer(
    cols = c(k_mean_sd, T90_mean_sd),
    names_to = "metric",
    values_to = "mean_sd"
  )

# Adjust temperature to numeric (remove "C")
summary_long$temperature <- as.numeric(gsub("C", "", summary_long$temperature))

# Reshape to wide format: site + metric
final_table <- summary_long %>%
  unite(site_metric, site, metric, sep = "_") %>%
  pivot_wider(
    names_from = site_metric,
    values_from = mean_sd
  ) %>%
  arrange(gene, temperature)

# Export to Excel
write_xlsx(list("Summary" = final_table), "marker_site_summary_k_T90.xlsx")

cat("Excel file 'marker_site_summary_k_T90.xlsx' generated with mean ± SD for k and T90!\n")

#### PCA ####
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
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(vegan)
library(grid)
library(gridExtra)
library(readr)
library(FactoMineR)
library(Factoshiny)
library(missMDA)
library(FactoInvestigate)
library(car)
library(factoextra)

# Load libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpubr)  # For ggarrange()

# Read data
lnraw <- read_csv("lnraw.csv", col_types = cols(
  temperature = col_character(),
  site = col_character(),
  time = col_character(),
  pBI = col_number(), crass = col_number(),
  gyrB = col_number(), HE183 = col_number(),
  intl = col_number(), sul = col_number(),
  invA = col_number(), sfmD = col_number(),
  HAdVs = col_number(), NoVGII = col_number(),
  NoVGI = col_number(), TC = col_number(),
  EC = col_number()
))

# Convert to factors
lnraw$time <- factor(lnraw$time, levels = unique(lnraw$time))
lnraw$temperature <- factor(lnraw$temperature)

# Select genes of interest
df_selected <- lnraw[, c("time", "temperature", "pBI", "crass", "gyrB", "HE183", "intl", "sul")]

# PCA
pca_res <- PCA(df_selected[, c("pBI", "crass", "gyrB", "HE183", "intl", "sul")], graph = FALSE)

# PCA plot for individuals
ind_plot <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  col.ind = df_selected$time,          # Color: time (factor)
  shape.ind = df_selected$temperature, # Shape: temperature
  pointsize = 2,
  addEllipses = TRUE, ellipse.type = "confidence",
  palette = "Dark2",
  legend.title = list(fill = "Day", shape = "Temperature")
) +
  theme(
    legend.text = element_text(size = 7, color = "black"),
    legend.title = element_text(size = 7, face = "bold", color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, face = "bold", color = "black")
  )

# PCA plot for variables
var_plot <- fviz_pca_var(
  pca_res,
  col.var = "contrib",                 # Gradient for contributions
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
) +
  theme(
    legend.text = element_text(size = 7, color = "black"),
    legend.title = element_text(size = 7, face = "bold", color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, face = "bold", color = "black")
  )

# Combine plots
final_plot <- ggarrange(ind_plot, var_plot, ncol = 2, widths = c(1, 1))

# Save and show
ggsave("PCAplot_combined.jpeg", plot = final_plot, width = 10, height = 5, dpi = 300)
print(final_plot)

# Load libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpubr)

# Read your data
lnraw <- read_csv("lnraw.csv", col_types = cols(
  temperature = col_character(),
  site = col_character(),
  time = col_character(),
  pBI = col_number(), crass = col_number(),
  gyrB = col_number(), HE183 = col_number(),
  intl = col_number(), sul = col_number(),
  invA = col_number(), sfmD = col_number(),
  HAdVs = col_number(), NoVGII = col_number(),
  NoVGI = col_number(), TC = col_number(),
  EC = col_number()
))

# Convert to factors
lnraw$temperature <- factor(lnraw$temperature)

# Select genes of interest
df_selected <- lnraw[, c("temperature", "pBI", "crass", "gyrB", "HE183", "intl", "sul")]

# Perform PCA
pca_res <- PCA(df_selected[, c("pBI", "crass", "gyrB", "HE183", "intl", "sul")], graph = FALSE)

# PCA plot for individuals (colored by temperature)
ind_plot <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  col.ind = df_selected$temperature,  # Color by temperature
  pointsize = 2,
  palette = "Dark2",
  addEllipses = TRUE, ellipse.type = "confidence",
  legend.title = list(fill = "Temperature")
) +
  theme(
    legend.text = element_text(size = 7, color = "black"),
    legend.title = element_text(size = 7, face = "bold", color = "black"),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, face = "bold", color = "black")
  )

# PCA plot for variables (genes)
var_plot <- fviz_pca_var(
  pca_res,
  col.var = "black",  # One color for all variables
  repel = TRUE
) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, face = "bold", color = "black")
  )

# Combine plots
final_plot <- ggarrange(ind_plot, var_plot, ncol = 2, widths = c(1, 1))

# Save and show
ggsave("PCAplot_genes_temperature.jpeg", plot = final_plot, width = 10, height = 5, dpi = 300)
print(final_plot)

####PCA temperature####
library(tidyverse)
library(FactoMineR)
library(factoextra)

lnraw <- read_csv("lnraw.csv", col_types = cols(
  temperature = col_character(),
  site = col_character(),
  time = col_character(),
  pBI143 = col_number(), CrAssphage = col_number(),
  gyrB = col_number(), HF183 = col_number(),
  intl1 = col_number(), sul1 = col_number(),
  invA = col_number(), sfmD = col_number(),
  HAdVs = col_number(), NoVGII = col_number(),
  NoVGI = col_number(), TC = col_number(),
  EC = col_number()
))

lnraw$temperature <- factor(lnraw$temperature)
####Biplot
df_selected <- lnraw[, c("temperature", "pBI143", "CrAssphage", "gyrB", "HF183", "intl1", "sul1")]
pca_res <- PCA(df_selected[, c("pBI143", "CrAssphage", "gyrB", "HF183", "intl1", "sul1")], graph = FALSE)
fviz_pca_var(pca_res, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
             )

p2 <- fviz_pca_biplot(pca_res,
                     geom.ind = "point",
                     fill.ind = df_selected$temperature,
                     col.ind = "black",
                     pointshape = 21, pointsize = 2,
                     palette = "startrek",
                     addEllipses = TRUE, ellipse.type = "confidence",
                     col.var = "contrib",
                     repel = TRUE,
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     legend.title = list(fill = "temperature", color = "Contrib")
                     )
p2

####Biplot poster
df_selected <- lnraw[, c("temperature", "pBI143", "CrAssphage", "gyrB", "HF183")]
pca_res <- PCA(df_selected[, c("pBI143", "CrAssphage", "gyrB", "HF183")], graph = FALSE)
fviz_pca_var(pca_res, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

p2 <- fviz_pca_biplot(pca_res,
                      geom.ind = "point",
                      fill.ind = df_selected$temperature,
                      col.ind = "black",
                      pointshape = 21, pointsize = 2,
                      palette = "startrek",
                      addEllipses = TRUE, ellipse.type = "confidence",
                      col.var = "contrib",
                      repel = TRUE,
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      legend.title = list(fill = "temperature", color = "Contrib")
)
p2

px2 <- p2 +
  xlab("PC1 (83.8%)") + labs(fill = "Temperature") + ylab("PC2 (10%)") + labs(title = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
px2
ggsave(file="PCAplottemp2.jpeg", px2, width= 190, height = 160, units = "mm", dpi=600)


####PCA temp####
lnraw <- read_csv("lnraw.csv", col_types = cols(
  temperature = col_character(),
  site = col_character(),
  time = col_character(),
  pBI143 = col_number(), crAssphage = col_number(),
  gyrB = col_number(), HF183 = col_number(),
  intl1 = col_number(), sul1 = col_number(),
  invA = col_number(), sfmD = col_number(),
  HAdVs = col_number(), NoVGII = col_number(),
  NoVGI = col_number(), TC = col_number(),
  EC = col_number()
))

head(lnraw)

lnraw$temperature <- factor(lnraw$temperature)
lnraw$temperature <- factor(lnraw$temperature, levels = c("4C", "15C", "25C"))
df_selected <- lnraw[, c("temperature", "pBI143", "crAssphage", "gyrB", "HF183")]
pca_res <- PCA(df_selected[, c("pBI143", "crAssphage", "gyrB", "HF183")], graph = FALSE)

fviz_pca_var(pca_res, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

p1 <- fviz_pca_biplot(pca_res,
                      geom.ind = "point",
                      fill.ind = df_selected$temperature,
                      col.ind = "black",
                      pointshape = 21, pointsize = 3,
                      palette = "startrek",
                      addEllipses = TRUE, ellipse.type = "confidence",
                      col.var = "contrib",
                      repel = TRUE,
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      labelsize = 6,  # <-- Increase gene label size
                      legend.title = list(fill = "Temperature", color = "Contrib")
)
p1

px1 <- p1 +
  scale_shape_manual(values = c(21, 22, 24)) +  # assign shapes to 4C, 15C, 25C
  xlab("PC1 (88.2%)") + labs(fill = "Temperature", shape = "Temperature") +
  ylab("PC2 (7.2%)") +
  labs(title = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 20, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 20, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 20), #detail site label
        legend.title = element_text(color = "black", size = 20, face = "bold"), #site detail label
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20),
        strip.text.y = element_text(color = "black", size = 20, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
px1
ggsave(file="PCAplottemperature1.jpeg", px1, width= 220, height = 200, units = "mm", dpi=600)


library(ggrepel)

px1 <- p1 +
  geom_text_repel(aes(label = df_selected$time), size = 5, color = "black", max.overlaps = 50) +  # เพิ่ม label เป็นวัน
  scale_shape_manual(values = c(21, 22, 24)) +
  xlab("PC1 (88.2%)") + labs(fill = "Temperature", shape = "Temperature") +
  ylab("PC2 (7.2%)") +
  labs(title = NULL) +
  theme(
    axis.title.x = element_text(color = "black", size = 20, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20, face = "bold"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),
    strip.text.y = element_text(color = "black", size = 20, face = "bold"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
px1

#### PCA time ####

lnraw <- read_csv("lnraw.csv", col_types = cols(
  temperature = col_character(),
  site = col_character(),
  time = col_character(),
  pBI143 = col_number(), CrAssphage = col_number(),
  gyrB = col_number(), HF183 = col_number(),
  intl1 = col_number(), sul1 = col_number(),
  invA = col_number(), sfmD = col_number(),
  HAdVs = col_number(), NoVGII = col_number(),
  NoVGI = col_number(), TC = col_number(),
  EC = col_number()
))

lnraw$time <- factor(lnraw$time)
df_selected <- lnraw[, c("time", "pBI143", "crAssphage", "gyrB", "HE183", "intl1", "sul1")]
pca_res <- PCA(df_selected[, c("pBI143", "crAssphage", "gyrB", "HE183", "intl1", "sul1")], graph = FALSE)

fviz_pca_var(pca_res, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

p <- fviz_pca_biplot(pca_res,
                     geom.ind = "point",
                     fill.ind = df_selected$time,
                     col.ind = "black",
                     pointshape = 21, pointsize = 2,
                     palette = "startrek",
                     addEllipses = TRUE, ellipse.type = "confidence",
                     col.var = "contrib",
                     repel = TRUE,
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     legend.title = list(fill = "time", color = "Contrib")
)
p

px <- p +
  xlab("PC1 (83.8%)") + labs(fill = "time") + ylab("PC2 (10%)") + labs(title = NULL) +
  theme(axis.title.x = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนx
        axis.title.y = element_text(color = "black", size = 7, face = "bold"), #ชื่อแกนy
        legend.text = element_text(color = "black", size = 7), #detail site label
        legend.title = element_text(color = "black", size = 7, face = "bold"), #site detail label
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(color = "black", size = 7),
        strip.text.y = element_text(color = "black", size = 5, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "white", size = 0.1),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
px
save_plot("PCAplottime.jpeg", px)


####NDMS####

# Load required libraries
library(vegan)
library(tidyverse)

# Step 1: Select gene columns and metadata
gene_data <- lnraw[, c("pBI143", "CrAssphage", "gyrB", "HF183", "intl1", "sul1")]
meta_time <- lnraw$time

# Step 2: Remove rows with NA
complete_cases <- complete.cases(gene_data)
gene_data_clean <- gene_data[complete_cases, ]
meta_time_clean <- meta_time[complete_cases]

# Step 3: Remove rows with all 0 (needed for Bray-Curtis)
nonzero_rows <- rowSums(gene_data_clean) > 0
gene_data_clean <- gene_data_clean[nonzero_rows, ]
meta_time_clean <- meta_time_clean[nonzero_rows]

# Step 4: Compute Bray-Curtis distance matrix
dist_matrix <- vegdist(gene_data_clean, method = "bray")

# Step 5: Perform NMDS
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100)

# Step 6: Prepare NMDS result for ggplot
nmds_df <- as.data.frame(nmds_result$points)
nmds_df$time <- factor(meta_time_clean)

# Step 7: Plot
ggplot(nmds_df, aes(MDS1, MDS2, color = time)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS Plot Grouped by Time", color = "Time") +
  theme(text = element_text(size = 12, face = "bold"))



####NDMs FINAL####
library(ggrepel)
library(vegan)
library(tidyverse)

head(lnraw)
# Step 1: Subset data
Wastewater <- lnraw[1:84, ]
Downstream <- lnraw %>% filter(site == "Downstream")
Upstream <- lnraw %>% filter(site == "Upstream")

# Step 2: Add new column to each subset
Wastewater$site_2 <- "Wastewater"
Downstream$site_2 <- "Downstream"
Upstream$site_2 <- "Upstream"

# Step 3: Combine into one dataframe
lnraw <- bind_rows(Wastewater, Downstream, Upstream)


# Step 1: Select gene and metadata
gene_data <- lnraw[, c("pBI143", "CrAssphage", "gyrB", "HF183", "intl1", "sul1")]
meta_temp <- lnraw$temperature
meta_site2 <- lnraw$site_2  # new column for label

# Step 2: Filter out missing values
complete_cases <- complete.cases(gene_data, meta_temp, meta_site2)
gene_data_clean <- gene_data[complete_cases, ]
meta_temp_clean <- meta_temp[complete_cases]
meta_site2_clean <- meta_site2[complete_cases]

# Step 3: Remove rows with all zero
nonzero_rows <- rowSums(gene_data_clean) > 0
gene_data_clean <- gene_data_clean[nonzero_rows, ]
meta_temp_clean <- meta_temp_clean[nonzero_rows]
meta_site2_clean <- meta_site2_clean[nonzero_rows]

# Step 4: Compute Bray-Curtis distance
dist_matrix <- vegdist(gene_data_clean, method = "bray")

# Step 5: Run NMDS
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100)

# Step 6: Prepare NMDS result for plotting
nmds_df <- as.data.frame(nmds_result$points)
nmds_df$temperature <- factor(meta_temp_clean, levels = c("4C", "15C", "25C"))
nmds_df$site_2 <- meta_site2_clean

# Step 7: Fit gene vectors
fit <- envfit(nmds_result, gene_data_clean, perm = 999)
arrow_df <- as.data.frame(scores(fit, "vectors"))
arrow_df$gene <- rownames(arrow_df)

# Step: Calculate centroids for each group
centroids <- nmds_df %>%
  group_by(site_2) %>%
  summarize(MDS1 = mean(MDS1), MDS2 = mean(MDS2))

# NMDS Plot with PERMANOVA annotation
NMDS <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, shape = temperature, color = site_2)) +
  geom_point(aes(fill = site_2, shape = temperature),
             color = "black", size = 4, stroke = 1) +
  stat_ellipse(aes(group = site_2, fill = site_2),
               geom = "polygon", alpha = 0.2, color = NA) +
  geom_text_repel(data = centroids,
                  aes(x = MDS1, y = MDS2, label = NA, color = site_2),
                  fontface = "bold", size = 5, inherit.aes = FALSE,
                  nudge_y = 0.2) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = arrow_df,
                  aes(x = NMDS1, y = NMDS2, label = NA),
                  size = 4, fontface = "bold", inherit.aes = FALSE) +
  scale_color_manual(values = c(
    "Wastewater" = "red",
    "Upstream" = "green",
    "Downstream" = "orange"
  )) +
  scale_fill_manual(values = c(
    "Wastewater" = "red",
    "Upstream" = "green",
    "Downstream" = "orange"
  )) +
  scale_shape_manual(values = c("4C" = 21, "15C" = 22, "25C" = 24)) +
  theme_minimal() +
  labs(
    color = "Water Source",
    fill = "Water Source",
    shape = "Temperature"
  ) +
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text.y = element_text(color = "black", size = 10, face = "bold"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "mm"),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.1),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )


NMDS
ggsave(file="NMDS.jpeg", NMDS, width= 200, height = 150, units = "mm", dpi=800)


##### PERMANOVA####
library(vegan)
library(dplyr)
library(tidyr)
library(knitr)

# Step 1: Define gene names
genes <- c("pBI143", "CrAssphage", "gyrB", "HF183", "intl1", "sul1")

# Step 2: Initialize empty list
results_list <- list()

# Step 3: Loop over genes and run PERMANOVA
for (gene in genes) {

  # Extract relevant data
  gene_vector <- lnraw[[gene]]
  meta_df <- data.frame(
    Gene = gene_vector,
    Temperature = lnraw$temperature,
    Matrix = lnraw$site_2
  )

  # Remove missing and zero rows
  complete_rows <- complete.cases(meta_df) & meta_df$Gene > 0
  meta_clean <- meta_df[complete_rows, , drop = FALSE]

  if (nrow(meta_clean) > 2) {
    # Create dissimilarity matrix
    dist_matrix <- vegdist(data.frame(value = meta_clean$Gene), method = "bray")

    # PERMANOVA
    adonis_result <- adonis2(dist_matrix ~ Temperature * Matrix, data = meta_clean, permutations = 999)

    # Format output
    df <- data.frame(
      Term = rownames(adonis_result),
      R2 = round(100 * adonis_result$R2, 1),
      p_value = ifelse(adonis_result$`Pr(>F)` < 0.0001, "<0.0001", signif(adonis_result$`Pr(>F)`, 3)),
      Gene = gene
    )

    # Filter only Temperature, Matrix, and Interaction rows
    df <- df %>% filter(Term %in% c("Temperature", "Matrix", "Temperature:Matrix"))

    results_list[[gene]] <- df
  }
}

# Step 4: Combine all results
combined_df <- bind_rows(results_list)

# Step 5: Reshape to match desired output
summary_table <- combined_df %>%
  pivot_wider(
    names_from = Gene,
    values_from = c(R2, p_value),
    names_glue = "{Gene} ({.value})"
  ) %>%
  rename(`Source of variation` = Term)

# Step 6: Display table
kable(summary_table, caption = "PERMANOVA: % Variation and p-values by Gene for Temperature, Matrix, and Interaction")

####test####
# Load required libraries
