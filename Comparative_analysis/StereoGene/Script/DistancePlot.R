### Distance Plots script
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

####Instructions:
## Change the path and modality here
setwd("/data/StereoGene/Dmod_v_Dmod")
modality = "Dmod_Dmod"


dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("results_Avg_", dirs)]
all_data <- list()

for (d in dirs) {
  dist_file <- list.files(d, pattern = "\\.dist$", full.names = TRUE)
  
  if (length(dist_file) == 1) {
    df <- read.table(dist_file, header = TRUE, sep = "\t")
    base_name <- sub("^results_Avg_", "", basename(d))
    base_name <- sub("Avg.","",base_name)
    label <- base_name
    df_long <- df %>%
      select(dist, Bkg, Fg) %>%
      mutate(dist = dist/1000, sample = label) %>%
      pivot_longer(cols = c("Bkg", "Fg"), names_to = "type", values_to = "value")
    all_data[[length(all_data)+1]] <- df_long
  }
}

plot_data <- bind_rows(all_data)

###############################################
##   Signal Distance Plot (Facet-Wrap)  ###
###############################################

plot_data$sample_modified <- plot_data$sample

plot_data$sample_modified <- gsub(
  "_WT_bin1000_sigma1500|_bin1000_sigma1500|\\.neglog10PValue_sign\\.bedGraph|_WT_KO\\.bedGraph",
  "",
  plot_data$sample_modified
)

plot_data$sample_modified <- gsub(
  "_Delta_RNA\\.bedGraph",
  "_RNA",
  plot_data$sample_modified
)

plot_data <- plot_data %>% filter(dist >= -10 & dist <= 10)


plot_data <- plot_data %>%
  mutate(value_log2 = log2(abs(value) + 1) * sign(value))

min_value <- min(plot_data$value_log2, na.rm = TRUE)
max_value <- max(plot_data$value_log2, na.rm = TRUE)

unique_samples <- unique(plot_data$sample_modified)

if (length(unique_samples) > 9) {
  colorMap <- setNames(
    colorRampPalette(brewer.pal(9, "Dark2"))(length(unique_samples)),
    unique_samples
  )
} else {
  colorMap <- setNames(
    brewer.pal(n = max(3, length(unique_samples)), "Dark2")[seq_along(unique_samples)],
    unique_samples
  )
}

Distance <- ggplot(plot_data, aes(
  x = dist,
  y = value_log2,
  color = sample_modified,
  linetype = type
)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.3) +
  facet_wrap(~sample_modified, scales = "free_y") +
  labs(
    x = "Distance (kb)",
    y = "log2(Density)",
    linetype = "Signal"
  ) +
  scale_color_manual(values = colorMap, guide = "none") +  # <-- remove color legend
  scale_linetype_manual(values = c("Fg" = "solid", "Bkg" = "dashed")) +
  scale_y_continuous(limits = c(min_value, max_value)) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

Distance
ggsave( paste0(modality, "_distance.pdf"),
        plot = Distance,
        device = cairo_pdf,
        width = 9.51, height = 5.3,   # in inches
        units = "in")


