library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
 
################################################
################ Density Plots #################
################################################
 
#Chane the path and modality here -->
setwd("/data/Stereogene/Dmod_v_Dmod/")
modality = "Dmod_Dmod"
 
dirs <- list.dirs(".", recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("^./results_Avg_", dirs)]

all_data <- list()
 
for (d in dirs) {
  bkg_file <- list.files(d, pattern = "\\.bkg$", full.names = TRUE)
  fg_file  <- list.files(d, pattern = "\\.fg$", full.names = TRUE)
   
   if (length(bkg_file) == 1 && length(fg_file) == 1) {
     bkg <- read.table(bkg_file, sep = "\t")[,1]   # 1
     fg  <- read.table(fg_file, sep="\t")[,4]    # 4
     
     
     base_name <- sub("^results_Avg_", "", basename(d))
     base_name <- sub("Avg.","",base_name)
     label <- base_name
     
     all_data[[length(all_data)+1]] <- data.frame(
       value = c(bkg, fg),
       type  = rep(c("bkg", "fg"), c(length(bkg), length(fg))),
       sample = label
     )
   }
 }
 
 plot_data <- bind_rows(all_data)
 
 plot_data$sample_modified <- gsub(
   "_WT_bin1000_sigma1500|_bin1000_sigma1500|\\.neglog10PValue_sign\\.bedGraph|\\_WT\\_KO\\.bedGraph|",
   "",
   plot_data$sample
 )
 
 plot_data$sample_modified <- plot_data$sample_modified <- gsub(
   "\\_Delta\\_RNA\\.bedGraph",
   "_RNA",
   plot_data$sample_modified
 )
 
 
#### --------------------------------- ###
#### combined density plot violin plot ###
#### --------------------------------- ###
plot_data$sample[plot_data$type == "bkg"] <- "Background"
min_value <-  min(plot_data$value)
max_value <- max(plot_data$value)
bkg_data <- plot_data[ which( plot_data$type == "bkg"),]
bkg_data <- na.omit(bkg_data)
fg_data <- plot_data[ which( plot_data$type == "fg"),]
  
all_samples <- unique(c(fg_data$sample_modified, bkg_data$sample_modified))
my_colors <- setNames(brewer.pal(length(all_samples), "Dark2")[1:length(all_samples)], all_samples)
 
p_fg <- ggplot(fg_data, aes(x = sample_modified, y = value, fill = sample_modified)) +
   geom_violin(alpha = 0.3, trim = FALSE) +
   geom_boxplot(width = 0.1, alpha = 0.3, outlier.size = 0.5) +
   geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.3) +
   scale_y_continuous(limits = c(min_value, max_value)) + 
   scale_fill_manual(values = my_colors) + 
   theme_bw() + 
   theme(
     legend.position = "none",
     axis.text.x = element_text(angle = 45, hjust = 1)
   )
 
p_bkg <- ggplot(bkg_data, aes(x = "Background", y = value, fill = sample_modified, color = sample_modified)) +
   geom_violin(position = "identity", alpha = 0.3, trim = FALSE, size = 0.4) +
   geom_boxplot(width = 0.1, position = position_identity(),
                alpha = 0.3, outlier.size = 0.5, color = "black") +
   geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.3) +
   scale_y_continuous(limits = c(min_value, max_value)) +
   scale_fill_manual(values = my_colors) +
   scale_color_manual(values = my_colors) +
   theme_bw() +
   theme(
     axis.text.x = element_text(angle = 45, hjust = 1),
     legend.position = "none"
   ) +
   labs(x = "", y = "Value")
 
combined <- plot_grid(p_bkg, p_fg, rel_widths = c(1, 2), ncol = 2, align = "h")
combined
 
ggsave(paste0(modality,"_density.pdf"),
        plot = combined,
        device = cairo_pdf,
        width = 9.51, height = 5.3,   # in inches
       units = "in")
