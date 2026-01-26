library(ggplot2)

#set working directory to source file directory as all paths are relative to it
gwdir<-getwd() #getSrcDirectory(function(){})[1]
setwd(gwdir)
getwd()

#read in GV/GO DKO results
edgeR_result_files<-dir("../../RNA-seq/DEG/data/",pattern="*FC2_fdr05.tsv",full.names=TRUE)
edgeR_result_files<-edgeR_result_files[!grepl("MLL",edgeR_result_files)]
edgeR_resList<-lapply(edgeR_result_files,function(X)read.table(X,sep="\t",quote="",header=TRUE))
names(edgeR_resList)<-unlist(lapply(strsplit(basename(edgeR_result_files),split="\\."),function(X)X[2]))

#read in single and double KO results
edgeR_result_files<-dir("../../RNA-seq/DEG/data/",pattern="*FC2_fdr05.tsv",full.names=TRUE)
edgeR_result_files<-edgeR_result_files[grepl("MLL",edgeR_result_files)]
edgeR_resList<-c(edgeR_resList,lapply(edgeR_result_files,function(X)read.table(X,sep="\t",quote="",header=TRUE)))
names(edgeR_resList)[3:5]<-paste0("GV_",gsub("\\.FC2","",unlist(lapply(strsplit(basename(edgeR_result_files),split="_"),function(X)X[7]))))


set_custom_wd<-function(label){
  wdir<-file.path(gwdir,"MA_Volcano_plots",label)
  system(paste0('mkdir -p ', wdir))
  setwd(wdir)
  getwd()
}

recalc_change_cats<-function(table){
  table$change_recalc<-"none"
  table$change_recalc[table$change %in% "UP"]<-"LFC>=1 and FDR<0.05"
  table$change_recalc[table$change %in% "DOWN"]<-"LFC<=-1 and FDR<0.05"
  table$change_recalc[table$logFC>0 & table$logFC<1 & table$FDR<0.05]<-"1>LFC>0 and FDR<0.05"
  table$change_recalc[table$logFC<0 & table$logFC>-1 & table$FDR<0.05]<-"-1<LFC<0 and FDR<0.05"
  return(table)
}

PLOT_BASE_SIZE <- 18
PLOT_THEME <- function(){
  theme_classic(base_size = PLOT_BASE_SIZE) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.margin = margin(6, 6, 6, 6),
      axis.title = element_text(size = PLOT_BASE_SIZE, face = "bold"),
      axis.text = element_text(size = PLOT_BASE_SIZE-2),
      plot.title = element_text(hjust = 0.5, size = PLOT_BASE_SIZE+2, face = "bold"),
      plot.margin = margin(10, 10, 10, 10),
      aspect.ratio = 1,
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      axis.line = element_blank()
    )
}
CHANGE_LEVELS <- c(
  "LFC>=1 and FDR<0.05",
  "LFC<=-1 and FDR<0.05",
  "1>LFC>0 and FDR<0.05",
  "-1<LFC<0 and FDR<0.05",
  "none"
)
get_change_color_map <- function(){
  c(
    'LFC>=1 and FDR<0.05' = '#D15FEE',
    'LFC<=-1 and FDR<0.05' = '#4575B4',
    '1>LFC>0 and FDR<0.05' = colorspace::lighten('#D15FEE',amount=0.5),
    '-1<LFC<0 and FDR<0.05' = colorspace::lighten("#4575B4",amount=0.5),
    'none' = 'grey'
  )
}

add_group_counts <- function(p, table, x_pos){
  colmap <- get_change_color_map()
  cnt <- table(table$change_recalc)
  up_key   <- "LFC>=1 and FDR<0.05" # leaving only those
  down_key <- "LFC<=-1 and FDR<0.05" # leaving only those
  n_up   <- ifelse(up_key   %in% names(cnt), as.integer(cnt[up_key]),   0)
  n_down <- ifelse(down_key %in% names(cnt), as.integer(cnt[down_key]), 0)
  p <- p +
    annotate("text",
             x = x_pos, y = 10,
             label = paste("Up:", n_up),
             color = colmap[up_key],
             hjust = 0, vjust = 1, size = 5.5) +
    annotate("text",
             x = x_pos, y = -9,
             label = paste("Down:", n_down),
             color = colmap[down_key],
             hjust = 0, vjust = 1, size = 5.5)
  return(p)
}


plot_custom_MA<-function(table,label){
  table<-recalc_change_cats(table)
  table$change_recalc <- factor(table$change_recalc, levels = CHANGE_LEVELS)
  
  if(grepl("MLL",label)){
    if(grepl("MLL34",label)){
      label_to_title<-"Mll3/4 cKO"
    }else{label_to_title<-sprintf("%s cKO",stringr::str_to_title(gsub("GV_","",label)))}
  }else{
    label_to_title<-sprintf("%s stage oocyte",gsub("_.+","",label))
  }
  colmap <- get_change_color_map()
  p.ma <- ggplot(table, aes(logCPM, logFC)) +
    geom_point(aes(color = change_recalc), show.legend = TRUE, size = 1, alpha = 1) +  # Adjust point size and transparency
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at log2FoldChange = 0
    scale_color_manual(values = colmap, drop = FALSE, guide = guide_legend(ncol=1)) +
    labs(x = 'log2 CPM', y = 'log2 Fold Change (cKO/WT)') +
    ylim(-11, 11) +
    PLOT_THEME() + 
    ggtitle(sprintf("MA Plot\nDifferentially expressed genes in \n%s",label_to_title)) +
    geom_hline(yintercept=1,linetype = "dashed", color = "grey80", linewidth = 0.3)+
    geom_hline(yintercept=(-1),linetype = "dashed", color = "grey80", linewidth = 0.3)
  
  # Save the plot to a file
  p.ma <- add_group_counts(p.ma, table, x_pos = 0)
  ggsave(paste0(label, '.ma_plot.png'), p.ma, width = 8, height = 6)
  ggsave(paste0(label, '.ma_plot.pdf'), p.ma, width = 8, height = 6)
}

plot_custom_Volcano<-function(table,label){
  table<-recalc_change_cats(table)
  table$change_recalc <- factor(table$change_recalc, levels = CHANGE_LEVELS)
  
  if(grepl("MLL",label)){
    if(grepl("MLL34",label)){
      label_to_title<-"Mll3/4 cKO"
    }else{label_to_title<-sprintf("%s cKO",stringr::str_to_title(gsub("GV_","",label)))}
  }else{
    label_to_title<-sprintf("%s stage oocyte",gsub("_.+","",label))
  }
  colmap <- get_change_color_map()
  p.volcano <- ggplot(table, aes(logFC, -log10(FDR))) +
    geom_point(aes(color = change_recalc), show.legend = TRUE, size = 1, alpha = 1) +  # Adjust size and transparency
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical line at x = 0
    scale_color_manual(values = colmap, drop = FALSE, guide = guide_legend(ncol=1)) +  # Muted red, blue, and gray
    labs( x = 'log2 Fold Change (cKO/WT)', y = '-log10(FDR)') +
    xlim(-11, 11) +
    ylim(0, 11) +
    PLOT_THEME() + 
    ggtitle(sprintf("Volcano Plot\nDifferentially expressed genes in \n%s",label_to_title))+
    geom_vline(xintercept=1,linetype = "dashed", color = "grey80", linewidth = 0.3)+geom_vline(xintercept=(-1),linetype = "dashed", color = "grey80", linewidth = 0.3)+
    geom_hline(yintercept=-log10(0.05),linetype = "dashed", color = "grey80", linewidth = 0.3)
  
  # Save the plot to a file
  ggsave(paste0(label, '.volcano_plot.png'), p.volcano, width = 8, height = 6)
  ggsave(paste0(label, '.volcano_plot.pdf'), p.volcano, width = 8, height = 6)
  
}  


produce_results<-function(table,label){
  set_custom_wd(label)
  plot_custom_MA(table,label)
  plot_custom_Volcano(table,label)
}

mapply(FUN=function(X,Y)produce_results(X,Y),edgeR_resList,names(edgeR_resList))

setwd(gwdir)
getwd()

sink("MA_Volcano_plots/MA_Volcano_sessionInfo.txt")
sessionInfo()
sink()



###
