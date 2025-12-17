library(ggplot2)

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
  wdir<-file.path("RNA_manuscript_figures",label)
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


plot_custom_MA<-function(table,label){
  table<-recalc_change_cats(table)
  if(grepl("MLL",label)){
    if(grepl("MLL34",label)){
      label_to_title<-"Mll3/4 cKO"
    }else{label_to_title<-sprintf("%s cKO",stringr::str_to_title(gsub("GV_","",label)))}
  }else{
    label_to_title<-sprintf("%s stage oocyte",gsub("_.+","",label))
  }
  p.ma <- ggplot(table, aes(logCPM, logFC)) +
    geom_point(aes(color = change_recalc), show.legend = TRUE, size = 1, alpha = 1) +  # Adjust point size and transparency
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Horizontal line at log2FoldChange = 0
    scale_color_manual(values = c('LFC>=1 and FDR<0.05' = '#D15FEE', 'LFC<=-1 and FDR<0.05' = '#4575B4',
                                  '1>LFC>0 and FDR<0.05'=colorspace::lighten('#D15FEE',amount=0.5),'-1<LFC<0 and FDR<0.05'=colorspace::lighten("#4575B4",amount=0.5), 'none' = 'grey'),
                       guide = guide_legend(nrow=3,ncol=2))+
    labs(x = 'logCPM', y = 'log2FoldChange') +
    ylim(-12, 12) +
    theme_classic() +
    theme(
      legend.position = "top",  # Place legend at the top
      legend.title = element_blank(),  # No legend title
      legend.key = element_blank(),  # No background for legend keys
      legend.margin = margin(6, 6, 6, 6),
      axis.title = element_text(size = 12),  # Larger axis titles
      axis.text = element_text(size = 10),  # Larger axis labels
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered, bold title
      plot.margin = margin(10, 10, 10, 10), # Add margins to the plot
      aspect.ratio=1
    ) + ggtitle(sprintf("MA Plot: Differentially expressed genes \n in %s samples.",label_to_title))+
    geom_hline(yintercept=1,linetype = "dashed", color = "grey80", linewidth = 0.3)+geom_hline(yintercept=(-1),linetype = "dashed", color = "grey80", linewidth = 0.3)
  
  # Save the plot to a file
    ggsave(paste0(label, '.ma_plot.png'), p.ma, width = 6, height = 6)
    ggsave(paste0(label, '.ma_plot.pdf'), p.ma, width = 6, height = 6)
}

plot_custom_Volcano<-function(table,label){
  table<-recalc_change_cats(table)
  if(grepl("MLL",label)){
    if(grepl("MLL34",label)){
      label_to_title<-"Mll3/4 cKO"
    }else{label_to_title<-sprintf("%s cKO",stringr::str_to_title(gsub("GV_","",label)))}
  }else{
    label_to_title<-sprintf("%s stage oocyte",gsub("_.+","",label))
  }
    p.volcano <- ggplot(table, aes(logFC, -log10(FDR))) +
    geom_point(aes(color = change_recalc), show.legend = TRUE, size = 1, alpha = 1) +  # Adjust size and transparency
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical line at x = 0
    scale_color_manual(values = c('LFC>=1 and FDR<0.05' = '#D15FEE', 'LFC<=-1 and FDR<0.05' = '#4575B4',
                                  '1>LFC>0 and FDR<0.05'=colorspace::lighten('#D15FEE',amount=0.5),'-1<LFC<0 and FDR<0.05'=colorspace::lighten("#4575B4",amount=0.5), 'none' = 'grey'),
                       guide = guide_legend(nrow=3,ncol=2)) +  # Muted red, blue, and gray
    labs( x = 'log2 Fold Change', y = '-log10(FDR)') +
    theme_classic() +
    xlim(-12, 12) +
    ylim(0, 50) +
    theme(
      legend.position = "top",  # Place the legend at the top
      legend.title = element_blank(),  # No title for the legend
      legend.key = element_blank(),  # No background for the legend keys
      legend.margin = margin(6, 6, 6, 6),
      axis.title = element_text(size = 12),  # Larger axis titles
      axis.text = element_text(size = 10),  # Larger axis labels
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered and bold title
      plot.margin = margin(10, 10, 10, 10),  # Add margins to the plot
      aspect.ratio=1
    ) + ggtitle(sprintf("Volcano Plot: Differentially expressed genes \n in %s samples.",label_to_title))+
      geom_vline(xintercept=1,linetype = "dashed", color = "grey80", linewidth = 0.3)+geom_vline(xintercept=(-1),linetype = "dashed", color = "grey80", linewidth = 0.3)+
      geom_hline(yintercept=-log10(0.05),linetype = "dashed", color = "grey80", linewidth = 0.3)
  
  p.volcano 
  
  # Save the plot to a file
  ggsave(paste0(label, '.volcano_plot.png'), p.volcano, width = 6, height = 6)
  ggsave(paste0(label, '.volcano_plot.pdf'), p.volcano, width = 6, height = 6)
  
}  


produce_results<-function(table,label){
  set_custom_wd(label)
  plot_custom_MA(table,label)
  plot_custom_Volcano(table,label)
}

mapply(FUN=function(X,Y)produce_results(X,Y),edgeR_resList,names(edgeR_resList))

sink("RNA_manuscript_figures/MA_Volcano_sessionInfo.txt")
sessionInfo()
sink()



###
