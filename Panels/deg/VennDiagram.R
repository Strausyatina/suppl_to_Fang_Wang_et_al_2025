library(ggvenn)
library(ggplot2)
library(eulerr)
library(gridExtra)

#set working directory to source file directory as all paths are relative to it
gwdir<-getwd()#getSrcDirectory(function(){})[1]
setwd(gwdir)
getwd()

main_output_folder<-"VennDiagramms"
subfolder_names<-c("GV_GO","MLL3_MLL4_MLL3.4")

top_list<-vector("list",length(subfolder_names))
names(top_list)<-subfolder_names

#specify file order manually:
top_list[["GV_GO"]]<-c("../../../../RNA-seq/DEG/data/DE_edgeR.GV_pooled.WT_cKO.FC2_fdr05.tsv","../../../../RNA-seq/DEG/data/DE_edgeR.GO_p14_pooled.WT_cKO.FC2_fdr05.tsv")

top_list[["MLL3_MLL4_MLL3.4"]]<-c("../../../../RNA-seq/DEG/data/DE_edgeR.GV.Mll_3_4_34.WT_cKO_MLL3.FC2_fdr05.tsv","../../../../RNA-seq/DEG/data/DE_edgeR.GV.Mll_3_4_34.WT_cKO_MLL4.FC2_fdr05.tsv","../../../../RNA-seq/DEG/data/DE_edgeR.GV.Mll_3_4_34.WT_cKO_MLL34.FC2_fdr05.tsv")

MLL_values<-grepl("MLL",names(top_list))


##define workflow steps###########################################################################
set_custom_wd<-function(label){
  wdir<-file.path(gwdir,main_output_folder,label)
  system(paste0('mkdir -p ', wdir))
  setwd(wdir)
  getwd()
}

get_genes<-function(file_set,MLL_value,direction_value){
  res<-lapply(file_set,function(X)read.table(X,sep="\t",quote="",header=TRUE))
  if(MLL_value){
    names(res)<-gsub("\\.FC2","",unlist(lapply(strsplit(basename(file_set),split="_"),function(X)X[7])))
  }else{
    names(res)<-unlist(lapply(strsplit(basename(file_set),split="\\."),function(X)X[2]))
  }
  id_sub<-lapply(res,function(X)X$id[X$change %in% direction_value])
  print(id_sub)
  return(id_sub)

}

draw_Euler<-function(gene_id_list,fill_values,border_values,direction_value,pdf_file){
  eu <- eulerr::euler(gene_id_list)
  # centers:
  ell <- as.data.frame(stats::coef(eu)) 
  xmin <- min(ell$h - ell$a, na.rm = TRUE)
  xmax <- max(ell$h + ell$a, na.rm = TRUE)
  x_npc <- (ell$h - xmin) / (xmax - xmin)
  # Euler plot:
  g_top <- plot(
    eu, fill = fill_values, stroke = 0.5, edges = list(col = border_values),
    quantities = list(cex = 2.5), # quant. text size
    labels = FALSE, alpha = 0.5
  )
  # labels: 
  g_bottom <- grid::grobTree(
    grid::textGrob(
      label = names(gene_id_list),
      x = grid::unit(x_npc, "npc"),
      y = grid::unit(0.5, "npc"),
      gp = grid::gpar(cex = 3) # label size
    )
  )
  # total:
  combined <- gridExtra::arrangeGrob(
    g_top,
    g_bottom,
    ncol = 1,
    heights = grid::unit.c(grid::unit(1, "npc") - grid::unit(0.18, "npc"),
                           grid::unit(0.18, "npc"))
  )
  pdf(pdf_file, width = 6, height = 3.5)
  grid::grid.newpage()
  grid::grid.draw(combined)
  dev.off()
  return(invisible(eu))
}

draw_Venn_diagram<-function(gene_id_list,MLL_value,label,direction_value){
  if(MLL_value){
    cat_names<-c("Mll3 cKO","Mll4 cKO","Mll3/4 cKO")
    fill_values<-c("gold","lightsteelblue3","plum3")
    border_values<-c("gold4","lightsteelblue4","plum4")
    cat_to_title<-sprintf("%s, %s and %s samples",cat_names[1],cat_names[2],cat_names[3])
    }else{
    cat_names<-c("GV","GO")
    fill_values<-c("plum3","gold")
    border_values<-c("plum4","gold4")
    cat_to_title<-sprintf("%s and %s samples",cat_names[1],cat_names[2])
    }
  dir_to_title<-c("UP"="Upregulated","DOWN"="Downregulated")
  names(gene_id_list)<-cat_names
  #
  # 1. Venn (size-negligent scheme):
  p.venn<-ggvenn(gene_id_list, show_percentage = FALSE, fill_color = fill_values, # dafault alpha = 0.5
                 stroke_color=border_values, stroke_alpha=0.5, stroke_size=0.5,
                 set_name_size = 9, text_size = 9)
  ggsave(paste0(label,"_",direction_value,"_VennDiagram.pdf"), p.venn, width = 6, height = 6)
  # ggsave(paste0(label,"_",direction_value,"_VennDiagram.png"), p.venn, width = 6, height = 6)
  #
  # 2. Euler (size-udjusted), running only for pairs:
  if(length(gene_id_list) == 2){
    # -- version with labels on top --
    # p.euler<-plot(euler(gene_id_list), fill = fill_values,stroke=0.5, alpha = 0.5, 
    #               cex.main=1, labels = list(cex = 2), quantities = list(cex = 2), edges=list(border_values), 
    #               mar=c(2,2,2,2))
    # ggsave(paste0(label,"_",direction_value,"_EulerDiagram.pdf"), p.euler, width = 6, height = 6)
    # ggsave(paste0(label,"_",direction_value,"_EulerDiagram.png"), p.euler, width = 6, height = 6)
    # 
    # -- version with labels below --
    draw_Euler(
      gene_id_list,fill_values,border_values,direction_value,
      pdf_file = paste0(label,"_",direction_value,"_EulerDiagram.pdf")
    )
  }
}

##################################################################################################

produce_results<-function(file_list,MLL_value){
  if(MLL_value){
    label<-"MLL3_MLL4_MLL3.4"
  }else{
    label<-"GV_GO"
  }
  set_custom_wd(label) 
  for(chdir in c("UP","DOWN")){
  g<-get_genes(file_list,MLL_value,chdir)
  draw_Venn_diagram(g,MLL_value,label,chdir)
  }
}

mapply(FUN=function(X,Y)produce_results(X,Y),top_list,MLL_values)

#gwdir<-'.' #getSrcDirectory(function(){})[1]
setwd(gwdir)

sink("VennDiagramms/VennDiagram_sessionInfo.txt")
sessionInfo()
sink()

