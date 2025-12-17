#library(ggVennDiagram)
library(ggvenn)
library(ggplot2)
library(eulerr)

main_output_folder<-"RNA_manuscript_figures"
subfolder_names<-c("GV_GO","MLL3_MLL4_MLL3.4")

top_list<-vector("list",length(subfolder_names))
names(top_list)<-subfolder_names

#specify file order manually
top_list[["GV_GO"]]<-c("../../RNA-seq/DEG/data/DE_edgeR.GV_pooled.WT_cKO.FC2_fdr05.tsv","../../RNA-seq/DEG/data/DE_edgeR.GO_p14_pooled.WT_cKO.FC2_fdr05.tsv")

top_list[["MLL3_MLL4_MLL3.4"]]<-c("../../RNA-seq/DEG/data/DE_edgeR.GV.Mll_3_4_34.WT_cKO_MLL3.FC2_fdr05.tsv","../../RNA-seq/DEG/data/DE_edgeR.GV.Mll_3_4_34.WT_cKO_MLL4.FC2_fdr05.tsv","../../RNA-seq/DEG/data/DE_edgeR.GV.Mll_3_4_34.WT_cKO_MLL34.FC2_fdr05.tsv")

MLL_values<-grepl("MLL",names(top_list))


##define workflow steps###########################################################################

set_custom_wd<-function(label){
  wdir<-file.path(main_output_folder,label)
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

draw_Venn_diagram<-function(gene_id_list,MLL_value,label,direction_value){
  if(MLL_value){
    cat_names<-c("Mll3 cKO","Mll4 cKO","Mll3/4 cKO")
    fill_values<-c("gold2","lightsteelblue3","plum3")
    cat_to_title<-sprintf("%s, %s and %s samples",cat_names[1],cat_names[2],cat_names[3])
    }else{
    cat_names<-c("GV","GO")
    fill_values<-c("plum3","gold2")
    cat_to_title<-sprintf("%s and %s samples",cat_names[1],cat_names[2])
    }
  dir_to_title<-c("UP"="Upregulated","DOWN"="Downregulated")
  names(gene_id_list)<-cat_names
  p.venn<-ggvenn(gene_id_list, show_percentage = FALSE, fill_color = fill_values,stroke_color="plum4",stroke_alpha=0.5,stroke_size=0.5)
  ggsave(paste0(label,"_",direction_value,"_VennDiagram.png"), p.venn, width = 6, height = 6)
  ggsave(paste0(label,"_",direction_value,"_VennDiagram.pdf"), p.venn, width = 6, height = 6)
  p.euler<-plot(euler(gene_id_list),fill = fill_values,stroke=0.5,main=sprintf("Euler diagram: \n %s genes comparison between \n %s .",dir_to_title[direction_value],cat_to_title),
                cex.main=1,labels = list(cex = 2),edges=list("plum4"),mar=c(2,2,2,2))
  ggsave(paste0(label,"_",direction_value,"_EulerDiagram.png"), p.euler, width = 6, height = 6)
  ggsave(paste0(label,"_",direction_value,"_EulerDiagram.pdf"), p.euler, width = 6, height = 6)
  
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

sink("RNA_manuscript_figures/VennDiagram_sessionInfo.txt")
sessionInfo()
sink()
