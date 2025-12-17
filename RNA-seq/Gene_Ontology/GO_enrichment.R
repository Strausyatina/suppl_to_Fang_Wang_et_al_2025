library("clusterProfiler")
library("org.Mm.eg.db")
library(enrichplot)
library(RColorBrewer)

#read in GV/GO DKO results
edgeR_result_files<-dir("../DEG/data/",pattern="*FC2_fdr05.tsv",full.names=TRUE)
edgeR_result_files<-edgeR_result_files[!grepl("MLL",edgeR_result_files)]
edgeR_resList<-lapply(edgeR_result_files,function(X)read.table(X,sep="\t",quote="",header=TRUE))
names(edgeR_resList)<-unlist(lapply(strsplit(basename(edgeR_result_files),split="\\."),function(X)X[2]))

set_custom_wd<-function(label){
  wdir<-file.path(label)
  system(paste0('mkdir -p ', wdir))
  setwd(wdir)
  getwd()
}


calculate_GO<-function(table,label){

##process the downregulated branch
ego_down <- enrichGO(gene     = table$id[table$change %in% "DOWN"],
                OrgDb         = org.Mm.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

ego_down_gs<-DOSE::setReadable(x=ego_down,OrgDb=org.Mm.eg.db)


#write to table
ego_down.tab = ego_down_gs@result
write.table(ego_down.tab,paste0(label,"_DOWN_eGO.table.tsv"),sep = "\t", quote = F, 
            row.names = F, col.names = T)

barplot(ego_down, showCategory=15) + viridis::scale_fill_viridis(option = "G")
ggplot2::ggsave(paste0(label,"_DOWN_eGO.barplot.png"),width=9,height=7)
ggplot2::ggsave(paste0(label,"_DOWN_eGO.barplot.pdf"),width=9,height=7)

##process the upregulated branch
ego_up <- enrichGO(gene     = table$id[table$change %in% "UP"],
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)

ego_up_gs<-DOSE::setReadable(x=ego_up,OrgDb=org.Mm.eg.db)


#write to table
ego_up.tab = ego_up_gs@result
write.table(ego_up.tab,paste0(label,"_UP_eGO.table.tsv"),sep = "\t", quote = F, 
            row.names = F, col.names = T)

barplot(ego_up, showCategory=15) + viridis::scale_fill_viridis(option = "G")
ggplot2::ggsave(paste0(label,"_UP_eGO.barplot.png"),width=9,height=7)
ggplot2::ggsave(paste0(label,"_UP_eGO.barplot.pdf"),width=9,height=7)

}


produce_results<-function(table,label){
  set_custom_wd(label)
  calculate_GO(table,label)
}

mapply(FUN=function(X,Y)produce_results(X,Y),edgeR_resList,names(edgeR_resList))

sink("GO_sessionInfo.txt")
sessionInfo()
sink()
