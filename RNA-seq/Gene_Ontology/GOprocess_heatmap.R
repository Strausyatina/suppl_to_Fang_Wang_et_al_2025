library(pheatmap)
library(GenomicRanges)
library(ggplot2)

wdir<-"/data/iovino/mouse_oocytes/20250602_manuscript_figures/GV_pooled"
system(paste0('mkdir -p ', wdir))
setwd(wdir)
getwd()

cpm<-read.table("/data/iovino/mouse_oocytes/DE_edgeR/cpm.id_v.GV_pooled.WT_cKO.tsv",header=TRUE)
logcpm<-log2(cpm[,-1]+0.1)
rownames(logcpm)<-cpm[,1]

GOterm<-system("grep 'regulation\\ of\\ lipid\\ metabolic\\ process' /data/iovino/mouse_oocytes/20250602_manuscript_figures/GV_pooled/GV_pooled_DOWN_eGO.table.tsv",intern=TRUE)
GOclean<-unlist(strsplit(unlist(strsplit(GOterm[1],split="\t"))[8],split="/"))

#for the list of gene symbols belonging to the chosen GO term, get ensembl gene IDs from the gtf
gtf<-rtracklayer::import("/data/repository/organisms/GRCm38_ensembl/gencode/m19/genes.gtf")
GO_eID<-unique(gtf$gene_id[match(GOclean,gtf$gene_name)])
sum(is.na(GO_eID))

heatmap_data<-logcpm[match(GO_eID,rownames(logcpm)),]
sum(is.na(heatmap_data))
rownames(heatmap_data)<-GOclean

set.seed(123)
heatmap_plot<-pheatmap(heatmap_data,
             cluster_rows = TRUE,
             lustering_method = "average",
             treeheight_row = 0,
             show_rownames = TRUE,
             cluster_cols = FALSE,
             scale = "row",
             color = colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(255),
             main = "Gene ontology: Regulation of metabolic process")$gtable

ggsave("RegulationOfLipidMetabolicProcess_heatmap.png",plot=heatmap_plot)
ggsave("RegulationOfLipidMetabolicProcess_heatmap.pdf",plot=heatmap_plot)


sink("/data/iovino/mouse_oocytes/20250602_manuscript_figures/Katarzyna_Rscripts/GOprocess_heatmap_sessionInfo.txt")
sessionInfo()
sink()



