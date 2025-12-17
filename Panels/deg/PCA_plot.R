library(ggplot2)

#set working directory to source file directory as all paths are relative to it
gwdir<-getSrcDirectory(function(){})[1]
setwd(gwdir)
getwd()

# for GO and GV stage oocyte, produce one PCA plot,each
#for MLL3/4 single and double knockouts should be merged into one matrix and put on one PCA plot
# all cpm matrices should be log2-transformed with pseudocount of 0.1

cpm_files<-dir("../../RNA-seq/DEG/data/",pattern="cpm.id_v*",full.names=TRUE)
cpm_files<-cpm_files[!grepl("MLL",cpm_files)]
cpm_List<-lapply(cpm_files,function(X)read.table(X,sep="\t",quote="",header=TRUE))
names(cpm_List)<-unlist(lapply(strsplit(basename(cpm_files),split="\\."),function(X)X[3]))
#sanity check if there are any NAs
lapply(cpm_List,function(X)sum(is.na(X)))

MLL_cpm_files<-dir("../../RNA-seq/DEG/data/",pattern="cpm.id_v*",full.names=TRUE)
MLL_cpm_files<-MLL_cpm_files[grepl("MLL",MLL_cpm_files)]
MLL_cpm_List<-lapply(MLL_cpm_files,function(X)read.table(X,sep="\t",quote="",header=TRUE))
names(MLL_cpm_List)<-paste0("GV_",gsub("\\.WT","",unlist(lapply(strsplit(basename(MLL_cpm_files),split="_"),function(X)X[4]))))


#mapply is completely useless in this case
MLL_cpm_renamed_cols<-vector("list",length(MLL_cpm_List))
names(MLL_cpm_renamed_cols)<-names(MLL_cpm_List)
for(i in seq_along(MLL_cpm_List)){
  MLL_cpm_renamed_cols[[i]]<-MLL_cpm_List[[i]]
  colnames(MLL_cpm_renamed_cols[[i]])[-1]<-paste0(names(MLL_cpm_List)[i],colnames(MLL_cpm_renamed_cols[[i]])[-1])
}

MLL_cpm_merged<-Reduce(function(...) merge(..., all=T,by="id"), MLL_cpm_renamed_cols)
#kick out rows with NA values as they will cause issues for PC calculation
MLL_cpm_merged_noNA<-MLL_cpm_merged[complete.cases(MLL_cpm_merged),]

cpm_List[[3]]<-MLL_cpm_merged_noNA
names(cpm_List)[3]<-"MLL3_MLL4_MLL3.4"

log_cpm_List<-lapply(cpm_List,function(X){
                              log2(X[,-1]+0.1)
})

var_List<-lapply(log_cpm_List,function(X)as.numeric(apply(X,1,var)))
#sanity check if all variance values are defined
lapply(var_List,function(X)sum(is.na(X)))

var_List_top500<-lapply(var_List,function(X){
                  names(X)<-1:length(X)
                  names(sort(X,decreasing = TRUE)[1:500])
})

lapply(var_List_top500,function(X)sum(is.na(X)))

pca_input<-mapply(FUN=function(X,Y){
  rownames(X)<-1:nrow(X)
  X[Y,]
},log_cpm_List,var_List_top500)
lapply(pca_input,function(X)sum(is.na(X)))

#define a vector with number of groups to map PCA plotting function to
ngroupv<-c(2,2,4)

set_custom_wd<-function(label){
  wdir<-file.path(gwdir,"PCA_plots",label)
  system(paste0('mkdir -p ', wdir))
  setwd(wdir)
  getwd()
}

calculate_plot_PCA_2groups<-function(table,label){
  set.seed(123)
  res.pca = prcomp(t(table),center=TRUE,scale=FALSE)
  plot_data<-as.data.frame(res.pca$x)
  plot_data$Group_tmp<-gsub("_.+","",rownames(plot_data))
  plot_data$Group<-"Control"
  plot_data$Group[plot_data$Group_tmp %in% "KO"]<-"Mll3/4 cKO"
  plot_data$Group<-factor(plot_data$Group,levels=unique(plot_data$Group))
  z<-summary(res.pca)
  perc_var_explained<-z$importance['Proportion of Variance',1:2]
  ggplot(plot_data)+geom_point(aes(x=PC1,y=PC2,colour=Group),size=5)+
  xlab(paste0("PC1: ",round(perc_var_explained[1]*100),"% variance"))+ylab(paste0("PC2: ",round(perc_var_explained[2]*100),"% variance"))+
  scale_colour_manual(values = c("black", "#D15FEE"))+
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + # Horizontal dotted line
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") + # Vertical dotted line
  theme_classic(base_size = 12) + # Clean background
  theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      aspect.ratio=1,plot.title = element_text(hjust = 0.5)
    )+ ggtitle (sprintf("PCA plot on log2 CPMs from  \n Control and Mll3/4 cKO %s stage samples.",gsub("_.+","",label)))
  ggsave(paste0(label,"_PCA.png"))
  ggsave(paste0(label,"_PCA.pdf"))
  
}

calculate_plot_PCA_4groups<-function(table,label){
  set.seed(123)
  res.pca = prcomp(t(table),center=TRUE,scale=FALSE)
  plot_data<-as.data.frame(res.pca$x)
  plot_data$Group_tmp<-gsub("_rep[1-3]","",rownames(plot_data))
  plot_data$Group_tmp<-gsub("GV_","",plot_data$Group_tmp)
  plot_data$Group<-"Control"
  plot_data$Group[plot_data$Group_tmp %in% "MLL3KO"]<-"Mll3 cKO"
  plot_data$Group[plot_data$Group_tmp %in% "MLL4KO"]<-"Mll4 cKO"
  plot_data$Group[plot_data$Group_tmp %in% "MLL34KO"]<-"Mll3/4 cKO"
  plot_data$Group<-factor(plot_data$Group,levels=unique(plot_data$Group))
  z<-summary(res.pca)
  perc_var_explained<-z$importance['Proportion of Variance',1:2]
  ggplot(plot_data)+geom_point(aes(x=PC1,y=PC2,colour=Group),size=5)+
    xlab(paste0("PC1: ",round(perc_var_explained[1]*100),"% variance"))+ylab(paste0("PC2: ",round(perc_var_explained[2]*100),"% variance"))+
    scale_colour_manual(values = c("black","orange","darkblue", "#D15FEE"))+
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") + # Horizontal dotted line
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") + # Vertical dotted line
    theme_classic(base_size = 12) + # Clean background
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      aspect.ratio=1,plot.title = element_text(hjust = 0.5)
    ) + ggtitle ("PCA plot on log2 CPMs from Control, \n Mll3 cKO, Mll4 cKO and Mll3/4 cKO samples.")
  ggsave(paste0(label,"_PCA.png"))
  ggsave(paste0(label,"_PCA.pdf"))
  
}



produce_results<-function(table,label,ngroups){
  set_custom_wd(label)
  if(ngroups==2){calculate_plot_PCA_2groups(table,label)
                 }else{calculate_plot_PCA_4groups(table,label)}
}

mapply(FUN=function(X,Y,Z)produce_results(X,Y,Z),pca_input,names(pca_input),ngroupv)

#set working directory to source file directory as all paths are relative to it
gwdir<-getSrcDirectory(function(){})[1]
setwd(gwdir)
getwd()

sink("PCA_plots/PCA_sessionInfo.txt")
sessionInfo()
sink()

