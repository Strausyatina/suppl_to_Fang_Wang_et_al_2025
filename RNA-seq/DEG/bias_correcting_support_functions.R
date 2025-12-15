## =========
## GET DATA:
## =========
get_sample_genes_counts <- function(counts_df, spike_name = "ERCC"){
  ### F: Removes genes which names start with a particular word (spike_name); 
  ### counts_df should have 1st column with gene_id. 
  return(counts_df[!startsWith(counts_df[, 1], spike_name), ])
}
counts_paired_pooling <- function(counts_df){
  ### This funclion is needed only if you have a paired structure (like initial seq and deep re-seq)
  ### names should be "GROUP_REP_SEMIGROUP", it is expected that semigroups will be in pairs (like seq and re-seq)
  data.frame(
    counts_df$gene_id,
    do.call(cbind, lapply(1:floor(ncol(counts_df)/2), function(i){
      rowSums(counts_df[, (2*i):(2*i+1)])
    }))
  ) %>% setNames(c("gene_id", 
                   names(counts_df)[seq(2,ncol(counts_df),2)] %>% strsplit("_") %>% 
                     lapply(., function(x) paste(x[1:2], collapse='_')) %>% unlist()
  ))
}
get_dgList <- function(counts_df, group_grep = c("WT", "KO"), thr_cpm = 1, n_cpm = NA){
  if(is.na(n_cpm)) n_cpm = floor((ncol(counts_df)-1) / 2)
  sampleType = rep(0, ncol(counts_df) -1) 
  sampleType[grep(group_grep[2], colnames(counts_df[, -1]))] = 1
  dgList = DGEList(counts=counts_df[, -1], genes=counts_df[, 1],
                    group = sampleType)
  countsPerMillion = cpm(dgList)
  ### We keep only those genes that are represented at least thr_cpm reads in at least n_cpm samples (default, n/2):
  countCheck = countsPerMillion > thr_cpm
  dgList = dgList[which(rowSums(countCheck) >= n_cpm),]
  return(list(dgList=dgList, cpm=cpm(dgList)))
}
#
## =========
## FIT DATA:
## =========
select_housekeeping <- function(dgList, cpm_mean_thr = 100){
  ### F: Select genes with low variance.
  cpm_stats_df = data.frame(gene_id = dgList$genes$genes,
                            cpm_mean = cpm(dgList) %>% rowMeans(),
                            cpm_variance = apply(cpm(dgList), 1, var)) 
  cpm_stats_df$low_var = (cpm_stats_df$cpm_mean > cpm_mean_thr & cpm_stats_df$cpm_mean > cpm_stats_df$cpm_variance)
  genes_stable = cpm_stats_df[cpm_stats_df$low_var, "gene_id"]
  return(list(cpm_stats_df=cpm_stats_df, genes_stable=genes_stable))
}
fit_housekeeping <- function(dgList, cpm_mean_thr = 100, genes_stable = NA, 
                             correction_mode = "housekeeping_geommean", robustGLM = F){
  ### (!) correction_mode %in% c("housekeeping_size", "housekeeping_geommean")
  if(anyNA(genes_stable)) {
    selection = select_housekeeping(dgList, cpm_mean_thr)
    genes_stable = selection$genes_stable
    cpm_stats_df = selection$cpm_stats_df
  } else {
    selection = select_housekeeping(dgList, cpm_mean_thr)
    cpm_stats_df = selection$cpm_stats_df
  }
  # Adjusting:
  libsizes = dgList$counts  %>% colSums()
  if(correction_mode == "housekeeping_size"){
    psALL_H = (dgList$counts[dgList$genes$genes %in% genes_stable, ]  %>% colSums()) / libsizes
    house_coefficient = psALL_H / ((psALL_H %>% sort())[round(ncol(dgList$counts)/2)])
    dgList$samples$norm.factors = house_coefficient
  } else if(correction_mode == "housekeeping_geommean"){
    sampleMeans = edgeR::cpm(dgList, log = TRUE, prior.count = 3) %>% colMeans()
    overallMean = mean(sampleMeans)
    ref_column = which.min(abs(sampleMeans - overallMean))
    #
    counts_stable = dgList$counts[dgList$genes$genes %in% genes_stable, ]
    log2_sums = counts_stable %>% log2() %>% colSums()
    n = counts_stable %>% nrow()
    x_term = - 1/n * log2_sums[ref_column] + log2(libsizes[ref_column])
    house_coefficient = 2**(1/n * log2_sums - log2(libsizes) + x_term)
    dgList$samples$norm.factors = house_coefficient
  } else{
    return()
  }
  #
  designMat <- model.matrix(~ dgList$samples$group)
  if(!robustGLM) {
    dgList = estimateDisp(dgList, designMat)
  } else {
    dgList = estimateGLMRobustDisp(dgList, designMat)
  # ^ additional variance shrinkage, if needed
  }
  res_fit = glmQLFit(dgList, designMat)
  res_lrt = glmQLFTest(res_fit) 
  # ^ reasoning: we have small number of samples, so glmQLFTest(); otherwise limma+voom might be better
  #
  pltMDS = plotMDS(dgList)
  pltBCF = plotBCV(dgList)
  #
  return(list(
    dgList = dgList,
    fit = res_fit,
    lrt = res_lrt,
    genes_low_var = genes_stable,
    cpm_stats_df = cpm_stats_df,
    plots = list(pltMDS=pltMDS, pltBCF=pltBCF)
  ))
}
#
## =====
## PLOT:
## =====
table_for_de_plots <- function(res_fit, fdr_thr = 0.001, log2FC_thr = 1){
  res = res_fit$lrt$table
  res$neglog10P = -log10(res$PValue)
  res$FDR = p.adjust(res$PValue, method = "BH")
  res$FCthreshold = 2**log2FC_thr
  res$significance = fdr_thr
  res$significant = res$FDR <= fdr_thr & abs(res$logFC) >= log2FC_thr
  res$change = "NONE"
  res$change[res$logFC > 0 & res$significant] = "UP"
  res$change[res$logFC < 0 & res$significant] = "DOWN"
  return(data.frame(id = sapply(strsplit(res_fit$lrt$genes$genes, "[.]"), function(x) x[1]), 
                    id_v = res_fit$lrt$genes$genes, 
                    res))
}
plot_de_plots <- function(tab_de, thr_cpm = 1){
  nonsignificant_stats =
    do.call(rbind, list(tab_de[!tab_de$significant, "logFC"] %>% summary(),
                        tab_de[tab_de$logCPM > log2(thr_cpm) & !tab_de$significant, "logFC"] %>% summary(),
                        tab_de[tab_de$logCPM > log2(thr_cpm) & !tab_de$significant & abs(tab_de$logFC) < tab_de$FCthreshold, "logFC"] %>% summary()))
  
  gg_hist = tab_de %>% 
    ggplot(aes(logFC, group = significant, fill = significant), alpha = 0.5) +
    geom_histogram(binwidth = 0.05, alpha = 1) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c("black", "red"))+
    facet_grid(significant ~ .) +
    scale_x_continuous(breaks = seq(-5,5,1), limits = c(-5,5)) +
    labs(x = "log2 (FC)")
  
  gg_volc = tab_de %>% 
    ggplot(aes(logFC, -log10(FDR), col = significant), alpha = 0.5) +
    geom_point() +
    geom_hline(yintercept = -log10(tab_de$significance[1]), lty = "dashed") +
    geom_vline(xintercept = c(-1,1) * log2(tab_de$FCthreshold[1]), 
               lty = "dashed") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("black", "red"))+
    labs(y = "-log10 (FDR)", x = "log2 (FC)")
  
  gg_pine = tab_de %>% 
    ggplot(aes(logCPM, logFC)) +
    geom_point(size=0.5, aes(col = significant)) +
    ylim(-max(abs(tab_de$logFC)),max(abs(tab_de$logFC))) +
    theme_minimal() +
    scale_color_manual(values = c("black", "red")) +
    theme(legend.position = "bottom")
  
  gg_pine_dense = cowplot::plot_grid(
    tab_de %>% 
      ggplot(aes(logCPM, logFC)) +
      geom_pointdensity(size=0.5) +
      xlim(-0.5, max(tab_de$logCPM)) +
      ylim(-max(abs(tab_de$logFC)),max(abs(tab_de$logFC))) +
      theme_minimal() +
      labs(col = "point density", y = "log2 (FC)", x = "log2 (CPM)") +
      theme(legend.position = "bottom"),
    tab_de %>% 
      ggplot(aes(logFC, group = significant, col = significant)) +
      stat_bin(aes(y=..density..), geom="step", binwidth = 0.1, position="identity") +
      xlim(-max(abs(tab_de$logFC)),max(abs(tab_de$logFC))) +
      labs(x = "log2 (FC)") +
      theme_minimal() +
      labs(col = "significance") +
      theme(legend.position = "bottom") +
      # facet_grid(~ (FDR<0.05)) +
      coord_flip() +
      scale_color_manual(values = c("black", "red")),
    rel_widths = c(1, 0.4), align = "h"
  )
  
  library(car)
  gg_hist_pval = tab_de %>% ggplot(aes(PValue)) + geom_histogram(binwidth = 0.02) + theme_minimal()
  gg_QQ_pval   = qqPlot(tab_de$PValue)
  
  change_stats = table(tab_de[, c("significant", "change")])
  
  return(list(
    nonsignificant_stats = nonsignificant_stats,
    change_stats = change_stats,
    gg_hist = gg_hist,
    gg_volc = gg_volc,
    gg_pine = gg_pine,
    gg_pine_dense = gg_pine_dense,
    gg_hist_pval = gg_hist_pval,
    gg_QQ_pval = gg_QQ_pval
  ))
}
#
## ==============================
## ONE FUNCTION TO RULE THEM ALL:
## ==============================
perform_de_grouped_norm <- function(counts_df_lst, 
                             spike_name = "ERCC", group_grep = c("WT", "KO"), 
                             thr_cpm = 1, n_cpm = NA, cpm_mean_thr = 100,
                             fdr_thr = 0.001, log2FC_thr = 1, 
                             correction_mode = "housekeeping_geommean",
                             pair_pool = F, robustGLM = F){
  ### (!) should be: length(group_grep) == length(counts_df_lst)
  if (pair_pool){
    pooled_counts_lst = lapply(counts_df_lst, function(counts_df)
      counts_df %>% get_sample_genes_counts(spike_name=spike_name) %>% counts_paired_pooling()
    )
  } else {
    pooled_counts_lst = lapply(counts_df_lst, function(counts_df)
      counts_df %>% get_sample_genes_counts(spike_name=spike_name) 
    )
  }  
  dgList_lst = lapply(pooled_counts_lst, function(pooled_counts) 
    pooled_counts %>% get_dgList(group_grep=group_grep, thr_cpm=thr_cpm, n_cpm=n_cpm)
  )
  
  genes_stable = lapply(dgList_lst, function(dgList) select_housekeeping(dgList$dgList, cpm_mean_thr=cpm_mean_thr))
  genes_stable_merged  = Reduce(intersect, lapply(genes_stable, function(x) x$genes_stable))
  
  res_data_lst = lapply(dgList_lst, function(dgList){
    fit_data = dgList$dgList %>% fit_housekeeping(cpm_mean_thr=cpm_mean_thr, genes_stable=genes_stable_merged, 
                                                  correction_mode=correction_mode, robustGLM=robustGLM)
    res_tab = fit_data %>% table_for_de_plots(fdr_thr=fdr_thr, log2FC_thr=log2FC_thr) 
    res_plt = res_tab %>% plot_de_plots()
    return(list(
      fit_data=fit_data,
      res_tab=res_tab,
      res_plt=res_plt
    ))
  })
  
  return(lapply(1:length(counts_df_lst), function(i)
    list(counts_df = counts_df_lst[[i]],
         counts_df_sample = pooled_counts_lst[[i]],
         cpm_df= data.frame(id = res_data_lst[[i]]$res_tab$id, dgList_lst[[i]]$cpm),
         genes_stable = genes_stable[[i]]$genes_stable,
         genes_stable_merged = genes_stable_merged,
         data_fit = res_data_lst[[i]]$fit_data,
         tab_plt = res_data_lst[[i]]$res_tab,
         res = res_data_lst[[i]]$res_plt)
  ))
}
