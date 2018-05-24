library(ggplot2)

file <- "/media/DATAPART1/lncRNAs_Human/base_data/expression_data/GTEX/GTEx_depthmatched_onlyENSG.gct"
GTEX_file = read.table(file, header = TRUE, row.names=1, check.names=FALSE)

exp_data = data.matrix(GTEX_file)

list_lncRNA_mapping = list("NORAD" = "ENSG00000260032", "CRNDE"="ENSG00000245694" , "BANCR"="ENSG00000278910", "FENDRR"="ENSG00000268388", "DANCR"="ENSG00000226950", "NEAT1"="ENSG00000245532", "TUG1"="ENSG00000253352", "HOTAIR"="ENSG00000228630", "PVT1"="ENSG00000249859", "RMST"="ENSG00000255794", "ST8SIA3"="ENSG00000267225", "TINCR"="ENSG00000223573", "RP11-686F15.3"="ENSG00000260030", "RP3-337D23.3"="ENSG00000196634", "H19"="ENSG00000130600" , "GAS5"="ENSG00000234741", "SNHG19"="ENSG00000260260", "RP11"="ENSG00000275216", "CTB"="ENSG00000266469")
list_RBP_mapping = list("ELAVL1"="ENSG00000066044" , "HNRNPC"="ENSG00000092199", "PUM1"="ENSG00000134644" , "PUM2"="ENSG00000055917", "IGF2BP2"="ENSG00000073792", "IGF2BP1"="ENSG00000159217", "TIA1"="ENSG00000116001", "KHSRP"="ENSG00000088247")

tissue_list = readLines("/media/DATAPART1/lncRNAs_Human/base_data/expression_data/GTEX/GTEX_tissue_subset_names.txt")

rbp = "PUM2"
lncRNA = "NORAD"

if (rbp == "PUM2"){
  target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/targets_removedENCzeros.txt" , sep = ""))
  nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/nontargets_removedENCzeros.txt" , sep = ""))
} else {
  target_list = readLines(paste("", sep = ""))
  nontarget_list = readLines(paste("", sep = ""))
}

count = 0
list_dataframes = list()
for (tissue in tissue_list){
  count = count + 1
  print (paste(count, tissue, sep = " "))
  file_tissue <- paste ("/media/DATAPART1/lncRNAs_Human/base_data/expression_data/GTEX/tissue_specific_samples/", tissue , ".txt", sep = "")
  tissue_samples = readLines(file_tissue)
  common_samples = intersect(tissue_samples, colnames(exp_data))
  
  if (length(common_samples) != 0){
    
    exp_data_tissue = subset(exp_data, select=c(common_samples))
    
    existing_target_list = intersect(target_list, row.names(exp_data_tissue))
    num_targets = length(existing_target_list)
    
    existing_nontarget_list = intersect(nontarget_list, row.names(exp_data_tissue))
    num_nontargets = length(existing_nontarget_list)
    
    lncRNA_data = exp_data_tissue[list_lncRNA_mapping[[lncRNA]],2:ncol(exp_data_tissue)]
    
    target_corrs = list()
    for (target in existing_target_list){
      target_data = exp_data_tissue[target,2:ncol(exp_data_tissue)]
      corr_value = cor(target_data, lncRNA_data, method = 'spearman')
      target_corrs <- c(target_corrs, corr_value)
    }
    
    nontarget_corrs = list()
    for (nottarget in existing_nontarget_list){
      nontarget_data = exp_data_tissue[nottarget,2:ncol(exp_data_tissue)]
      corr_value = cor(nontarget_data, lncRNA_data, method = 'spearman')
      nontarget_corrs <- c(nontarget_corrs, corr_value)
    }
    
    a = do.call(rbind.data.frame, target_corrs)
    b = do.call(rbind.data.frame, nontarget_corrs)
    dataframe_a <- data.frame(group = "Targets", tissue, value = a[,1])
    dataframe_b <- data.frame(group = "Background", tissue, value = b[,1])
    plot.data.tissue <- rbind(dataframe_a,dataframe_b)
    list_dataframes[[count]] <- plot.data.tissue
    #list_dataframes[[count]] <- dataframe_a   
    
  }
}

plot.data.total = do.call(rbind, list_dataframes)

plot.data.total$tissue<-reorder(plot.data.total$tissue,-plot.data.total$value)
ggplot(data = plot.data.total, aes(x=tissue, y=value)) + geom_boxplot(aes(fill=group)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#ggplot(plot.data, aes(x=group, y=value, fill= group)) + scale_fill_brewer(palette="Set1") + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.5, label= signif(pvalue, digits = 6))

filename = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_lncRNAdata/tissue_specific/all_in_one/', lncRNA, '_vs_PUM_target_nontarget_GTEX_Alltissues.jpeg' , sep= '')
ggsave(filename, plot = last_plot())
