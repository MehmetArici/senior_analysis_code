library(ggplot2)

file_2706 <- "/media/DATAPART1/lncRNAs_Human/base_data/expression_data/E-MTAB-2706-query-results_removedcomments.tsv"
EMTAB2706_file = read.csv(file_2706, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")

exp_data = data.matrix(EMTAB2706_file)

list_lncRNA_mapping = list("NORAD" = "ENSG00000260032", "CRNDE"="ENSG00000245694" , "BANCR"="ENSG00000278910", "FENDRR"="ENSG00000268388", "DANCR"="ENSG00000226950", "NEAT1"="ENSG00000245532", "TUG1"="ENSG00000253352", "HOTAIR"="ENSG00000228630", "PVT1"="ENSG00000249859", "RMST"="ENSG00000255794", "ST8SIA3"="ENSG00000267225", "TINCR"="ENSG00000223573", "RP11-686F15.3"="ENSG00000260030", "RP3-337D23.3"="ENSG00000196634", "H19"="ENSG00000130600" , "GAS5"="ENSG00000234741", "SNHG19"="ENSG00000260260", "RP11"="ENSG00000275216", "CTB"="ENSG00000266469", "LRRC75A"="ENSG00000175061", "MALAT1"="ENSG00000251562", "SNHG1"="ENSG00000255717", "SNHG6"="ENSG00000245910", "SNHG19"="ENSG00000260260", "SNHG19"="ENSG00000260260")
list_RBP_mapping = list("ELAVL1"="ENSG00000066044" , "HNRNPC"="ENSG00000092199", "PUM1"="ENSG00000134644" , "PUM2"="ENSG00000055917", "IGF2BP2"="ENSG00000073792", "IGF2BP1"="ENSG00000159217", "TIA1"="ENSG00000116001", "KHSRP"="ENSG00000088247")

rbp = "ELAVL1"

if (rbp == "PUM2"){
  target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/targets_removedENCzeros.txt" , sep = ""))
  nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/nontargets_removedENCzeros.txt" , sep = ""))
} else {
  target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v11_approach/", rbp, "_target_genes.txt" , sep = ""))
  nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v11_approach/", rbp, "_nontarget_genes.txt" , sep = ""))
}


existing_target_list = intersect(target_list, row.names(exp_data))
num_targets = length(existing_target_list)

existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
num_nontargets = length(existing_nontarget_list)

#lncRNA_list = list("H19", "MALAT1", "SNHG1", "SNHG6", "LRRC75A", "TUG1", "SNHG19")
#lncRNA = "NEAT1"
lncRNA = "TUG1"

for (lncRNA in lncRNA_list){
  print (lncRNA)
  lncRNA_data = exp_data[list_lncRNA_mapping[[lncRNA]],2:ncol(exp_data)]
  
  target_corrs = list()
  for (target in existing_target_list){
    target_data = exp_data[target,2:ncol(exp_data)]
    corr_value = cor(target_data, lncRNA_data, method = 'spearman')
    target_corrs <- c(target_corrs, corr_value)
  }
  
  nontarget_corrs = list()
  for (nottarget in existing_nontarget_list){
    nontarget_data = exp_data[nottarget,2:ncol(exp_data)]
    corr_value = cor(nontarget_data, lncRNA_data, method = 'spearman')
    nontarget_corrs <- c(nontarget_corrs, corr_value)
  }
  
  a = do.call(rbind.data.frame, target_corrs)
  b = do.call(rbind.data.frame, nontarget_corrs)
  dataframe_a <- data.frame(group = "Targets" , value = a[,1])
  dataframe_b <- data.frame(group = "Background" , value = b[,1])
  plot.data <- rbind(dataframe_a,dataframe_b)
  
  result = wilcox.test(a[,1], b[,1])
  pvalue <- result$p.value
  pvalue
  
  ggplot(plot.data, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
  filename = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_lncRNAdata/', rbp, '_targets_nontargets_', 'lncRNAdata_', lncRNA, '_E-MTAB-2706' , '.jpeg' , sep= '')
  ggsave(filename, plot = last_plot())
  
}
