library(ggplot2)
library(DAAG)
# devtools::install_github('gokceneraslan/DAAG')

setwd('/Users/User/Dropbox/SpongeWebsite/analysis_code/regression/')

corr_RBP_lncRNA <- function(existing_list, rbp_data, lncRNA_data){
  list_corrs_RBP_CV = list();
  list_corrs_RBPlncRNA_CV = list();
  for (gene in existing_list){
    gene_data = exp_data[gene,2:ncol(exp_data)];
    
    if (rbp == "PUM2")
    {
      data_RBP = data.frame(gene_data, rbp_data, PUM1_data);
      cvres_RBP = cv.lm(data_RBP, formula(gene_data ~ rbp_data + PUM1_data), m=10, plotit = FALSE, printit = FALSE);
      data_RBP_lncRNA = data.frame(gene_data, rbp_data, lncRNA_data, PUM1_data);
      cvres_RBP_lncRNA = cv.lm(data_RBP_lncRNA, formula(gene_data ~ rbp_data + PUM1_data + lncRNA_data), m=10,plotit = FALSE, printit = FALSE);
  
    }
    else
    {
       data_RBP = data.frame(gene_data, rbp_data)  ;
       cvres_RBP = cv.lm(data_RBP, formula(gene_data ~ rbp_data), m=10, plotit = FALSE, printit = FALSE);
       data_RBP_lncRNA = data.frame(gene_data, rbp_data, lncRNA_data);
       cvres_RBP_lncRNA = cv.lm(data_RBP_lncRNA, formula(gene_data ~ rbp_data + lncRNA_data), m=10, plotit = FALSE, printit = FALSE);
    }
    corr_RBP_CV = cor(cvres_RBP$cvpred, gene_data, method = "spearman");
    list_corrs_RBP_CV <- c(list_corrs_RBP_CV, corr_RBP_CV);
    corr_RBPlncRNA_CV = cor(cvres_RBP_lncRNA$cvpred, gene_data, method = "spearman");
    list_corrs_RBPlncRNA_CV <- c(list_corrs_RBPlncRNA_CV, corr_RBPlncRNA_CV);
  }
  
  return (list(list_corrs_RBP_CV, list_corrs_RBPlncRNA_CV))
}



file_2770 <- "/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/E-MTAB-2770-query-results_removedcomments.tsv"
file_2706 <- "/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/E-MTAB-2706-query-results_removedcomments.tsv"

setwd('/Users/User/Dropbox/SpongeWebsite/analysis_code/regression/')
EMTAB2770_file = read.csv(file_2770, header = TRUE, row.names=1, check.names=FALSE, sep = "\t");
EMTAB2706_file = read.csv(file_2706, header = TRUE, row.names=1, check.names=FALSE, sep = "\t");

file = "2770";
if (file == "2770"){
  exp_data = data.matrix(EMTAB2770_file);
} else {
  exp_data = data.matrix(EMTAB2706_file);
}


rbp = "PUM2";

PUM1_data = exp_data[list_RBP_mapping[["PUM1"]],2:ncol(exp_data)];

#if (rbp == "PUM2"){
#  target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/targets_removedENCzeros.txt" , sep = ""))
#  nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/nontargets_removedENCzeros.txt" , sep = ""))
#} else {
#  target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v11_approach/", rbp, "_target_genes.txt" , sep = ""))
#  nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v11_approach/", rbp, "_nontarget_genes.txt" , sep = ""))
#}
target_list = readLines(paste('/Users/User/Dropbox/SpongeWebsite/Table2_target_notarget/results/PUM2_target_genes.txt', sep = ''))
nontarget_list = readLines(paste('/Users/User/Dropbox/SpongeWebsite/Table2_target_notarget/results/PUM2_nontarget_genes.txt', sep=''))

existing_target_list = intersect(target_list, row.names(exp_data));
num_targets = length(existing_target_list);

existing_nontarget_list = intersect(nontarget_list, row.names(exp_data));
num_nontargets = length(existing_nontarget_list);

rbp_data = exp_data[list_RBP_mapping[[rbp]],2:ncol(exp_data)];


#lncRNA = "GAS5"
#lncRNA_list = list("H19", "MALAT1", "SNHG1", "SNHG6", "LRRC75A", "SNHG19")
lncRNA_list = list("NORAD");
list_lncRNA_mapping = list("NORAD" = "ENSG00000260032", "CRNDE"="ENSG00000245694" , "BANCR"="ENSG00000278910", "FENDRR"="ENSG00000268388", "DANCR"="ENSG00000226950", "NEAT1"="ENSG00000245532", "TUG1"="ENSG00000253352", "HOTAIR"="ENSG00000228630", "PVT1"="ENSG00000249859", "RMST"="ENSG00000255794", "ST8SIA3"="ENSG00000267225", "TINCR"="ENSG00000223573", "RP11-686F15.3"="ENSG00000260030", "RP3-337D23.3"="ENSG00000196634", "H19"="ENSG00000130600" , "GAS5"="ENSG00000234741", "SNHG19"="ENSG00000260260", "RP11"="ENSG00000275216", "CTB"="ENSG00000266469", "LRRC75A"="ENSG00000175061", "MALAT1"="ENSG00000251562", "SNHG1"="ENSG00000255717", "SNHG6"="ENSG00000245910", "SNHG19"="ENSG00000260260", "SNHG19"="ENSG00000260260");
list_RBP_mapping = list("ELAVL1"="ENSG00000066044" , "HNRNPC"="ENSG00000092199", "PUM1"="ENSG00000134644" , "PUM2"="ENSG00000055917", "IGF2BP2"="ENSG00000073792", "IGF2BP1"="ENSG00000159217", "TIA1"="ENSG00000116001", "KHSRP"="ENSG00000088247");

for (lncRNA in lncRNA_list){
  
  print (lncRNA);
  lncRNA_data = exp_data[list_lncRNA_mapping[[lncRNA]],2:ncol(exp_data)];
  #lncRNA_data_scaled = scale(lncRNA_data, center = TRUE)
  
  #alpha = 0.14 #scaled data
  #formula_scaled = PUM2_data_scaled - (alpha * lncRNA_data_scaled)
  
  result_target = corr_RBP_lncRNA(existing_target_list, rbp_data, lncRNA_data);
  target_corrs_RBP = result_target[[1]];
  target_corss_RBPlncRNA = result_target[[2]];
  
  
  result_nontarget = corr_RBP_lncRNA(existing_nontarget_list, rbp_data, lncRNA_data);
  nontarget_corrs_RBP = result_nontarget[[1]];
  nontarget_corss_RBPlncRNA = result_nontarget[[2]];
  

  #targets VS RBPlncRNA
  a1 = do.call(rbind.data.frame, target_corrs_RBP);
  b1 = do.call(rbind.data.frame, target_corss_RBPlncRNA);
  dataframe_a1 <- data.frame(group = "Targets RBP" , value = a1[,1]);
  dataframe_b1 <- data.frame(group = "Targets RBP+lncRNA" , value = b1[,1]);
  plot.data1 <- rbind(dataframe_a1,dataframe_b1);
  
  result = wilcox.test(a1[,1], b1[,1]);
  pvalue <- result$p.value;
  
  ggplot(plot.data1, aes(x=group, y=value, fill= group)) + scale_fill_brewer(palette="Dark2") + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
  filename1 = paste('./results/', lncRNA, '_vs_', rbp, '_target_10foldCV_', file , '_top3_kmers.jpeg' , sep= '');
  ggsave(filename1, plot = last_plot())
  
  
  
  #Targets RBPlncRNA VS. NonTargets RBPlncRNA 
  a2 = do.call(rbind.data.frame, target_corss_RBPlncRNA)
  b2 = do.call(rbind.data.frame, nontarget_corss_RBPlncRNA)
  dataframe_a2 <- data.frame(group = "Targets RBP+lncRNA" , value = a2[,1])
  dataframe_b2 <- data.frame(group = "Background RBP+lncRNA" , value = b2[,1])
  plot.data2 <- rbind(dataframe_a2,dataframe_b2)
  
  result = wilcox.test(a2[,1], b2[,1])
  pvalue <- result$p.value
  
  ggplot(plot.data2, aes(x=group, y=value, fill=group)) + scale_fill_brewer(palette="Accent") + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6));
  filename2 = paste('./results/', lncRNA, '_vs_', rbp, '_target_nontarget_10foldCV_', file , '_top3_kmers.jpeg' , sep= '')
  ggsave(filename2, plot = last_plot())
  
  
  #Targets
  a3 = do.call(rbind.data.frame, target_corrs_RBP);
  b3 = do.call(rbind.data.frame, target_corss_RBPlncRNA);
  c3 = do.call(rbind.data.frame, nontarget_corrs_RBP);
  d3 = do.call(rbind.data.frame, nontarget_corss_RBPlncRNA);
  dataframe_a3 <- data.frame(group = "Targets vs RBP" , value = a3[,1]);
  dataframe_b3 <- data.frame(group = "Targets vs RBP&lncRNA" , value = b3[,1]);
  dataframe_c3 <- data.frame(group = "Background vs RBP" , value = c3[,1]);
  dataframe_d3 <- data.frame(group = "Background vs RBP&lncRNA" , value = d3[,1]);
  plot.data3 <- rbind(dataframe_a3, dataframe_b3, dataframe_c3, dataframe_d3);
  
  ggplot(plot.data3, aes(x=group, y=value, fill=group)) + scale_fill_brewer(palette="Accent") + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient");
  filename3 = paste('./results/', lncRNA, '_vs_', rbp, '_target_nontarget_10foldCV_4plots_', file , '_top3_kmers.jpeg' , sep= '');
  ggsave(filename3, plot = last_plot());
  
  cat("\014");
  graphics.off();
}
