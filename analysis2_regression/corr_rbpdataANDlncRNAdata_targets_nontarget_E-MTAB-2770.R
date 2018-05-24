library(ggplot2)
library(DAAG)


corr_PUMs_lncRNA <- function(existing_list, PUM1_data, PUM2_data, lncRNA_data){
  #list_corrs_withformula = list()
  #list_corrs_resid = list()
  #list_corrs_withrbpdata = list()
  list_corrs_PUMs_CV = list()
  #list_corrs_PUM2_CV = list()
  list_corrs_PUMslncRNA_CV = list()
  #list_corrs_PUM2lncRNA_CV = list()
  for (gene in existing_list){
    gene_data = exp_data[gene,2:ncol(exp_data)]
    #gene_data_scaled = scale(gene_data, center = TRUE)
    
    #data_PUM2 = data.frame(gene_data, PUM2_data
    #cvres_PUM2 = cv.lm(data_PUM2, formula(gene_data ~ PUM2_data), m=10)
    #corr_PUM2_CV = cor(cvres_PUM2$cvpred, gene_data, method = "spearman")
    #list_corrs_PUM2_CV <- c(list_corrs_PUM2_CV, corr_PUM2_CV)
    
    data_PUMs = data.frame(gene_data, PUM2_data, PUM1_data)  
    cvres_PUMs = cv.lm(data_PUMs, formula(gene_data ~ PUM2_data + PUM1_data), m=10)
    corr_PUMs_CV = cor(cvres_PUMs$cvpred, gene_data, method = "spearman")
    list_corrs_PUMs_CV <- c(list_corrs_PUMs_CV, corr_PUMs_CV)
    
    data_PUMs_lncRNA = data.frame(gene_data, PUM2_data, PUM1_data, lncRNA_data)
    cvres_PUMs_lncRNA = cv.lm(data_PUMs_lncRNA, formula(gene_data ~ PUM2_data + PUM1_data + lncRNA_data), m=10)
    corr_PUMslncRNA_CV = cor(cvres_PUMs_lncRNA$cvpred, gene_data, method = "spearman")
    list_corrs_PUMslncRNA_CV <- c(list_corrs_PUMslncRNA_CV, corr_PUMslncRNA_CV)
    
    #data_PUM2_lncRNA = data.frame(gene_data, PUM2_data, lncRNA_data)
    #cvres_PUM2_lncRNA = cv.lm(data_PUM2_lncRNA, formula(gene_data ~ PUM2_data + lncRNA_data), m=10)
    #corr_PUM2lncRNA_CV = cor(cvres_PUM2_lncRNA$cvpred, gene_data, method = "spearman")
    #list_corrs_PUM2lncRNA_CV <- c(list_corrs_PUM2lncRNA_CV, corr_PUM2lncRNA_CV)
    
    #corr_value_scaled_withformula = cor(target_data_scaled, formula_scaled, method = 'spearman') 
    #corr_value_scaled_withrbpdata = cor(target_data_scaled, PUM2_data_scaled, method = 'spearman')
    #lm =  lm(formula = PUM2_data_scaled ~  lncRNA_data_scaled)
    #corr_resid_value = cor(target_data_scaled, resid(lm), method = 'spearman' )
    
    #list_corrs_withrbpdata <- c(list_corrs_withrbpdata, corr_value_scaled_withrbpdata)  
    #list_corrs_resid <- c(list_corrs_resid, corr_resid_value)  
    #list_corrs_withformula <- c(list_corrs_withformula, corr_value_scaled_withformula)
  }
  
  return (list(list_corrs_PUMs_CV, list_corrs_PUMslncRNA_CV))
}



file_2770 <- "/media/DATAPART1/lncRNAs_Human/base_data/expression_data/E-MTAB-2770-query-results_removedcomments.tsv"
file_2706 <- "/media/DATAPART1/lncRNAs_Human/base_data/expression_data/E-MTAB-2706-query-results_removedcomments.tsv"

EMTAB2770_file = read.csv(file_2770, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")
EMTAB2706_file = read.csv(file_2706, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")


list_lncRNA_mapping = list("NORAD" = "ENSG00000260032", "CRNDE"="ENSG00000245694" , "BANCR"="ENSG00000278910", "FENDRR"="ENSG00000268388", "DANCR"="ENSG00000226950", "NEAT1"="ENSG00000245532", "TUG1"="ENSG00000253352", "HOTAIR"="ENSG00000228630", "PVT1"="ENSG00000249859", "RMST"="ENSG00000255794", "ST8SIA3"="ENSG00000267225", "TINCR"="ENSG00000223573", "RP11-686F15.3"="ENSG00000260030", "RP3-337D23.3"="ENSG00000196634", "H19"="ENSG00000130600" , "GAS5"="ENSG00000234741", "SNHG19"="ENSG00000260260", "RP11"="ENSG00000275216", "CTB"="ENSG00000266469", "LRRC75A"="ENSG00000175061")
list_RBP_mapping = list("ELAVL1"="ENSG00000066044" , "HNRNPC"="ENSG00000092199", "PUM1"="ENSG00000134644" , "PUM2"="ENSG00000055917", "IGF2BP2"="ENSG00000073792", "IGF2BP1"="ENSG00000159217", "TIA1"="ENSG00000116001", "KHSRP"="ENSG00000088247")


target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/targets_removedENCzeros.txt" , sep = ""))
existing_target_list = intersect(target_list, row.names(exp_data))
num_targets = length(existing_target_list)

nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/nontargets_removedENCzeros.txt" , sep = ""))
existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
num_nontargets = length(existing_nontarget_list)


file = "2706"

if (file == "2770"){
  exp_data = data.matrix(EMTAB2770_file)
} else {
  exp_data = data.matrix(EMTAB2706_file)
}



PUM2_data = exp_data[list_RBP_mapping[["PUM2"]],2:ncol(exp_data)]
#PUM2_data_scaled = scale(PUM2_data, center = TRUE)

PUM1_data = exp_data[list_RBP_mapping[["PUM1"]],2:ncol(exp_data)]
#PUM1_data_scaled = scale(PUM1_data, center = TRUE)

lncRNA = "NORAD"

lncRNA_data = exp_data[list_lncRNA_mapping[[lncRNA]],2:ncol(exp_data)]
#lncRNA_data_scaled = scale(lncRNA_data, center = TRUE)


#alpha = 0.14 #scaled data
#formula_scaled = PUM2_data_scaled - (alpha * lncRNA_data_scaled)

result_target = corr_PUMs_lncRNA(existing_target_list, PUM1_data, PUM2_data, lncRNA_data)
target_corrs_PUMs = result_target[[1]]
target_corss_PUMslncRNA = result_target[[2]]


result_nontarget = corr_PUMs_lncRNA(existing_nontarget_list, PUM1_data, PUM2_data, lncRNA_data)
nontarget_corrs_PUMs = result_nontarget[[1]]
nontarget_corss_PUMslncRNA = result_nontarget[[2]]



#targets VS PUMSlncRNA
a1 = do.call(rbind.data.frame, target_corrs_PUMs)
b1 = do.call(rbind.data.frame, target_corss_PUMslncRNA)
dataframe_a1 <- data.frame(group = "Targets PUMs" , value = a1[,1])
dataframe_b1 <- data.frame(group = "Targets PUMs+lncRNA" , value = b1[,1])
plot.data1 <- rbind(dataframe_a1,dataframe_b1)

result = wilcox.test(a1[,1], b1[,1])
pvalue <- result$p.value

ggplot(plot.data1, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
filename1 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_10foldCV/', 'target_PUMslncRNA_', lncRNA, '_10foldCV_', file , '.pdf' , sep= '')
ggsave(filename1, plot = last_plot())



#Targets PUMslncRNA VS. NonTargets PUMslncRNA 
a2 = do.call(rbind.data.frame, target_corss_PUMslncRNA)
b2 = do.call(rbind.data.frame, nontarget_corss_PUMslncRNA)
dataframe_a2 <- data.frame(group = "Targets PUMs+lncRNA" , value = a2[,1])
dataframe_b2 <- data.frame(group = "Background PUMs+lncRNA" , value = b2[,1])
plot.data2 <- rbind(dataframe_a2,dataframe_b2)

result = wilcox.test(a2[,1], b2[,1])
pvalue <- result$p.value

ggplot(plot.data2, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
filename2 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_10foldCV/', 'target_nontarget_PUMslncRNA_', lncRNA, '_10foldCV_', file , '.pdf' , sep= '')
ggsave(filename2, plot = last_plot())


# #plot diff between corrs of target and nontargets with formula 
# a1 = do.call(rbind.data.frame, target_corrs_withformula)
# b1 = do.call(rbind.data.frame, nontarget_corrs_withformula)
# dataframe_a1 <- data.frame(group = "Targets" , value = a1[,1])
# dataframe_b1 <- data.frame(group = "Background" , value = b1[,1])
# plot.data1 <- rbind(dataframe_a1,dataframe_b1)
# 
# result = wilcox.test(a1[,1], b1[,1])
# pvalue <- result$p.value
# 
# ggplot(plot.data1, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
# 
# filename1 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_rbpdata_minus_alpha_lncRNAdata/', rbp, '_targets_nontargets', '_rbpdataANDlncRNAdata_scaled_', 'E-MTAB-2770' , '.jpeg' , sep= '')
# ggsave(filename1, plot = last_plot())
# 
# 
# 
# #Plot diff between target corrs and residual target corrs
# a2 = do.call(rbind.data.frame, target_corrs_withrbpdata)
# b2 = do.call(rbind.data.frame, target_corrs_resid)
# dataframe_a2 <- data.frame(group = "Targets Corr RBPdata" , value = a2[,1])
# dataframe_b2 <- data.frame(group = "Targets Corr Resid" , value = b2[,1])
# plot.data2 <- rbind(dataframe_a2,dataframe_b2)
# 
# result = wilcox.test(a2[,1], b2[,1])
# pvalue <- result$p.value
# 
# ggplot(plot.data2, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
# filename2 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_rbpdata_minus_alpha_lncRNAdata/', rbp, '_targets_', 'rbpdataANDlncRNAdata_scaled_resid_', 'E-MTAB-2770' , '.jpeg' , sep= '')
# ggsave(filename2, plot = last_plot())
# 
# 
# #Plot diff between nontarget corrs and residual target corrs
# a3 = do.call(rbind.data.frame, nontarget_corrs_withrbpdata)
# b3 = do.call(rbind.data.frame, nontarget_corrs_resid)
# dataframe_a3 <- data.frame(group = "NonTargets Corr RBPdata" , value = a3[,1])
# dataframe_b3 <- data.frame(group = "NonTargets Corr Resid" , value = b3[,1])
# plot.data3 <- rbind(dataframe_a3,dataframe_b3)
# 
# result = wilcox.test(a3[,1], b3[,1])
# pvalue <- result$p.value
# 
# ggplot(plot.data3, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
# filename3 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_rbpdata_minus_alpha_lncRNAdata/', rbp, '_nontargets_', 'rbpdataANDlncRNAdata_scaled_resid_', 'E-MTAB-2770' , '.jpeg' , sep= '')
# ggsave(filename3, plot = last_plot())


