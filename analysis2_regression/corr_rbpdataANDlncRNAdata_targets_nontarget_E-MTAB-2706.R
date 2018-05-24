library(ggplot2)
library(DAAG)

file <- "/Volumes/Elements/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/E-MTAB-2706-query-results_removedcomments.tsv"

EMTAB2706_file = read.csv(file, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")

exp_data = data.matrix(EMTAB2706_file)

list_lncRNA_mapping = list("NORAD" = "ENSG00000260032", "CRNDE"="ENSG00000245694" , "BANCR"="ENSG00000278910", "FENDRR"="ENSG00000268388", "DANCR"="ENSG00000226950", "NEAT1"="ENSG00000245532", "TUG1"="ENSG00000253352", "HOTAIR"="ENSG00000228630", "PVT1"="ENSG00000249859", "RMST"="ENSG00000255794", "ST8SIA3"="ENSG00000267225", "TINCR"="ENSG00000223573", "RP11-686F15.3"="ENSG00000260030", "RP3-337D23.3"="ENSG00000196634")
list_RBP_mapping = list("ELAVL1"="ENSG00000066044" , "HNRNPC"="ENSG00000092199", "PUM1"="ENSG00000134644" , "PUM2"="ENSG00000055917", "IGF2BP2"="ENSG00000073792", "IGF2BP1"="ENSG00000159217", "TIA1"="ENSG00000116001", "KHSRP"="ENSG00000088247")

#Formula: rbp_data - alpha(lncRNA_data)
rbp = "PUM2"
lncRNA = "NORAD"

rbp_data = exp_data[list_RBP_mapping[[rbp]],2:ncol(exp_data)]
rbp_data_scaled = scale(rbp_data, center = TRUE)

PUM1_data = exp_data[list_RBP_mapping[["PUM1"]],2:ncol(exp_data)]
PUM1_data_scaled = scale(PUM1_data, center = TRUE)

lncRNA_data = exp_data[list_lncRNA_mapping[[lncRNA]],2:ncol(exp_data)]
lncRNA_data_scaled = scale(lncRNA_data, center = TRUE)

#alpha = 1
alpha = 0.14 #scaled data
formula_scaled = rbp_data_scaled - (alpha * lncRNA_data_scaled)

target_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/targets_removedENCzeros.txt" , sep = ""))
existing_target_list = intersect(target_list, row.names(exp_data))
num_targets = length(existing_target_list)

nontarget_list = readLines(paste("/media/DATAPART1/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/nontargets_removedENCzeros.txt" , sep = ""))
existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
num_nontargets = length(existing_nontarget_list)

target_corrs_withformula = list()
target_corrs_resid = list()
target_corrs_withrbpdata = list()
target_corrs_PUMs_CV = list()
target_corrs_PUM2_CV = list()
target_corrs_PUMslncRNA_CV = list()
target_corrs_PUM2lncRNA_CV = list()
for (target in existing_target_list){
  target_data = exp_data[target,2:ncol(exp_data)]
  target_data_scaled = scale(target_data, center = TRUE)
  
  data_PUM2 = data.frame(target_data, rbp_data)
  cvres_PUM2 = cv.lm(data_PUM2, formula(target_data ~ rbp_data), m=10)
  corr_PUM2_CV = cor(cvres_PUM2$cvpred, target_data, method = "spearman")
  target_corrs_PUM2_CV <- c(target_corrs_PUM2_CV, corr_PUM2_CV)
  
  data_PUMs = data.frame(target_data, rbp_data, PUM1_data)  
  cvres_PUMs = cv.lm(data_PUMs, formula(target_data ~ rbp_data + PUM1_data), m=10)
  corr_PUMs_CV = cor(cvres_PUMs$cvpred, target_data, method = "spearman")
  target_corrs_PUMs_CV <- c(target_corrs_PUMs_CV, corr_PUMs_CV)
    
  data_PUMs_lncRNA = data.frame(target_data, rbp_data, PUM1_data, lncRNA_data)
  cvres_PUMs_lncRNA = cv.lm(data_PUMs_lncRNA, formula(target_data ~ rbp_data + PUM1_data + lncRNA_data), m=10)
  corr_PUMslncRNA_CV = cor(cvres_PUMs_lncRNA$cvpred, target_data, method = "spearman")
  target_corrs_PUMslncRNA_CV <- c(target_corrs_PUMslncRNA_CV, corr_PUMslncRNA_CV)
  
  data_PUM2_lncRNA = data.frame(target_data, rbp_data, lncRNA_data)
  cvres_PUM2_lncRNA = cv.lm(data_PUM2_lncRNA, formula(target_data ~ rbp_data + lncRNA_data), m=10)
  corr_PUM2lncRNA_CV = cor(cvres_PUM2_lncRNA$cvpred, target_data, method = "spearman")
  target_corrs_PUM2lncRNA_CV <- c(target_corrs_PUM2lncRNA_CV, corr_PUM2lncRNA_CV)
  
  corr_value_scaled_withformula = cor(target_data_scaled, formula_scaled, method = 'spearman') 
  corr_value_scaled_withrbpdata = cor(target_data_scaled, rbp_data_scaled, method = 'spearman')
  lm =  lm(formula = rbp_data_scaled ~  lncRNA_data_scaled)
  corr_resid_value = cor(target_data_scaled, resid(lm), method = 'spearman' )
  
  target_corrs_withrbpdata <- c(target_corrs_withrbpdata, corr_value_scaled_withrbpdata)  
  target_corrs_resid <- c(target_corrs_resid, corr_resid_value)  
  target_corrs_withformula <- c(target_corrs_withformula, corr_value_scaled_withformula)
}

target_count_signif_withrbpdata = length(which(target_corrs_withrbpdata < 0.05))
target_count_signif_resid= length(which(target_corrs_resid < 0.05))
target_count_signif_withformula = length(which(target_corrs_withformula < 0.05))

target_count_signif_PUMs_CV = length (which(target_corrs_PUMs_CV < 0.05))
target_count_signif_PUMslncRNA_CV = length (which(target_corrs_PUMslncRNA_CV < 0.05))

#print(target_count_signif_withrbpdata)
#print(target_count_signif_resid)
#print(target_count_signif_withformula)




nontarget_corrs_withformula = list()
nontarget_corrs_resid = list()
nontarget_corrs_withrbpdata = list()
nontarget_corrs_PUMs_CV = list()
nontarget_corrs_PUMslncRNA_CV = list()
for (nottarget in existing_nontarget_list){
  nontarget_data = exp_data[nottarget,2:ncol(exp_data)]
  #x <- nontarget_data[ nontarget_data != 0 ]
  #if (length(x) != 0){
  nontarget_data_scaled = scale(nontarget_data, center = TRUE)

  data_PUMs = data.frame(nontarget_data, rbp_data_scaled, PUM1_data_scaled)  
  cvres_PUMs = cv.lm(data_PUMs, formula(nontarget_data ~ rbp_data_scaled + PUM1_data_scaled), m=10)
  corr_PUMs_CV = cor(cvres_PUMs$cvpred, nontarget_data, method = "spearman")
  nontarget_corrs_PUMs_CV <- c(nontarget_corrs_PUMs_CV, corr_PUMs_CV)
  
  data_PUMs_lncRNA = data.frame(nontarget_data, rbp_data_scaled, PUM1_data_scaled, lncRNA_data_scaled)
  cvres_PUMs_lncRNA = cv.lm(data_PUMs_lncRNA, formula(nontarget_data ~ rbp_data_scaled + PUM1_data_scaled + lncRNA_data_scaled), m=10)
  corr_PUMslncRNA_CV = cor(cvres_PUMs_lncRNA$cvpred, nontarget_data, method = "spearman")
  nontarget_corrs_PUMslncRNA_CV <- c(nontarget_corrs_PUMslncRNA_CV, corr_PUMslncRNA_CV)

  corr_value_scaled_withformula = cor(nontarget_data_scaled, formula_scaled, method = 'spearman') 
  corr_value_scaled_withrbpdata = cor(nontarget_data_scaled, rbp_data_scaled, method = 'spearman')
  lm =  lm(formula = rbp_data_scaled ~  lncRNA_data_scaled)
  corr_resid_value = cor(nontarget_data_scaled, resid(lm), method = 'spearman' )
  
  nontarget_corrs_withrbpdata <- c(nontarget_corrs_withrbpdata, corr_value_scaled_withrbpdata)  
  nontarget_corrs_resid <- c(nontarget_corrs_resid, corr_resid_value)  
  nontarget_corrs_withformula <- c(nontarget_corrs_withformula, corr_value_scaled_withformula)
  #}
}

nontarget_count_signif_withrbpdata = length(which(nontarget_corrs_withrbpdata < 0.05))
nontarget_count_signif_resid= length(which(nontarget_corrs_resid < 0.05))
nontarget_count_signif_withformula = length(which(nontarget_corrs_withformula < 0.05))

nontarget_count_signif_PUMs_CV = length (which(nontarget_corrs_PUMs_CV < 0.05))
nontarget_count_signif_PUMslncRNA_CV = length (which(nontarget_corrs_PUMslncRNA_CV < 0.05))

#print(nontarget_count_signif_withrbpdata)
#print(nontarget_count_signif_resid)
#print(nontarget_count_signif_withformula)


print(target_count_signif_PUMs_CV)
print(target_count_signif_PUMslncRNA_CV)
print(nontarget_count_signif_PUMs_CV)
print(nontarget_count_signif_PUMslncRNA_CV)


#Plot diff between corrs of target and nontarget with 10 fold cross validation:
a1 = do.call(rbind.data.frame, nontarget_corrs_PUMs_CV)
b1 = do.call(rbind.data.frame, nontarget_corrs_PUMslncRNA_CV)
dataframe_a1 <- data.frame(group = "NonTargets PUMs" , value = a1[,1])
dataframe_b1 <- data.frame(group = "NonTargets PUMs + lncRNA" , value = b1[,1])
plot.data1 <- rbind(dataframe_a1,dataframe_b1)

result = wilcox.test(a1[,1], b1[,1])
pvalue <- result$p.value

ggplot(plot.data1, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))





#plot diff between corrs of target and nontargets with formula 
a1 = do.call(rbind.data.frame, target_corrs_withformula)
b1 = do.call(rbind.data.frame, nontarget_corrs_withformula)
dataframe_a1 <- data.frame(group = "Targets" , value = a1[,1])
dataframe_b1 <- data.frame(group = "Background" , value = b1[,1])
plot.data1 <- rbind(dataframe_a1,dataframe_b1)

result = wilcox.test(a1[,1], b1[,1])
pvalue <- result$p.value

ggplot(plot.data1, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))

filename1 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_rbpdata_minus_alpha_lncRNAdata/', rbp, '_targets_nontargets', '_rbpdataANDlncRNAdata_scaled_', 'E-MTAB-2706' , '.jpeg' , sep= '')
ggsave(filename1, plot = last_plot())





#Plot diff between target corrs and residual target corrs
a2 = do.call(rbind.data.frame, target_corrs_withrbpdata)
b2 = do.call(rbind.data.frame, target_corrs_resid)
dataframe_a2 <- data.frame(group = "Targets Corr RBPdata" , value = a2[,1])
dataframe_b2 <- data.frame(group = "Targets Corr Resid" , value = b2[,1])
plot.data2 <- rbind(dataframe_a2,dataframe_b2)

result = wilcox.test(a2[,1], b2[,1])
pvalue <- result$p.value

ggplot(plot.data2, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
filename2 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_rbpdata_minus_alpha_lncRNAdata/', rbp, '_targets_', 'rbpdataANDlncRNAdata_scaled_resid_', 'E-MTAB-2706' , '.jpeg' , sep= '')
ggsave(filename2, plot = last_plot())


#Plot diff between nontarget corrs and residual target corrs
a3 = do.call(rbind.data.frame, nontarget_corrs_withrbpdata)
b3 = do.call(rbind.data.frame, nontarget_corrs_resid)
dataframe_a3 <- data.frame(group = "NonTargets Corr RBPdata" , value = a3[,1])
dataframe_b3 <- data.frame(group = "NonTargets Corr Resid" , value = b3[,1])
plot.data3 <- rbind(dataframe_a3,dataframe_b3)

result = wilcox.test(a3[,1], b3[,1])
pvalue <- result$p.value

ggplot(plot.data3, aes(x=group, y=value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank()) + annotate("text", x=1.5, y=-0.4, label= signif(pvalue, digits = 6))
filename3 = paste('/media/DATAPART1/PUM_analysis/correlation/results/corr_rbpdata_minus_alpha_lncRNAdata/', rbp, '_nontargets_', 'rbpdataANDlncRNAdata_scaled_resid_', 'E-MTAB-2706' , '.jpeg' , sep= '')
ggsave(filename3, plot = last_plot())


