rm(list=ls(all=TRUE))
library(ggplot2)
library(DAAG)
# devtools::install_github('gokceneraslan/DAAG')


corr_RBP_lncRNA <- function(existing_list, rbp_data, lncRNA_data){
  num_genes = length(existing_list)
  list_corrs_RBP_CV = numeric(length = num_genes)
  list_corrs_RBPlncRNA_CV = numeric(length = num_genes)
  i = 1
  for (gene in existing_list){
    gene_data = exp_data[gene,skipcols:ncol(exp_data)];

    if (RBP == "PUM2")
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
    list_corrs_RBP_CV[i] <- round(corr_RBP_CV, digits = 5);
    corr_RBPlncRNA_CV = cor(cvres_RBP_lncRNA$cvpred, gene_data, method = "spearman");
    list_corrs_RBPlncRNA_CV[i] <-  round(corr_RBPlncRNA_CV,5);
    i = i + 1;
  }

  return (list(list_corrs_RBP_CV, list_corrs_RBPlncRNA_CV))
}




args = commandArgs(trailingOnly=TRUE)
lncRNA = args[1]
RBP = args[2]
dataset = args[3]
job_id = args[4]

if(length(args) > 4) # user supplied target  file
{
  targetfile = args[4]
}

infofileh <- file(paste(job_id, '_analysis2_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis2_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

targetdir = '/var/www/html/analysis_code/datasets/table2_target_notarget/results/'
datadir = '/var/www/html/analysis_code/datasets/Table3_expression_datasets/'
outdir = '/var/www/html/senior/public/results/analysis2_regression/'
setwd('/var/www/html/analysis_code/analysis2_regression/')

#lncRNA = "ENSG00000260032" # we expect ENSG id
#RBP = "PUM2"# we expect gene symbol
#dataset = "EMTAB2770"
#dataset = "GTEX_Liver"

load(paste(datadir, 'EMTAB_GTEX_mapping_frame.RData', sep = ''))

if (!(RBP %in%  mapping_frame$symbol))
{
  print("ERROR: There is no expression data for this RBP")
} else if(!(lncRNA %in%  mapping_frame$ensg_id)){
  print("ERROR: There is no expression data for this lncRNA")
} else {

  if(dataset == 'EMTAB2706'){
    skipcols = 2
    filename = paste(datadir, 'exp_data_2706.RData', sep = '')
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }
  } else if(dataset == 'EMTAB2770'){
    skipcols = 2
    filename = paste(datadir, 'exp_data_2770.RData', sep = '')
    if (file.exists(filename)){
      load(filename)
    }else{
      print("ERROR: cannot open input data file")
    }

  } else if(startsWith(dataset, 'GTEX')){
    skipcols = 1
    filename = paste(datadir, dataset, '.RData', sep = '')
    if (file.exists(filename)){
      load(filename)
      exp_data = data.matrix(exp_data)
      cat('exp data is read')
    }else{
      print("ERROR: cannot open input data file")
    }

  }


  RBP_ENSG_id = toString(mapping_frame[mapping_frame$symbol == RBP,]$ensg_id)

  if (RBP == 'PUM2'){
    PUM1_ENSG_id = 'ENSG00000134644'
    PUM1_data = exp_data[PUM1_ENSG_id,skipcols:ncol(exp_data)];
  }


  if(length(args) > 4){
    if (file.exists(targetfile)){
      target_list = readLines(targetfile)
    }
    else{
      print("ERROR: user supplied target file do not exist")
    }
  } else {
    target_list = readLines(paste(targetdir, RBP ,'_target_genes.txt', sep = ''))
  }


  existing_target_list = intersect(target_list, row.names(exp_data))
  num_targets = length(existing_target_list)

  #existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
  #num_nontargets = length(existing_nontarget_list)
  lncRNA_data = exp_data[lncRNA,skipcols:ncol(exp_data)]
  rbp_data = exp_data[RBP_ENSG_id,skipcols:ncol(exp_data)];
  result_target = corr_RBP_lncRNA(existing_target_list, rbp_data, lncRNA_data);
  target_corrs_RBP = result_target[[1]];
  target_corrs_RBPlncRNA = result_target[[2]];


  #result_nontarget = corr_RBP_lncRNA(existing_nontarget_list, rbp_data, lncRNA_data);
  #nontarget_corrs_RBP = result_nontarget[[1]];
  #nontarget_corss_RBPlncRNA = result_nontarget[[2]];


  #targets VS RBPlncRNA
  dataframe_a1 <- data.frame(group = "RBP only" , value = target_corrs_RBP);
  dataframe_b1 <- data.frame(group = "RBP + lncRNA" , value = target_corrs_RBPlncRNA);
  plot.data1 <- rbind(dataframe_a1,dataframe_b1);

  data.print = data.frame('gene_id' = existing_target_list,'RBP_only_corr' = target_corrs_RBP,'RBP_and_lncRNA_corr' = target_corrs_RBPlncRNA)

  outfilekey = paste(outdir, RBP, '_', lncRNA, '_analysis2_', dataset, sep = '')

  write.table(data.print, file = paste(outfilekey, '.txt', sep = ''),  quote = FALSE, sep = "\t", row.names = FALSE)

  result = wilcox.test(target_corrs_RBP, target_corrs_RBPlncRNA);
  pvalue <- result$p.value;


  ggplot(plot.data1, aes(x=group, y=value, fill= group)) + scale_fill_brewer(palette="Dark2") + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") +  annotate("text", x=1.5, y=-0.2, label= paste('p-val: ', signif(pvalue, digits = 3)))

  filename1 = paste(outfilekey, '.jpeg' , sep= '');
  ggsave(filename1, plot = last_plot())

}
sink(type="message")
sink(type="output")




