rm(list=ls(all=TRUE))
library(ggplot2)
library(data.table)


args = commandArgs(trailingOnly=TRUE)
lncRNA = args[1]
RBP = args[2]
dataset = args[3]
job_id = args[4]

infofileh <- file(paste(job_id, '_analysis1_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis1_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

#rm(list=ls(all=TRUE))
targetdir = '/var/www/html/analysis_code/datasets/table2_target_notarget/results/'
datadir = '/var/www/html/analysis_code/datasets/Table3_expression_datasets/'

outdir = '/var/www/html/senior/public/results/analysis1_correlation/'
setwd('/var/www/html/analysis_code/analysis1_correlation/')

#lncRNA = "ENSG00000260032"
# "ENSG00000260032" # we expect ENSG id
#RBP = "PUM2"# we expect gene symbol
#dataset = 'GTEX_Liver'

load(paste(datadir, 'EMTAB_GTEX_mapping_frame.RData', sep = ''))

if(!(lncRNA %in%  mapping_frame$ensg_id)){
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
            print("file is loaded")
        }else{
            print("ERROR: cannot open input data file")
        }

    }

    if(length(args) == 5){
        target_list = readLines(args[4])
        nontarget_list = readLines(args[5])
    } else {
        target_list = readLines(paste(targetdir, RBP ,'_target_genes.txt', sep = ''))
        nontarget_list = readLines(paste(targetdir, RBP ,'_nontarget_genes.txt', sep = ''))
    }

    existing_target_list = intersect(target_list, row.names(exp_data))
    num_targets = length(existing_target_list)

    existing_nontarget_list = intersect(nontarget_list, row.names(exp_data))
    num_nontargets = length(existing_nontarget_list)
    lncRNA_data = exp_data[lncRNA,skipcols:ncol(exp_data)]

    target_corrs = numeric(length = num_targets)
    i = 1
    for (target in existing_target_list){
        target_data = exp_data[target,skipcols:ncol(exp_data)]
        corr_value = cor(target_data, lncRNA_data, method = 'spearman')
        target_corrs[i] <- round(corr_value, digits = 5)
        i = i + 1
    }
    nontarget_corrs = numeric(length = num_nontargets)
    i = 1
    for (nottarget in existing_nontarget_list){
        nontarget_data = exp_data[nottarget,skipcols:ncol(exp_data)]
        corr_value = cor(nontarget_data, lncRNA_data, method = 'spearman')
        nontarget_corrs[i] <- round(corr_value, digits = 5)
        i = i + 1
    }

    targets_df <- data.frame(group = "Targets" , corr_value = target_corrs)
    nontargets_df <- data.frame(group = "Background" , corr_value = nontarget_corrs)
    plot.data <- rbind(targets_df,nontargets_df)

    plot.data$gene_id = c(existing_target_list, existing_nontarget_list);
    outfilekey = paste(outdir, RBP, '_', lncRNA, '_analysis1_', dataset, sep = '')
    write.table(plot.data, file = paste(outfilekey, '.txt', sep = ''),  quote = FALSE, sep = "\t", row.names = FALSE)

    result = wilcox.test(target_corrs, nontarget_corrs)
    pvalue <- result$p.value


    ggplot(plot.data, aes(x=group, y=corr_value, fill= group)) + geom_boxplot() + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size = rel(1.75)), axis.title.y=element_text(size = rel(1.5), angle = 90)) + labs(y = "Spearman correlation coefficient") +
        annotate("text", x=1.5, y=-0.4, label= paste('p-val: ', signif(pvalue, digits = 3)))
    filename = paste(outfilekey, '.jpeg' , sep= '')
    ggsave(filename, plot = last_plot())


}
sink(type="message")
sink(type="output")
