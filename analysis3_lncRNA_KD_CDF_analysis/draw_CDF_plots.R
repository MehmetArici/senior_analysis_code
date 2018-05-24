library(ggplot)


# read the data
setwd('/Users/User/Dropbox/SpongeWebsite/analysis_code/lncRNA_KD_CDF_analysis/')

LFC_filekey = 'NORAD_Ulitsky_KD'
LFC_filename = paste('/Users/User/Dropbox/SpongeWebsite/Table5_lncRNA_KD_datasets/', LFC_filekey, '.txt', sep = '')
LFC_data = read.table(LFC_filename, header = TRUE, row.names =1, sep = "\t")



RBP = 'PUM2'
target_list = readLines(paste('/Users/User/Dropbox/SpongeWebsite/Table2_target_notarget/results/', RBP ,'_target_genes.txt', sep = ''))
nontarget_list = readLines(paste('/Users/User/Dropbox/SpongeWebsite/Table2_target_notarget/results/', RBP ,'_nontarget_genes.txt', sep=''))


#target_list = readLines(paste("/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/targets_removedENCzeros.txt" , sep = ""))
#nontarget_list = readLines(paste("/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/analysis/separate_target_notarget_genes/results/code_v10_approach/nontargets_removedENCzeros.txt" , sep = ""))

existing_target_list = intersect(target_list, row.names(LFC_data));
num_targets = length(existing_target_list);

existing_nontarget_list = intersect(nontarget_list, row.names(LFC_data));
num_nontargets = length(existing_nontarget_list);

target_data = LFC_data[existing_target_list, 1]
nontarget_data = LFC_data[existing_nontarget_list, 1]

result = wilcox.test(target_data, nontarget_data)
pvalue <- result$p.value
print(pvalue)

df <- data.frame(x = c(target_data, nontarget_data), ggg=factor(rep(1:2, c(length(target_data),length(nontarget_data)))))
ggplot(df, aes(x, colour = ggg)) + stat_ecdf()+ scale_colour_hue(name="my legend", labels=c(paste('target',length(existing_target_list)),paste('nontarget',length(existing_nontarget_list))))

filename = paste('./results/CDF_plot_', LFC_filekey , '_top3_kmers.jpeg', sep = '')
ggsave(filename, plot = last_plot())

