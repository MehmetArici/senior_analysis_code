setwd('~/Dropbox/SpongeWebsite/Table3_expression_datasets/')
filekey = "E-MTAB-2706"
filename<- paste("/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/", filekey, "-query-results_removedcomments.tsv", sep = "")
EMTAB_file = read.csv(filename, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")
exp_data = data.matrix(EMTAB_file)
save(exp_data,file ="exp_data_2706.RData")

mapping = read.table('EMTAB_mapping.txt', header = TRUE, check.names = FALSE, sep = '\t')
mapping_frame = data.table(mapping)
save(mapping_frame, file = "EMTAB_mapping_frame.RData")

mapping = read.table('EMTAB_GTEX_mapping.txt', header = TRUE, check.names = FALSE, sep = '\t')
mapping_frame = data.table(mapping)
save(mapping_frame, file = "EMTAB_GTEX_mapping_frame.RData")



filekey = "E-MTAB-2770"
filename =  paste("/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/", filekey, "_median_by_tissue.txt", sep = "")
EMTAB_file = read.csv(filename, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")
exp_data = data.matrix(EMTAB_file)
save(exp_data,file ="exp_data_2770_median.RData")


filekey = "GTEX"
filename =  paste("/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/", filekey, "_median_by_tissue.txt", sep = "")
EMTAB_file = read.csv(filename, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")
exp_data = data.matrix(EMTAB_file)
save(exp_data,file ="GTEX_median.RData")

filename =  "/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/GTEX/GTEX_mapping.txt"
GTEX_mapping = read.table(filename, header = TRUE, check.names = FALSE, sep = '\t')

filename = "/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/base_data/expression_data/GTEX/GTEx_depthmatched_onlyENSG_rounded.gct"
GTEX_full_data = read.csv(filename, header = TRUE, row.names=1, check.names=FALSE, sep = "\t")

unique_tissues = unique(GTEX_mapping$tissue)

for (tissue_name in unique_tissues){
    #filename = paste('GTEX_', tissue_name, '.RData', sep = '')
    samples = GTEX_mapping[GTEX_mapping$tissue == tissue_name,]$sample_id;  
    print(paste(tissue_name, length(samples)), sep = '\t')
    #exp_data = GTEX_full_data[,samples];
    #save(exp_data, file = filename);
}



