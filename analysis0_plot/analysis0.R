rm(list=ls(all=TRUE))
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
gene_id_input = args[1]
gene_id = gene_id_input
dataset = args[2]
job_id = args[3]

infofileh <- file(paste(job_id, '_analysis0_info.txt', sep = ''), open="wt")
errorfileh <- file(paste(job_id, '_analysis0_error.txt', sep = ''), open="wt")
sink(infofileh, type="output", append=TRUE) # if this is empty, we should display 'unexpected error occurred'
sink(errorfileh, type="message", append=TRUE)

#datadir = '/Users/User/Dropbox/SpongeWebsite/Table3_expression_datasets/'
#outdir = '/Users/User/Dropbox/SpongeWebsite/analysis_code/analysis0_plot/results/'
#setwd('~/Dropbox/SpongeWebsite/analysis_code/analysis0_plot/')

datadir = '/var/www/html/analysis_code/datasets/Table3_expression_datasets/'
outdir = '/var/www/html/senior/public/results/analysis0_plot/'
setwd('/var/www/html/analysis_code/analysis0_plot/')

#gene_id = "PUM2"
#dataset = "GTEX"

load(paste(datadir, 'EMTAB_GTEX_mapping_frame.RData', sep = ''))

if(!(gene_id %in%  mapping_frame$ensg_id) && !(gene_id %in%  mapping_frame$symbol)){
    print('ERROR: this gene do not exist in expression data')
} else {
    if(!(gene_id %in%  mapping_frame$ensg_id) && (gene_id %in%  mapping_frame$symbol) )
    gene_id = toString(mapping_frame[mapping_frame$symbol == gene_id,]$ensg_id)
    if(dataset == 'EMTAB2706'){
        skip_col = 2
        filename = paste(datadir, 'exp_data_2706_median.RData', sep = '')
        if (file.exists(filename)){
            load(filename)
        }else{
            print("ERROR: cannot open input data file")
        }
    } else if(dataset == 'EMTAB2770'){
        skip_col = 2
        filename = paste(datadir, 'exp_data_2770_median.RData', sep = '')
        if (file.exists(filename)){
            load(filename)
        }else{
            print("ERROR: cannot open input data file")
        }
    } else if(dataset == 'GTEX'){
        skip_col = 1
        filename = paste(datadir, 'GTEX_median.RData', sep = '')
        if (file.exists(filename)){
            load(filename)
        }else{
            print("ERROR: cannot open input data file")
        }
    }

    outfilekey = paste(outdir, gene_id_input, '_analysis0_', dataset, sep = '')

    data = exp_data[gene_id, skip_col:ncol(exp_data)]
    sorted_data = sort(names(data), index.return = TRUE)
    data = data[sorted_data$ix]
    ylim_max = max(data) + 5
    if(dataset == 'EMTAB2706' || dataset == "GTEX"){
        df <- data.frame(tissues = names(data), expression = data) #coord_flip() +
        ggplot(df, aes(tissues, expression)) + geom_col() +  theme(text=element_text(size=16)) + coord_flip() #theme(axis.text.x=element_text(angle=90, hjust=1))
        filename = paste(outfilekey, '_.jpeg' , sep= '')
        ggsave(filename, plot = last_plot())
    } else if(dataset == 'EMTAB2770'){
        df <- data.frame(tissues = names(data[1:80]), expression = data[1:80]) #coord_flip() +
        ggplot(df, aes(tissues, expression)) + geom_col() +  theme(text=element_text(size=11)) + ylim(0, 50) + coord_flip() #theme(axis.text.x=element_text(angle=90, hjust=1))
        filename = paste(outfilekey, '_part1.jpeg' , sep= '')
        ggsave(filename, plot = last_plot())

        df <- data.frame(tissues = names(data[81:length(data)]), expression = data[81:length(data)]) #coord_flip() +
        ggplot(df, aes(tissues, expression)) + geom_col() +  theme(text=element_text(size=11)) + ylim(0, 50) + coord_flip()
        filename = paste(outfilekey, '_part2.jpeg' , sep= '')
        ggsave(filename, plot = last_plot())

    }
}
sink(type="message")
sink(type="output")