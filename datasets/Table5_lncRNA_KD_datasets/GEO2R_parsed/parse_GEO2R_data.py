import glob


finput_mapping_ids = open ('../../Table6_mapping/ensembl_biomart_hg19_proteincoding_gene_refseq_mapping.txt')
finput_mapping_ids.readline()

fout = open ("../NORAD_mendell_KD_LFC.txt" , "w")

dict_mapping_ids = {}
for line in finput_mapping_ids:
    words = line.strip('\n').split()
    dict_mapping_ids[words[2]] = words[0]
print dict_mapping_ids
files = glob.glob('./*.txt')

for file in files:
    print file
    dict_data = {}
    dict_pvalues = {}
    fhinput = open(file)
    header = fhinput.readline()
    for line in fhinput:
      words = line.split()
      #print words
      if len(words) > 0 and 'previous version' not in line and "NA" not in words[3]  and "NA" not in words[1]:
       gene = words[4][1:-1]
       if len(gene) >= 1:
          adj_pval = float(words[1][1:-1])
          logFC = float(words[3][1:-1])
          if gene in dict_data:
              dict_data[gene].append(logFC)
              dict_pvalues[gene].append(adj_pval)
          else:
              dict_data[gene] = [logFC]
              dict_pvalues[gene] = [adj_pval]
    outfile =  '../' + file[0:file.rfind('GEO2R')] + 'LFCs.txt'          
    fhout = open(outfile, 'w')
    print >>fhout, 'ENSG_id\tlogFC\tadj_pval\tgene'          
    for gene in dict_data:
        if gene in dict_mapping_ids:
           ENSG_id = dict_mapping_ids[gene]
           
           val = sum(dict_data[gene]) / len(dict_data[gene])
           #print gene, dict_data[gene], val
           adj_pval = max(dict_pvalues[gene])
           print >>fhout, ENSG_id + '\t' + str(val) + '\t' + str(adj_pval) + '\t' + gene
           
    fhout.close()      
        
            
    