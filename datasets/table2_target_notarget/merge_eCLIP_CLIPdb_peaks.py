dir = '/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/'

list_files = [dir + 'eCLIP/RBP_CLIP_peaks_on_genes_eCLIP_sponge.txt',  
              dir + 'CLIPdb/RBP_CLIP_peaks_on_genes_CLIPdb_sponge.txt']

fout = open (dir  + 'RBP_CLIP_peaks_on_genes_eCLIP_and_CLIPdb_merged.txt' , 'w')
#fout_list_rbps = open (add + 'list_all_RBPs.txt' , 'w')

dict_sites = {}
list_transcript_ids = []
list_rbps = []


for file in list_files:
    finput = open (file)
    print(file)
    for line in finput:
        words = line.strip('\n').split('\t')
        enst_id = words[0]
        if enst_id not in dict_sites:
            dict_sites[enst_id] = []
            list_transcript_ids.append(enst_id)
        for i in range (1, len(words)):
            rbp = words[i].split(':')[0].split('_')[0]
            if rbp not in list_rbps:
                list_rbps.append(rbp)
            dict_sites[enst_id].append(words[i])


for enst in dict_sites:
    fstr = enst
    for item in dict_sites[enst]:
        fstr += '\t' + item
    fout.write(fstr + '\n')
fout.close()

'''
for rbp in list_rbps:
    fout_list_rbps.write(rbp + '\n')
fout_list_rbps.close()
'''