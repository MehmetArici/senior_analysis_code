import gzip
import os
import glob

def intersect_peak_gene(peak_spos, peak_epos, gene_spos, gene_epos, pstrand):			   
   peak_local_spos = -1
   peak_local_epos = -1
   if pstrand == '+':
		if peak_spos >= gene_spos and peak_epos <= gene_epos: # peak is completely inside the gene
			peak_local_spos = peak_spos - gene_spos 
			peak_local_epos = peak_epos - gene_spos - 1 				
		elif peak_spos < gene_spos and peak_epos < gene_epos and peak_epos > gene_spos:
			peak_local_spos = 0
			peak_local_epos = peak_epos - gene_spos - 1 
			
		elif peak_spos > gene_spos and peak_spos < gene_epos and peak_epos > gene_epos:
			peak_local_spos = peak_spos - gene_spos 
			# the line below had a bug before, peak_local_epos can't be set to gene_epos as we work with local positions
			peak_local_epos = gene_epos - gene_spos - 1 			
   elif pstrand == '-':
		if peak_spos >= gene_spos and peak_epos <= gene_epos:
			peak_local_spos = gene_epos - peak_epos
			peak_local_epos = gene_epos - peak_spos - 1 
				
		elif peak_spos < gene_spos and peak_epos < gene_epos and peak_epos > gene_spos:
			peak_local_spos = gene_epos - peak_epos 
			peak_local_epos = gene_epos - gene_spos - 1 
				
		elif peak_spos > gene_spos and peak_spos < gene_epos and peak_epos > gene_epos:
			peak_local_spos = 0
			peak_local_epos = gene_epos - peak_spos - 1 # not sure about +1 or -1 
   return (peak_local_spos, peak_local_epos)			


source = 'CLIPdb'
dir = '/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/' 

fgene_poses = open (dir + 'ensembl_grch37_protcoding_longest_enst_3utr_seqs_sponge.txt')
fgene_poses.readline()


finput_list_datasets = open (dir + source + '/list_human_datasets.txt')


finput_cell = open(dir + source + '/Summary_added_cell_info.tsv')
header = finput_cell.readline()
hp = header.strip('"').split(',')


fout = open(dir + source + '/RBP_CLIP_peaks_on_genes_CLIPdb_sponge.txt','w')


clip_path = dir + '/CLIPdb/data/'

peak_bed_files = [os.path.join(dirpath, f)
              for dirpath, dirnames, files in os.walk(clip_path)
              for f in files if f.endswith('.bed')]
    
files_list = []
dict_tech = {} # a dictionary from RBP to type of peak analyzer Piranha or PARalyzer
for item in peak_bed_files:
    rbp_tech = item.split('/')[-2]
    rbp = rbp_tech.split('_')[0]
    tech = rbp_tech.split('_')[1]
    if rbp not in dict_tech:
        dict_tech[rbp] = [tech]
    else:
        if tech not in dict_tech[rbp]:
            dict_tech[rbp].append(tech)
    if tech == 'Piranha' or tech == 'PARalyzer':
        files_list.append(item)



list_human_dataset = []
for line in finput_list_datasets:
    item = line.strip('\n').split()[0]
    list_human_dataset.append(item)

    
dict_cell = {}
for line in finput_cell:
    words = line.strip('\n').split('\t')
    #words = [x.strip('"') for x in words]
    if words[1] == 'Human' and 'AGO' not in words[0]:
        rbp = words[0]
        cell = words[6]
        datasets_temp1 = words[3].split(';')
        datasets_temp2 = []
        for item in datasets_temp1:
            if ',' in item:
                datasets_temp2.extend(item.split(','))
            else:
                datasets_temp2.append(item)
        dataset = '-'.join(datasets_temp2)
        dict_cell[dataset] = cell 



# a dictionary of all CLIP peaks
dict_clip = {}
for RBP_id in dict_tech:
    if "PUM" in RBP_id:
        peak_calling_methods = dict_tech[RBP_id]
        # if both PARalyzer and Piranha are present, we only use PARalyzer
        if 'PARalyzer' in peak_calling_methods:
            tech = 'PARalyzer'
        else:
            tech = 'Piranha'
        list_files_rbp = glob.glob(dir + source + '/data/' + RBP_id + '/' + RBP_id + '_' + tech + '/*.bed')
        for file_name in list_files_rbp:
            name_part = file_name.split('/')[-1].split('_')
            dataset = name_part[1].split('.')[0]
            if dataset in list_human_dataset:
                cell = dict_cell[dataset]
                clip_file = open(file_name) 
                if RBP_id not in dict_clip:
                    dict_clip[RBP_id] = {}
                for line in clip_file:
                    words = line.strip('\n').split()
                    cchrom = words[0]
                    cspos = int(words[1])
                    cepos = int(words[2])
                    cstrand = words[5]
                    if tech == 'Piranha':
                        if len(words) >= 7:
                            sig_value = str(round(float(words[6]),3))
                        else:
                            sig_value = 'NA'
    
                    elif tech == 'PARalyzer':
                        sig_value = str(round(float(words[4]),3))
    
                    if cchrom not in dict_clip[RBP_id]:
                        dict_clip[RBP_id][cchrom] = [[cstrand,cspos,cepos,cell,sig_value]]
                    else:
                        dict_clip[RBP_id][cchrom].append([cstrand,cspos,cepos,cell,sig_value])


#genes = []
dict_sites = {}
count = 0
for line in fgene_poses:
   if '>' in line:
      words = line.strip('\n').split()
      enst_id = words[0][1:]
      count += 1
      #print(count)
      gchrom = words[3]
      gstrand = words[4]
      pos_part = words[5].split(',') # split to exons
      list_exons = []
      for pos_pair in pos_part:
          exon_spos = int(pos_pair.split(':')[0])
          exon_epos = int(pos_pair.split(':')[1])      
          list_exons.append([exon_spos, exon_epos]) 
      if gstrand == '-':
          list_exons = list_exons[::-1]    
          
      if enst_id not in dict_sites:
           dict_sites[enst_id] = []  
      for RBP in dict_clip:
          if gchrom in dict_clip[RBP]:
              for peak in dict_clip[RBP][gchrom]:        
                  if peak[0] == '.':
                    peak_strand = gstrand
                  else:
                    peak_strand = peak[0]
                  
                  if peak_strand == gstrand:
                      peak_spos = peak[1]
                      peak_epos = peak[2]
                      cell = peak[3]
                      sig_value = peak[4]
                      cumul_spos = 0	
                      for pos_pair in list_exons:
                          exon_spos = pos_pair[0]
                          exon_epos = pos_pair[1]
                          #print peak_spos, peak_epos, exon_spos, exon_epos

                          (peak_local_spos, peak_local_epos) = intersect_peak_gene(peak_spos, peak_epos, exon_spos, exon_epos, peak_strand)
                          if not (peak_local_spos == -1 and peak_local_epos == -1):
                              peak_str = RBP + ':' + str(cumul_spos + peak_local_spos) + ':' + str(cumul_spos + peak_local_epos) + ':' + source + ':' + cell + ':' + sig_value
                              #if enst_id == 'ENST00000215832':
                              #    print exon_spos, exon_epos, peak_spos, peak_epos
                              #    print peak_str
                              #    print cumul_spos
                              dict_sites[enst_id].append(peak_str)
                          cumul_spos = cumul_spos + exon_epos - exon_spos + 1 # check +1 -1    
      if dict_sites[enst_id] != []:
          fstr = enst_id + '\t' + '\t'.join(dict_sites[enst_id]) + '\n'                        
          fout.write(fstr)   
          fout.flush()       
      else:
          fout.write(enst_id + '\n')
          fout.flush()             


fout.close()
