import gzip
import os


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


source = 'eCLIP'
dir = '/media/mllab/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/'

fgene_poses = open (dir + 'ensembl_grch37_protcoding_longest_enst_3utr_seqs_sponge.txt')
fgene_poses.readline()


files_list = os.listdir(dir + "/eCLIP/data/fdr_corrected")

metadata_file = open(dir + '/eCLIP/metadata_added2RBPs.tsv','r')
header = metadata_file.readline()
header_parts = header.split('\t')
#for i in range (0,len(header_parts)):
    #print(i, header_parts[i])

fout = open('RBP_CLIP_peaks_on_genes_eCLIP_sponge.txt','w')


dict_metadata = {}
for line in metadata_file:
    words = line.split('\t')
    exp_id = words[0]
    print words
    cell = words[6]
    
    RBP = words[16].split('-')[0] + '_' + words[28] #biological replicate
    if words[16].split('-')[0]  == 'RBM5' or words[16].split('-')[0]  == 'HNRNPL':
       RBP = words[16].split('-')[0] + '_' + words[29] #biological replicate
    dict_metadata[exp_id] = [RBP, cell]


dict_clip = {}
for file_name in files_list:
    if '.bed' in file_name or 'bed.final' in file_name: #bed.final are those that we added after: 2 RBPs
        experiment_id = file_name[:file_name.find('.')]
        RBP_id = dict_metadata[experiment_id][0]
        
        if True: #RBP_id == 'PUM2_1' or RBP_id == 'PUM2_2': #True: # debugging
            print 'RBP id is ', RBP_id
            cell = dict_metadata[experiment_id][1]
            clip_file = open(dir + '/eCLIP/data/fdr_corrected/' + file_name ,'r')
            print clip_file 
            if RBP_id not in dict_clip:
                dict_clip[RBP_id] = {}
            for line in clip_file:
              words = line.split()
              if 'random' not in line and words[0] != 'chrUn' and not words[1].startswith('gl'):     
                cchrom = words[0]
                cspos = int(words[1])
                cepos = int(words[2])
                cstrand = words[5]
                cscore = str(round(float(words[6]),3))
                #cpvalue = float(words[7])
                cfdr = str(round(float(words[8]),3))
                #if cscore >= 3 and cfdr <= 0.05:
                if cchrom not in dict_clip[RBP_id].keys():
                    dict_clip[RBP_id][cchrom] = [[cstrand,cspos,cepos,cell,cscore,cfdr]]
                else:
                    dict_clip[RBP_id][cchrom].append([cstrand,cspos,cepos,cell,cscore,cfdr])
   
print dict_clip      

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
                  peak_strand = peak[0]
                  if peak_strand == gstrand:
                      peak_spos = peak[1]
                      peak_epos = peak[2]
                      cell = peak[3]
                      cscore = peak[4]
                      cfdr = peak[5]
                      cumul_spos = 0	
                      for pos_pair in list_exons:
                          exon_spos = pos_pair[0]
                          exon_epos = pos_pair[1]
                          #print peak_spos, peak_epos, exon_spos, exon_epos

                          (peak_local_spos, peak_local_epos) = intersect_peak_gene(peak_spos, peak_epos, exon_spos, exon_epos, peak_strand)
                          if not (peak_local_spos == -1 and peak_local_epos == -1):
                              peak_str = RBP + ':' + str(cumul_spos + peak_local_spos) + ':' + str(cumul_spos + peak_local_epos) + ':' + source + ':' + cell + ':' + cscore + ':' + cfdr
                              #print peak_str
                              if enst_id == 'ENST00000215832':
                                  print exon_spos, exon_epos, peak_spos, peak_epos
                                  print peak_str
                                  print cumul_spos
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
