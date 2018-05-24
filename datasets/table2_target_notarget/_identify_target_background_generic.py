import sys
import os.path
import random 
import numpy as np
import itertools
from scipy.stats import mannwhitneyu


'''
Targets
-should have expected LFC 
-should contain the motif within the CLIP peak
CLIP: here I can be very strict with CLIPDB

Nontargets
-should not contain the motif
-should have no change after KD
-should not have any CLIP peak
-length matched

'''

def read_kd_result(finput):
    dict = {}
    # read the header
    finput.readline()
    for line in finput:
        words = line.strip("\n").split("\t")
        ensg_id_kd = words[0][1:-1] # remove the quotes
        dict[ensg_id_kd] = False
        if words[2] != 'NA' and words[6] != 'NA':
          lfc = float(words[2])
          #print line
          adj_pvalue = float(words[6])
          if ensg_id_kd.startswith("ENSG"):
            if RBP_OF_INTEREST in list_repressor:
                if lfc > KD_LFC_THRESHOLD and adj_pvalue < KD_ADJ_PVAL_THRESHOLD:
                    dict[ensg_id_kd] = True
            elif RBP_OF_INTEREST in list_activator:
                if lfc < (-1 * KD_LFC_THRESHOLD) and adj_pvalue < KD_ADJ_PVAL_THRESHOLD:
                    dict[ensg_id_kd] = True
            else: # if we don't know about the regulatory activity of the RBP we should accept for upregulated and downregulated genes
                if (lfc > KD_LFC_THRESHOLD or lfc < (-1 * KD_LFC_THRESHOLD)) and adj_pvalue < KD_ADJ_PVAL_THRESHOLD:
                    dict[ensg_id_kd] = True         
    return dict

def getKey(item):
    return item[1]

def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [_int(_random() * n) for i in itertools.repeat(None, k)]

def sample_wor(population, k):
    "Chooses k random elements (with out replacement) from a population"
    list_samples = random.sample(population, k)
    sample_indicies = []
    for item in list_samples:
        index = population.index(item)
        if index not in sample_indicies:
            sample_indicies.append(index)
        else:
            index2 = population[index+1:].index(item)
            sample_indicies.append(index2+index+1)
    #print(sample_indicies)
    return sample_indicies

PIRANHA_PVAL_THRESHOLD = 0.05
PARAYZER_SCORE_THRESHOLD = 0.7
eCLIP_PVAL_THRESHOLD = 0.05
KD_LFC_THRESHOLD = 1 
KD_ADJ_PVAL_THRESHOLD = 0.05
pval_th = 0.05
score_th = 0.7


# the following file only contains CLIP peaks with an existing motif

RBP_OF_INTEREST = sys.argv[1]
dir = '/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/'

fclip = open (dir  + 'gene_RBP_FIMO_motifs_on_eCLIP_and_CLIPdb_peaks.txt')

ffimo = open(dir + '/FIMO_files/FIMO_all3UTRs_allmotifs_0thorder_trans_based.txt')

check_KDfile_HepG2 = '../DESeq2_results/HepG2/DESeq2_results_' + RBP_OF_INTEREST + '.txt'   
# fix this after running Deseq2
check_KDfile_K562 = '../DESeq2_results/HepG2/DESeq2_results_' + RBP_OF_INTEREST + '.txt'   

finput_len = open ("/Volumes/normal/UAU/Projects/lncrnas_transfer/lncrnas_eksik/ensembl_grch37_gene_transcript_length_sorted.txt")

finput_mapping = open ("/Volumes/normal/UAU/Projects/lncrnas_transfer/PUM_analysis/base_data/ensembl_biomart_hg19_proteincoding_gene_refseq_mapping.txt")
finput_mapping.readline()

list_piranha_rbps = ['DGCR8','PTBP1PTBP2','HNRNPU','HNRNPH','RBM10','HNRNPA2B1','PTBP1','HNRNPC','HNRNPF','HNRNPM','HNRNPA1','RBM47','TIA1','EIF4A3']
list_repressor = ["TIA1", "PUM2", "PUM1", "AUF1", "ZFP36", "KHSRP", "HNRNPC"]
list_activator = ["HNRNPL", "IGF2BP2", "IGF2BP3", "ELAVLV"]



# mapping between ENSG - ENST id
dict_mapping_nmid_to_ensg = {}
dict_mapping_enst_to_ensg = {}
dict_mapping_gene_to_ensg = {}
for line in finput_mapping:
    words = line.strip("\n").split("\t")
    ensg = words[0]
    enst = words[1]
    gene = words[2]
    dict_mapping_gene_to_ensg[gene.upper()] = ensg    
    NMid = words[3]
    if NMid  != "":
        if NMid not in dict_mapping_nmid_to_ensg:
            dict_mapping_nmid_to_ensg[NMid] = ensg
    dict_mapping_enst_to_ensg[enst] = ensg

# read the length
max_len = 0
dict_length = {}
for line in finput_len:
    words = line.strip("\n").split()
    gene = words[0]
    longest_enst = words[1].split(':')[0]
    length = int(words[1].split(':')[1])
    if length > max_len:
        max_len = length
    dict_length[gene] = length



CLIP_file_all_ENSG_ids = []
dict_count_signif_peaks = {}
dict_count_motifs_in_peaks = {}
#list_rbps = []
list_ensg_ids = []
dict_count_any_peak  = {}

count = 0
for line in fclip:
    count += 1
    if count % 1000 ==  1:
        print(count)
    words = line.strip("\n").split("\t")
    enst_id = words[0]
    peaks = words[1:]
    #print line
    #print enst_id
    ensg_id = dict_mapping_enst_to_ensg[enst_id]
    dict_count_any_peak[ensg_id] = False
    if ensg_id not in CLIP_file_all_ENSG_ids:
        CLIP_file_all_ENSG_ids.append(ensg_id)
        dict_count_signif_peaks[ensg_id] = 0
        dict_count_motifs_in_peaks[ensg_id] = 0
        if ensg_id == 'ENSG00000213927':
            print 'ENST id ', enst_id
    if RBP_OF_INTEREST in line:     
        dict_count_any_peak[ensg_id] = True
        for i in range(1,len(peaks)):
            peaks_words = peaks[i].split(":")
            if '_' in peaks_words[0]:
                rbp_id = peaks_words[0].split('_')[0]
            else:
                rbp_id =  peaks_words[0]  
            if rbp_id == RBP_OF_INTEREST:
                   flag = False #initially this CLIP peak is not considered to be true positive
                   clip_info = peaks_words[3]
                   if clip_info == 'eCLIP':
                        pval = float(peaks_words[6])
                        if pval <= eCLIP_PVAL_THRESHOLD:
                           flag = True
                   elif clip_info == 'CLIPdb':
                        if rbp_id in list_piranha_rbps:
                            if peaks_words[5] != "NA":
                                pval = float(peaks_words[5])
                                if pval <= PIRANHA_PVAL_THRESHOLD:
                                    flag = True
                        else:
                            sig_score = float(peaks_words[5])
                            if sig_score >= PARAYZER_SCORE_THRESHOLD:
                                flag = True

                   if flag == True: # this is a True positive CLIP peak
                       num_motifs =  int(peaks_words[-1])  
                       if ensg_id not in dict_count_signif_peaks:
                            dict_count_signif_peaks[ensg_id] = 1
                            dict_count_motifs_in_peaks[ensg_id] = num_motifs
                       else:
                            dict_count_signif_peaks[ensg_id] += 1
                            dict_count_motifs_in_peaks[ensg_id] += num_motifs
 
print 'read the CLIP file' 
   
dict_motif = {} 
for line in ffimo:   
   words = line.split('\t')
   enst_id = words[0]
   ensg_id = dict_mapping_enst_to_ensg[enst_id]
   dict_motif[ensg_id] = 0
   for word in words[1:]:
       word_splitted = word.split(',')
       if word_splitted[0] == RBP_OF_INTEREST:
          num_sites = len(word_splitted) - 1
          dict_motif[ensg_id] = num_sites   
           
print 'read the motif file'          
     


if os.path.isfile(check_KDfile_HepG2) and os.path.isfile(check_KDfile_K562):
	dict_KD_K562 = {}
	fKD_K562 = open(check_KDfile_K562)
	dict_KD_K562 = read_kd_result(fKD_K562)
	fKD_K562.close()

	dict_KD_HepG2 = {}
	fKD_HepG2 = open(check_KDfile_HepG2)
	dict_KD_HepG2 = read_kd_result(fKD_HepG2)
	fKD_HepG2.close()		

print 'read the KD files' 

# nontarget: no peak, no motif (whether it's in the peak or not) at all, signif LFC
# target: at least one signif peak, at least one motif in the signif peak, no signif LFC


target_ids = []
nontarget_ids = []

target_ids_lens = []
nontarget_ids_lens = []

for ensg_id in CLIP_file_all_ENSG_ids:
   if ensg_id in dict_length:
        
        if dict_count_signif_peaks[ensg_id] > 0 and dict_count_motifs_in_peaks[ensg_id] > 0:
           if ensg_id in dict_KD_K562 and dict_KD_K562[ensg_id] and ensg_id in dict_KD_HepG2 and dict_KD_HepG2[ensg_id]: 
               target_ids.append(ensg_id)
               target_ids_lens.append(dict_length[ensg_id])
        if dict_count_any_peak[ensg_id] == 0 and (ensg_id not in dict_motif or dict_motif[ensg_id] == 0): # no peak,  no motif (whether it's in the peak or not) at all..
           if ensg_id in dict_KD_K562 and not dict_KD_K562[ensg_id] and ensg_id in dict_KD_HepG2 and not dict_KD_HepG2[ensg_id]: # no signif LFC change
               nontarget_ids.append(ensg_id)
               nontarget_ids_lens.append(dict_length[ensg_id])



binsize = 100
num_first_bins = (max_len // binsize) + 1
start = 0




#list_rbps = ["IGF2BP2", "IGF2BP3", "HNRNPL", "AUF1", "ZFP36", "KHSRP", "HNRNPC"]
          
        


fout_target = open ('./results/' + RBP_OF_INTEREST + '_target_genes.txt', 'w')

fout_nontarget = open ('./results/'  + RBP_OF_INTEREST + '_nontarget_genes.txt', 'w')

list_targets_sorted = []

print rbp, len(target_ids), len(nontarget_ids) 
num_targets = [0] * num_first_bins
num_nontargets = [0] * num_first_bins

dict_bins_targets = {}
for i in range (len(target_ids_lens)):
	length = target_ids_lens[i]
	bin = length // binsize
	num_targets[bin] += 1
	if bin not in dict_bins_targets:
		dict_bins_targets[bin] = [i]
	else:
		dict_bins_targets[bin].append(i)
		
		
dict_bins_nontargets = {}
for i in range (len(nontarget_ids_lens)):
	length = nontarget_ids_lens[i]
	bin = length // binsize
	num_nontargets[bin] += 1
	if bin not in dict_bins_nontargets:
		dict_bins_nontargets[bin] = [i]
	else:
		dict_bins_nontargets[bin].append(i)

aa = 0       
bb = 0
cc = 0
dd = 0
ee = 0
count_aa = 0
count_bb = 0
count_cc = 0
count_dd = 0
count_ee = 0

selected_target_indices = []
selected_nontarget_indices = []
for i in range (start, num_first_bins):
	target_lens = []
	nontarget_lens = []
	if num_targets[i] != 0:
		aa += 1
		count_aa += num_targets[i]
	if num_targets[i] != 0 and num_nontargets[i] == 0:
		bb += 1
		count_bb += num_targets[i]
	if num_targets[i] != 0 and num_nontargets[i] != 0:
		cc += 1
		count_cc += num_targets[i]
		if num_targets[i] <= num_nontargets[i]:
			dd += 1
			count_dd += num_targets[i]
			sample = sample_wor(dict_bins_nontargets[i], len(dict_bins_targets[i]))
			global_indices = []
			for s in sample:
				global_index = dict_bins_nontargets[i][s]
				nontarget_lens.append(nontarget_ids_lens[global_index])
				global_indices.append(global_index)
			selected_nontarget_indices.extend(global_indices)
			selected_target_indices.extend(dict_bins_targets[i])
			for index in dict_bins_targets[i]:
				target_lens.append(target_ids_lens[index])
		elif num_targets[i] > num_nontargets[i]:
			ee += 1
			count_ee += num_targets[i]
			sample = sample_wor(dict_bins_targets[i] , len(dict_bins_nontargets[i]))
			global_indices = []
			for s in sample:
				global_index = dict_bins_targets[i][s]
				target_lens.append(target_ids_lens[global_index])
				global_indices.append(global_index)
			selected_target_indices.extend(global_indices)
			selected_nontarget_indices.extend(dict_bins_nontargets[i])
			for index in dict_bins_nontargets[i]:
				nontarget_lens.append(nontarget_ids_lens[index])

selected_target_lens = []
selected_target_ids = []

selected_nontarget_lens = []
selected_nontarget_ids = []

for i in range(len(selected_target_indices)):
	target_id = target_ids[selected_target_indices[i]]
	nontarget_id = nontarget_ids[selected_nontarget_indices[i]]
	
	selected_target_ids.append(target_id)
	selected_nontarget_ids.append(nontarget_id)
	
	selected_target_lens.append(target_ids_lens[selected_target_indices[i]])
	selected_nontarget_lens.append(nontarget_ids_lens[selected_nontarget_indices[i]])


a, pvalue_len = mannwhitneyu(selected_target_lens, selected_nontarget_lens)
print 'p-value mannwhitneyu in terms of length ', rbp , pvalue_len

#if count_cc != len(selected_target_ids):
	#print(count_aa,count_bb,count_cc, len(selected_target_ids))

for item in selected_nontarget_ids:
	fout_nontarget.write(item + "\n")
fout_nontarget.close()

for item in selected_target_ids:
	fout_target.write(item + "\n")
fout_target.close()

