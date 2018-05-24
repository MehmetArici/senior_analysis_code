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
    dict_reverse = {}
    # read the header
    finput.readline()
    for line in finput:
        words = line.strip("\n").split("\t")
        ensg_id_kd = words[0][1:-1] # remove the quotes
        dict[ensg_id_kd] = False
        dict_reverse[ensg_id_kd] = False
        if words[2] != 'NA' and words[6] != 'NA':
          lfc = float(words[2])
          #print line
          adj_pvalue = float(words[6])
          if ensg_id_kd.startswith("ENSG"):
            if RBP_OF_INTEREST in list_repressor:
                if lfc >= KD_LFC_THRESHOLD and adj_pvalue <= KD_ADJ_PVAL_THRESHOLD:
                    dict[ensg_id_kd] = True
                if lfc <= (-1 * KD_LFC_THRESHOLD) and adj_pvalue <= KD_ADJ_PVAL_THRESHOLD:
                    dict_reverse[ensg_id_kd] = True
            elif RBP_OF_INTEREST in list_activator:
                if lfc <= (-1 * KD_LFC_THRESHOLD) and adj_pvalue <= KD_ADJ_PVAL_THRESHOLD:
                    dict[ensg_id_kd] = True
                if lfc >= KD_LFC_THRESHOLD and adj_pvalue <= KD_ADJ_PVAL_THRESHOLD:
                    dict_reverse[ensg_id_kd] = True
            else: # if we don't know about the regulatory activity of the RBP we should accept for upregulated and downregulated genes
                if (lfc >= KD_LFC_THRESHOLD or lfc <= (-1 * KD_LFC_THRESHOLD)) and adj_pvalue <= KD_ADJ_PVAL_THRESHOLD:
                    dict[ensg_id_kd] = True         
    return dict, dict_reverse

def getKey(item):
    return item[1]

def sample_wor(population, k):
    "Chooses k random elements (with out replacement) from a population"
    list_samples = random.sample(population, k)
    sample_indices = []
    for item in list_samples:
        index = population.index(item)
        if index not in sample_indices:
            sample_indices.append(index)
        else:
            index2 = population[index+1:].index(item)
            sample_indices.append(index2+index+1)
    #print(sample_indices)
    return sample_indices
    
def select_sample(population, k, type):
    len_list = []

    if type == 'nontarget':
        for index in population:
              len_list.append(nontarget_ids_lens[index])
        sorted_indices = sorted(range(len(len_list)), key=lambda k: len_list[k], reverse = True)
    elif type == 'target':
        for index in population:
              len_list.append(target_ids_lens[index])
        sorted_indices = sorted(range(len(len_list)), key=lambda k: len_list[k])    
    
    sample = []
    for i in range(k):
       sample.append(population[sorted_indices[i]])
    return sample
    

PIRANHA_PVAL_THRESHOLD = 0.05
PARAYZER_SCORE_THRESHOLD = 0.7
eCLIP_PVAL_THRESHOLD = 0.05
KD_LFC_THRESHOLD = 0
KD_ADJ_PVAL_THRESHOLD = 0.05


# the following file only contains CLIP peaks with an existing motif

RBP_OF_INTEREST = sys.argv[1]
dir = '/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/'

fclip = open (dir + 'gene_RBP_kmer_motifs_on_eCLIP_and_CLIPdb_peaks.txt')

#ffimo = open('FIMO_all3UTRs_allmotifs_0thorder_trans_based.txt')
ffimo = open('all_RBPs_top3kmer_matches.txt')

check_KDfile_HepG2 = '../DESeq2_results/HepG2/DESeq2_results_' + RBP_OF_INTEREST + '_smaller.txt'   
# fix this after running Deseq2
check_KDfile_K562 = '../DESeq2_results/K562/DESeq2_results_' + RBP_OF_INTEREST + '_smaller.txt'   

finput_len = open ("ensembl_grch37_gene_transcript_length_sorted.txt")

finput_mapping = open ("ensembl_biomart_hg19_proteincoding_gene_refseq_mapping.txt")
finput_mapping.readline()

list_piranha_rbps = ['DGCR8','PTBP1PTBP2','HNRNPU','HNRNPH','RBM10','HNRNPA2B1','PTBP1','HNRNPC','HNRNPF','HNRNPM','HNRNPA1','RBM47','TIA1','EIF4A3']
list_repressor = ["TIA1", "PUM2", "PUM1", "AUF1", "ZFP36", "KHSRP", "HNRNPC"]
list_activator = ["HNRNPL", "IGF2BP2", "IGF2BP3", "ELAVL1"]

fout_target = open ('./results/' + RBP_OF_INTEREST + '_target_genes.txt', 'w')

fout_nontarget = open ('./results/'  + RBP_OF_INTEREST + '_nontarget_genes.txt', 'w')


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
dict_count_motifs_in_signif_peaks = {}
dict_count_motifs_in_all_peaks = {}
#list_rbps = []
list_ensg_ids = []
dict_count_any_peak  = {}

count = 0
num_genes_with_signif_peak = 0
num_genes_with_any_peak = 0
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
        dict_count_motifs_in_signif_peaks[ensg_id] = 0
        dict_count_motifs_in_all_peaks[ensg_id] = 0
        if ensg_id == 'ENSG00000213927':
            print 'ENST id ', enst_id
    if RBP_OF_INTEREST in line:     
        dict_count_any_peak[ensg_id] = True
        #num_genes_with_any_peak += 1
        for i in range(1,len(peaks)):
            peaks_words = peaks[i].split(":")
            if '_' in peaks_words[0]:
                rbp_id = peaks_words[0].split('_')[0]
            else:
                rbp_id =  peaks_words[0]  
            if rbp_id == RBP_OF_INTEREST:
                   #print 'HERE'
                   flag = False #initially this CLIP peak is not considered to be true positive
                   clip_info = peaks_words[3]
                   num_motifs =  int(peaks_words[-1])  
                   dict_count_motifs_in_all_peaks[ensg_id] += num_motifs
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
                       dict_count_signif_peaks[ensg_id] += 1
                       #if num_motifs > 0:
                       #     print "HERE"
                       dict_count_motifs_in_signif_peaks[ensg_id] += num_motifs
 
#simple check  
for ensg_id in dict_count_motifs_in_signif_peaks:
   if dict_count_motifs_in_signif_peaks[ensg_id] > 0:
        num_genes_with_signif_peak += 1
   if dict_count_motifs_in_all_peaks[ensg_id] > 0:
        num_genes_with_any_peak += 1
   
print 'num_genes_with_signif_peak and any_peak ', num_genes_with_signif_peak, num_genes_with_any_peak
 
 
print 'read the CLIP file' 
   
dict_motif = {} 
for line in ffimo:   
   words = line.strip().split('\t')
   enst_id = words[0]
   ensg_id = dict_mapping_enst_to_ensg[enst_id]
   dict_motif[ensg_id] = 0
   for word in words[1:]:
       word_splitted = word.split(',')
       if word_splitted[0] == RBP_OF_INTEREST:
          num_sites = len(word_splitted) - 1
          dict_motif[ensg_id] = num_sites   
           
print 'read the motif file'          
     
K562_not_exist = False
dict_KD_K562 = {}
if  os.path.isfile(check_KDfile_K562):
	fKD_K562 = open(check_KDfile_K562)
	dict_KD_K562, dict_KD_K562_reverse = read_kd_result(fKD_K562)
	fKD_K562.close()
else:
    K562_not_exist = True 
	
HepG2_not_exist = False
dict_KD_HepG2 = {}	
if os.path.isfile(check_KDfile_HepG2):
	fKD_HepG2 = open(check_KDfile_HepG2)
	dict_KD_HepG2, dict_KD_HepG2_reverse  = read_kd_result(fKD_HepG2)
	fKD_HepG2.close()	
else:
    HepG2_not_exist = True 	

print 'read the KD files' 

# nontarget: no peak, no motif (whether it's in the peak or not) at all, signif LFC
# target: at least one signif peak, at least one motif in the signif peak, no signif LFC

select_oneortwo_targets = False
mult_factor_init = 1.2
target_ids = []
nontarget_ids = []

target_ids_lens = []
nontarget_ids_lens = []

count_nolen = 0
count_first_pass = 0
count_second_pass = 0
step = 'None'
for ensg_id in CLIP_file_all_ENSG_ids:
   if ensg_id in dict_length:
        
        motif_check_target = True
        motif_check_target = dict_count_motifs_in_all_peaks[ensg_id] > 0  
        
        motif_check_target = dict_count_motifs_in_signif_peaks[ensg_id] > 0    
        KD_check_target = True
        #KD_check_target = (ensg_id in dict_KD_K562 and dict_KD_K562[ensg_id]) or (ensg_id in dict_KD_HepG2 and dict_KD_HepG2[ensg_id])    
        if dict_count_signif_peaks[ensg_id] > 0 and motif_check_target and KD_check_target:
              count_first_pass += 1
              target_ids.append(ensg_id)
              target_ids_lens.append(dict_length[ensg_id])
              
        motif_check_nontarget = True      
        motif_check_nontarget = ensg_id not in dict_motif or dict_motif[ensg_id] == 0 
        KD_check_nontarget = True
        KD_check_nontarget = (K562_not_exist or (ensg_id in dict_KD_K562 and not dict_KD_K562[ensg_id])) or (HepG2_not_exist or (ensg_id in dict_KD_HepG2 and not dict_KD_HepG2[ensg_id]))        
        KD_check_nontarget = (K562_not_exist or (ensg_id in dict_KD_K562 and not dict_KD_K562[ensg_id])) and (HepG2_not_exist or (ensg_id in dict_KD_HepG2 and not dict_KD_HepG2[ensg_id]))        

        if dict_count_any_peak[ensg_id] == 0 and motif_check_nontarget and KD_check_nontarget:
               nontarget_ids.append(ensg_id)
               nontarget_ids_lens.append(dict_length[ensg_id])
   else:
       count_nolen += 1


print 'len target list ', len(target_ids_lens)
print 'len nontarget list ', len(nontarget_ids_lens)

binsize = 120
num_first_bins = (max_len // binsize) + 1
start = 0




#list_rbps = ["IGF2BP2", "IGF2BP3", "HNRNPL", "AUF1", "ZFP36", "KHSRP", "HNRNPC"]
          
        

list_targets_sorted = []

print RBP_OF_INTEREST, len(target_ids), len(nontarget_ids) 
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
		
print 'num_first_bins ', num_first_bins		
		
count_check = 0	
for bin in dict_bins_targets:
   count_check += len(dict_bins_targets[bin])
print 'count check ', count_check

	
dict_bins_nontargets = {}
for i in range (len(nontarget_ids_lens)):
	length = nontarget_ids_lens[i]
	bin = length // binsize
	num_nontargets[bin] += 1
	if bin not in dict_bins_nontargets:
		dict_bins_nontargets[bin] = [i]
	else:
		dict_bins_nontargets[bin].append(i)

num_targets_nozero = 0       
num_tarnozero_nontarzero = 0
num_both_nozero = 0
num_more_targets = 0
num_more_nontargets = 0


count_targets_nozero = 0
count_tarnozero_nontarzero = 0
count_both_nozero = 0
count_more_targets = 0
count_more_nontargets = 0 


ff = 0

gg = 0

selected_target_indices = []
selected_nontarget_indices = []
for i in range (start, num_first_bins):
	target_lens = []
	nontarget_lens = []
	if num_targets[i] != 0:
		num_targets_nozero += 1
		count_targets_nozero += num_targets[i]
	if num_targets[i] != 0 and num_nontargets[i] == 0:
		num_tarnozero_nontarzero += 1
		count_tarnozero_nontarzero += num_targets[i]
		if select_oneortwo_targets:
		   if num_targets[i] <= 2:
		       selected_target_indices.extend(dict_bins_targets[i])
	if num_targets[i] != 0 and num_nontargets[i] != 0:
		num_both_nozero += 1
		count_both_nozero += num_targets[i]
		if num_targets[i] <= num_nontargets[i]:
			mult_factor = 1
			num_more_nontargets += 1
			count_more_nontargets += num_targets[i]

			if num_targets[i] * 2 <= num_nontargets[i]:
			    gg += 1
			    mult_factor = mult_factor_init
			#sample = sample_wor(dict_bins_nontargets[i], int(len(dict_bins_targets[i]) * mult_factor))
			#sample = np.random.choice(dict_bins_nontargets[i], int(len(dict_bins_targets[i]) * mult_factor), replace=False)
			sample = select_sample(dict_bins_nontargets[i], int(len(dict_bins_targets[i]) * mult_factor), 'nontarget')
			global_indices = []
			selected_nontarget_indices.extend(sample)
			selected_target_indices.extend(dict_bins_targets[i])


		elif num_targets[i] > num_nontargets[i]:
			#if num_targets[i] - num_nontargets[i] <= 2:
			#	selected_target_indices.extend(dict_bins_targets[i])
			#	selected_nontarget_indices.extend(dict_bins_nontargets[i])
		 
			num_more_targets += 1
			count_more_targets += num_targets[i]
			#print 'targets more than nontargets', len(dict_bins_targets[i]), len(dict_bins_nontargets[i])
			#sample = sample_wor(dict_bins_targets[i] , len(dict_bins_nontargets[i]))
			sample = select_sample(dict_bins_targets[i], len(dict_bins_nontargets[i]), 'target')
			selected_target_indices.extend(sample)
			selected_nontarget_indices.extend(dict_bins_nontargets[i])



print 'num_first_bins ', num_more_nontargets, gg

selected_target_lens = []
selected_target_ids = []

selected_nontarget_lens = []
selected_nontarget_ids = []

for i in range(len(selected_target_indices)):
	target_id = target_ids[selected_target_indices[i]]
	selected_target_ids.append(target_id)
	selected_target_lens.append(target_ids_lens[selected_target_indices[i]])


for i in range(len(selected_nontarget_indices)):
	nontarget_id = nontarget_ids[selected_nontarget_indices[i]]
	selected_nontarget_ids.append(nontarget_id)
	selected_nontarget_lens.append(nontarget_ids_lens[selected_nontarget_indices[i]])


a, pvalue_len = mannwhitneyu(selected_target_lens, selected_nontarget_lens)
print 'p-value mannwhitneyu in terms of length ', RBP_OF_INTEREST , pvalue_len, np.median(selected_target_lens), np.median(selected_nontarget_lens)

#if count_cc != len(selected_target_ids):
	#print(count_aa,count_bb,count_cc, len(selected_target_ids))
print RBP_OF_INTEREST, len(selected_target_ids) , len(selected_nontarget_ids)
for item in selected_nontarget_ids:
	fout_nontarget.write(item + "\n")
fout_nontarget.close()

for item in selected_target_ids:
	fout_target.write(item + "\n")
fout_target.close()

