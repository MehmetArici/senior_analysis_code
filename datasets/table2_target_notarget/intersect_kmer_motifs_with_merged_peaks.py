import sys


dir = '/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/'

#ffimo_matches = open(dir + '/FIMO_files/FIMO_all3UTRs_allmotifs_0thorder.txt')
#fpeaks = open (dir  + 'RBP_CLIP_peaks_on_genes_eCLIP_and_CLIPdb_merged.txt')
#fhout = open(dir + 'gene_RBP_FIMO_motifs_on_eCLIP_and_CLIPdb_peaks.txt', 'w')

# CHANGED NOW ON 20 MARCH 2018 from FIMO to top3 kmers
fkmer_matches = open('all_RBPs_top3kmer_matches.txt')
fpeaks = open (dir  + 'RBP_CLIP_peaks_on_genes_eCLIP_and_CLIPdb_merged.txt')
fhout = open(dir + 'gene_RBP_kmer_motifs_on_eCLIP_and_CLIPdb_peaks.txt', 'w')


fkmer_matches.readline()
dict_trans_motifs = {}
for line in fkmer_matches: 
   words = line.split()
   enst_id = words[0]
   if enst_id not in dict_trans_motifs:
       dict_trans_motifs[enst_id] = {}
   for word in words[1:]:
       word_splitted = word.split(',')
       RBP_id = word_splitted[0]
       if RBP_id not in dict_trans_motifs[enst_id]:
        dict_trans_motifs[enst_id][RBP_id] = []
        for motif in word_splitted[1:]:
              motif_splitted = motif.split(':')
              spos = int(motif_splitted[0])
              epos = int(motif_splitted[1])
              dict_trans_motifs[enst_id][RBP_id].append([spos, epos])
   
#print dict_trans_motifs

#count_trans = 0
#holds the peaks on a transcript
dict_trans_peaks = {}
for line in fpeaks:
   words = line.split()
   enst_id = words[0]
   forprint =  enst_id + '\t' 
   peaks = words[1:]
   for peak in peaks:
       peak_words = peak.split(':')
       if 'eCLIP' in peak:
            RBP_id = peak_words[0].split('_')[0]
       else:
            RBP_id = peak_words[0]
       #print RBP_id     
       peak_spos = int(peak_words[1])
       peak_epos = int(peak_words[2])
       count = 0
       if enst_id in dict_trans_motifs and RBP_id in dict_trans_motifs[enst_id]:
             #print 'inside if '
             for motif_match in dict_trans_motifs[enst_id][RBP_id]:
                 motif_spos = motif_match[0]
                 motif_epos = motif_match[1]
                 #print peak_spos, motif_spos, motif_epos, peak_epos
                 if  peak_spos <= motif_spos and motif_epos <= peak_epos: # motif is inside the peak
                     count += 1
                     #if count >= 1:
                     #    print 'count is ', count
                 elif (motif_spos <= peak_spos  and peak_spos<= motif_epos) or (motif_spos <= peak_epos and peak_epos<=motif_epos):
                     #print motif_spos, motif_epos, peak_spos, peak_epos, 'INTERESTING OVERLAP'
                     count += 1
#             if count >= 1:
#                  count_trans += 1
       peak = peak + ':' + str(count) # we add it even if it's 0 
       forprint =  forprint + peak + '\t'    
   print >>fhout, forprint              
     
#print count_trans                
             






