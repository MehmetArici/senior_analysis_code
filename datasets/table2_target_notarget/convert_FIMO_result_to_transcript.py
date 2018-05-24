
dir = '/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/'


fhinput = open(dir + '/FIMO_files/FIMO_all3UTRs_allmotifs_0thorder.txt')

fhout = open(dir + '/FIMO_files/FIMO_all3UTRs_allmotifs_0thorder_trans_based.txt', 'w')


fhinput.readline()
dict_trans_motifs = {}
for line in fhinput: 
   words = line.split()
   if words[1] == 'RBPmap' or words[1] == 'RBNS':
       RBP_id = words[0]
   else:
       RBP_id = words[1]

   enst_id = words[2]
   spos = words[3]
   epos = words[4]
   pval = float(words[7])
   if enst_id not in dict_trans_motifs:
       dict_trans_motifs[enst_id] = {}
   if RBP_id not in dict_trans_motifs[enst_id]:
        dict_trans_motifs[enst_id][RBP_id] = []
   dict_trans_motifs[enst_id][RBP_id].append([spos, epos])
   


for enst_id in dict_trans_motifs:
   forprint = enst_id + '\t'
   if len(dict_trans_motifs[enst_id]) > 0:
         RBPs =  dict_trans_motifs[enst_id]
         for RBP in RBPs:
               sites = dict_trans_motifs[enst_id][RBP]
               sites_print = RBP + ','
               for site in sites:
                    sites_print = sites_print + site[0] + ':' + site[1] + ','
               
               sites_print = sites_print[:-1]
               
               forprint = forprint + sites_print + '\t'
   print >>fhout, forprint
