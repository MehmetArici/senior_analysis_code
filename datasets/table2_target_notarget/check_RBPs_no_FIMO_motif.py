fRBPs = open('../Table0_list_of_RBPs/SpongeWebsite_RBP_list.txt')

ffimo = open('FIMO_all3UTRs_allmotifs_0thorder_trans_based.txt')

fseqs = open('')


dict_RBP_total_motif_count = {}
dict_RBP_total_gene_count = {}
RBPs = []
for line in fRBPs:
  RBP = line.split()[0]
  RBPs.append(RBP)
  dict_RBP_total_motif_count[RBP] = 0
  dict_RBP_total_gene_count[RBP] = []


for line in ffimo:   
   words = line.split('\t')
   enst_id = words[0]
   for word in words[1:]:
       word_splitted = word.split(',')
       RBP = word_splitted[0]
       if RBP in RBPs:
           dict_RBP_total_motif_count[RBP] += 1
           if enst_id not in dict_RBP_total_gene_count[RBP]:
                dict_RBP_total_gene_count[RBP].append(enst_id)
                
for RBP in RBPs:
   print RBP, dict_RBP_total_motif_count[RBP], len(dict_RBP_total_gene_count[RBP])             
