#!/usr/bin/env python2.3
#
#
#  Created by Hilal Kosucu on 30/03/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#


import sys
import math
import itertools
from operator import itemgetter


def generate_all_possible_kmers(width):
        kmers = []
        for elem in itertools.product('ACGU', repeat=width):
                kmers.append(''.join(elem))
        return kmers

def list_to_str(lst):
    ''' Given a list, return the string of that list with tab separators 
    '''
    return reduce( (lambda s, f: s + '\t' + str(round(f,3))), lst, '')

def convert(base):
	if(base == 'A'):
		#print "yes"
		return 0
	if(base == 'C'):
		#print "yes"
		return 1
	if(base == 'G'):
		#print "yes"
		return 2
	if(base == 'U'):
		#print "yes"
		return 3	
	if(base == 'T'):
		#print "yes"
		return 3

def logistic(x):
	return 1 / (1+math.exp(-1 *x))

def generate_all_possible_kmers(width):
        kmers = []
        for elem in itertools.product('ACGU', repeat=width):
                kmers.append(''.join(elem))
        return kmers
	
class kmerClass:
        def __init__(self, seqid, pos, seqscore, strscore, totalscore, seq):
			self.seqid = seqid
			self.pos = pos
			self.seqscore = seqscore
			self.strscore = strscore
			self.totalscore = totalscore 
			self.seq = seq
        def __repr__(self):
			return repr((self.seqid, self.pos, self.score))




def return_all_kmers(motif_width):
   if(motif_width == 4):
       return all_kmers_4
   elif(motif_width == 5):
       return all_kmers_5
   elif(motif_width == 6):
       return all_kmers_6
   elif(motif_width == 7):
       return all_kmers_7  
   elif(motif_width == 8):
       return generate_all_possible_kmers(8)             
   elif(motif_width == 9):
       return generate_all_possible_kmers(9)

fhparams = open('/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/FIMO_files/all_motifs_for_FIMO.txt') #parameters found b

fhspongeRBPs = open('../Table0_list_of_RBPs/SpongeWebsite_RBP_list.txt')

fseqs = open('/Volumes/normal/UAU/Projects/lncrnas_transfer/lncRNAs_Human/CLIP_datasets/ensembl_grch37_protcoding_longest_enst_3utr_seqs_sponge.txt')

fhout = open('all_RBPs_top3kmer_matches.txt', 'w')

dict_seqs = {}
for line in fseqs:
    if '>' in line:
       words = line.split()
       id = words[0][1:]
    else:
       seq = line.strip()
       dict_seqs[id] = seq   
   
all_kmers_4 = generate_all_possible_kmers(4)
all_kmers_5 = generate_all_possible_kmers(5)
all_kmers_6 = generate_all_possible_kmers(6)
all_kmers_7 = generate_all_possible_kmers(7)




sponge_RBP_list = []
for line in fhspongeRBPs:
  sponge_RBP_list.append(line.split()[0])

#sponge_RBP_list.append('HuR')

dict_topkmers = {}
line = fhparams.readline()
dict_exist = {}
min_mw = 20
max_mw = 0
while line != '':
   if line.startswith('MOTIF'):
       words = line.split()
       entry2 = words[1]
       entry3 = words[2]
       if entry2.startswith('RNCMPT'):
           RBP_name = entry3
       else:
           RBP_name = entry2 
       line = fhparams.readline()
       line = fhparams.readline()
       motif_width = int(line.split()[5])
       if motif_width < min_mw:
           min_mw = motif_width
       if motif_width > max_mw:
           max_mw = motif_width    
       if RBP_name in dict_exist:
          dict_exist[RBP_name].append((entry2, entry3,  motif_width) )
       else:
          dict_exist[RBP_name] = [(entry2, entry3,  motif_width) ]
          
       params = []
       for i in range(motif_width):
          line = fhparams.readline()
          words = [float(x) for x in line.split()]
          params.append(words[0:len(words)])	
       if RBP_name in sponge_RBP_list:   
         #print params
         all_kmers = return_all_kmers(motif_width)
         print motif_width
         kmer_scores = {}
         for kmer in all_kmers:
           kmer_score = 1 
           for pos in range(0, motif_width):
                kmer_score *= params[int(pos)][int(convert(kmer[pos]))]
           if kmer_score > 0:     
               kmer_scores[kmer] = kmer_score
         sorted_kmer_seq_objects = sorted(kmer_scores, key=kmer_scores.get, reverse=True) #alternative sorting method
         print sorted_kmer_seq_objects[0:3]
         topKmers = []
         numtopkmer = 3
         for i in range(0, min(numtopkmer, len(kmer_scores))):
          kmerSeq = sorted_kmer_seq_objects[i]
          topKmers.append(kmerSeq)
         if RBP_name not in dict_topkmers:
            dict_topkmers[RBP_name] = topKmers
         else:
            dict_topkmers[RBP_name] = list(set(dict_topkmers[RBP_name]).union(set(topKmers)))     
         print RBP_name, dict_topkmers[RBP_name]   

   line = fhparams.readline()      
   
count = 0 
print 'min max ', min_mw, max_mw
for id in dict_seqs:
    seq = dict_seqs[id]
    forprint = id + '\t'
    dict_matches = {}
    if count % 100 == 1:
       print count
    for motif_width in range(min_mw, max_mw + 1):
      for pos in range(0, len(seq) - motif_width + 1):	
        kmer = seq[pos:pos+motif_width]
        for RBP in sponge_RBP_list:
               if kmer in dict_topkmers[RBP]:
                   if RBP in dict_matches:
                       dict_matches[RBP].append(str(pos) + ':' + str(pos + motif_width))
                   else:
                       dict_matches[RBP] =  [str(pos) + ':' + str(pos + motif_width)]
                       
    for RBP in dict_matches:
        forprint += RBP + ',' + ','.join(dict_matches[RBP]) + '\t'                    
    forprint = forprint[:-1] # get rid of the last tab  
    print >>fhout, forprint           
    count += 1	  
		   

        
'''
count_dup = 0
for RBP in RBP_list:
   if len(dict_exist[RBP]) > 1:
        print RBP, dict_exist[RBP] 
        count_dup += 1
        
print count_dup, len(dict_exist)
'''        
'''
       params = []
       for i in range(motif_width):
          line = fhparams.readline()
          words = [float(x) for x in line.split()]
          params.append(words[1:len(words)])	
       
       kmer_scores = {}
       for kmer in all_kmers:
	       kmer_score = 1 
	       for pos in range(0, mw):
		      kmer_score *= params[int(pos)][int(convert(kmer[pos]))]
           kmer_scores[kmer] = kmer_score
       sorted_kmer_seq_objects = sorted(kmer_scores, key=kmer_scores.get, reverse=True) #alternative sorting method
       topKmers = []
       numtopkmer = 10;
       for i in range(0, numtopkmer):
          kmerSeq = sorted_kmer_seq_objects[i]
          topKmers.append(kmerSeq)
       dict_topkmers[]
'''
       
       
'''
seqs = []
gene_ids = []
for i in range(0, len(fileseqs)):
	words = fileseqs[i].split()
	seq = words[1]
	seqs.append(seq)
	gene_ids.append(words[0])

for i in range(0, len(seqs)):
	seq = seqs[i]
	id = gene_ids[i]
	indices = []
	occur_indices = []
	for pos in range(0, len(seq) - mw + 1):	
		kmer = seq[pos:pos+mw]
		if kmer in topKmers:
			occur_index = topKmers.index(kmer)
			indices.append(pos)
			occur_indices.append(occur_index)
			
	if len(indices) > 0:
		print >>fout, rbp_name, id,
		for j in range(0, len(indices)):
			print >>fout, '%d:%d'%(indices[j],occur_indices[j]),
		print >>fout, '\n',
'''