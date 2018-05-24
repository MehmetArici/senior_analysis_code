
fhinput = open('all_RBPs_top3kmer_matches_.txt')

fhout = open('all_RBPs_top3kmer_matches.txt', 'w')


for line in fhinput:
 if 'PTBP1' in line:
   words = line.split()
   ensg_id = words[0]
   forprint = ensg_id + '\t'
   for word in words[1:]:
      if word.startswith('PTBP1'):
          motifs = word.split(',')
          PTB_print = 'PTBP1,' 
          for i in range(1,len(motifs), 2):
             PTB_print += motifs[i] + ','
          PTB_print = PTB_print[:-1] 
          forprint += PTB_print + '\t'
      else:
          forprint += word + '\t'        
   forprint = forprint[:-1]
   print >>fhout, forprint
 else:
    print >>fhout, line,                 