# This code is run from /Volumes/normal/UAU/Projects/lncrnas_transfer/CLIP_datasets/codes/

fhlist = open('../list_ENSTids_longest_transcripts_unique.txt') # these are the transcripts that I need

fhseqs = open('../ensembl_grch37_proteincoding_3utr_sequences_single.txt')

fhcoords = open('../ensembl_gencode_v25lift37_proteincoding_coordinates.txt')

fhout = open('../ensembl_grch37_protcoding_longest_enst_3utr_seqs_sponge.txt', 'w')


transcripts_needed = []
for line in fhlist:
    words = line.split()
    transcripts_needed.append(words[0])
    
dict_coords = {}
for line in fhcoords:
   words = line.strip().split()
   enst_id = words[0]
   chr = words[1]
   strand = words[2]
   exon = words[3]
   
   
   if enst_id in dict_coords:
      dict_coords[enst_id][2].append(exon)
   else:
      dict_coords[enst_id] = [chr, strand, [exon]]
      
      
need_print = False
for line in fhseqs:
   if '>' in line:
      words = line.strip().split('|')
      enst_id = words[2]
      need_print = False
      if enst_id in transcripts_needed:
         need_print = True
         ensg_id = words[0][1:]
         gene_symbol = words[1]
         coords = dict_coords[enst_id][2]
         coords_print = ''
         for coord in coords:
            coords_print += coord + ','
         coords_print = coords_print[:-1]   
         forprint = '>' + enst_id + '\t' + ensg_id + '\t' + gene_symbol + '\t' + dict_coords[enst_id][0] + '\t' + dict_coords[enst_id][1] + '\t' + coords_print
         print >>fhout, forprint
   else:
       if need_print:
           print >>fhout, line.strip()
    
    
    
    
            
          