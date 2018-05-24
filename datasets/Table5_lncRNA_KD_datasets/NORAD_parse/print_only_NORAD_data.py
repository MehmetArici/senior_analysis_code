
fhinput = open('NORAD_ulitsky_Pum_NORAD_KD_OE.txt')

fhoutKD = open('../NORAD_Ulitsky_KD.txt', 'w')

fhoutOE = open('../NORAD_Ulitsky_OE.txt', 'w')

line = fhinput.readline()
print >>fhoutKD, 'ENSG_id\tLFC'
print >>fhoutOE, 'ENSG_id\tLFC'

dictKD_data = {}
dictOE_data = {}
for line in fhinput:
      words = line.split()
      ENSG_ids = words[0]
      print ENSG_ids
      ENSG_id_list = []
      if '_' in ENSG_ids:
          ENSG_id_list = ENSG_ids.split('_')
      else: 
          ENSG_id_list = [ENSG_ids]
          
      for ENSG_id in ENSG_id_list:
         if len(words) > 5:
              val = float(words[5].replace(',', '.'))
              if ENSG_id not in dictKD_data:
                   dictKD_data[ENSG_id] = val
              else:
                   dictKD_data[ENSG_id] = (dictKD_data[ENSG_id] + val) / 2      
         if len(words) > 6:
              val = float(words[6].replace(',', '.'))
              if ENSG_id not in dictOE_data:
                   dictOE_data[ENSG_id] = val
              else:
                   dictOE_data[ENSG_id] = (dictOE_data[ENSG_id] + val) / 2  
              
              
for ENSG_id in dictKD_data:
   print >>fhoutKD, ENSG_id + '\t' + str(dictKD_data[ENSG_id])

for ENSG_id in dictOE_data:
   print >>fhoutOE, ENSG_id + '\t' + str(dictOE_data[ENSG_id])
           
fhoutKD.close()
fhoutOE.close()
               
 