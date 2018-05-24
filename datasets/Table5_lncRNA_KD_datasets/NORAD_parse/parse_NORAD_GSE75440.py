finput = open ("GSE75440_NORADedgeR.txt")
header = finput.readline()

fhout = open("NORAD_Mendell_KD.txt", 'w')
print >>fhout, 'ENSG_id\tLFC\tgene_symbol'
for line in finput:
    words = line.strip("\n").split("\t")
    ensg_id = words[0]
    lfc = words[1]
    fhout.write(ensg_id + "\t" + lfc + '\t' + words[-1] + "\n")

fhout.close()
    