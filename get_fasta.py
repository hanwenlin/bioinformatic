
ofile = open('mir.fasta','w')
with open('Differentially_Expressed_tRF.txt','r') as f:
	header = f.readline()
	for i in f:
		lines = i.split('\t')
		ofile.write('>'+lines[0]+'\n')
		print(lines[1])
		ofile.write(lines[1]+'\n')
		
ofile.close()
