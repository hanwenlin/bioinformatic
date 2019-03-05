#!/usr/bin/env python
import os


alldir = os.listdir(r'C:\mirTarget_results\2018111419A_CircRNAs_targetscan')
print(alldir[1:5])
mirseq = {}
with open('mature_hsa.fa','r') as f:
	for i in f:
		if i.startswith('>'):
			mirname = i.strip('>').strip()
			mirseq[mirname] = ''
		else:
			mirseq[mirname] += i

remain_fasta = open('remain_fasta.fa','w')			
for k in mirseq.keys():
	if k not in alldir:
		remain_fasta.write('>'+k+'\n')
		remain_fasta.write(mirseq[k])
remain_fasta.close()
