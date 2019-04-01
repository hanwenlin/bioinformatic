#!/usr/bin/env python
#-*- coding:utf8 -*-

seq = ''
with open('sequence.txt','r') as f:
	for i in f:
		seq += i.strip()
		
bases = 'ATCG'
seq_len = len(seq)
freqs = []

with open('statistic_bases.txt','w') as f:
	for base in bases:
		freq = seq.count(base)/float(seq_len)
		freqs.append(freq)
		f.write('%s is %.3f \n' %(base,freq))
		
	f.write('the maximum of freqency is %.3f' %(max(freqs)))
	

gc = (seq.count('G')+seq.count('C'))/float(seq_len)
print(gc)
