#!/usr/bin/env python
#-*- coding:utf8 -*-

codon = {'GGU':'G','GGC':'G','GGA':'G','GGG':'G','GCU':'A','GCC':'A','GCA':'A','GCG':'A',
		  'GUU':'V','GUC':'V','GUA':'V','GUC':'V','CUU':'L','CUC':'L','CUA':'L','CUC':'L','UUA':'L','UUG':'L',
		  'AUU':'I','AUC':'I','AUA':'I','CCU':'P','CCA':'P','CCG':'P','CCC':'P','UUU':'F','UUC':'F','UAU':'Y','UAC':'Y',
		  'UGG':'W','UCU':'S','UCA':'S','UCC':'S','UCG':'S','AGU':'S','AGC':'S','ACU':'T','ACC':'T','ACG':'T','ACA':'T',
		  'UGU':'C','UCC':'C','AUG':'M','AAU':'N','AAC':'N','CAA':'Q','CAG':'Q','GAU':'D','GAC':'D','GAA':'E','GAG':'E',
		  'AAA':'K','AAG':'K','CGU':'R','CGC':'R','CGG':'R','CGA':'R','AGA':'R','AGG':'R','CAU':'H','CAC':'H',
		  'UAG':'STOP','UGA':'STOP','UAA':'STOP'
		}
	


with open('A0666-RNA.fasta','r') as f:
	seq = ''
	for line in f:
		if not line.startswith('>'):
			seq += line.strip()
	for j in range(3):
		amon = ''	
		for i in range(j,len(seq)-2,3):
			codes = seq[i:i+3]
			if codon.get(codes):
				if codon.get(codes) == 'STOP':
					amon += '*'
				amon += codon[codes]
			else:
				amon +='-'
		print('Frame %d:\n' %(j+1))
		print(amon,'\n')
		

