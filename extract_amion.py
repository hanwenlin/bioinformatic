#!/usr/bin/env python
#-*- coding:utf-8 -*-

aa_codes ={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F',
			'GLY':'G','HIS':'H','LYS':'K','ILE':'I','LEU':'L',
			'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
			'SER':'S','THR':'T','VAL':'V','TYR':'Y','TRP':'W',
		  }
		  
seq =''
with open('1td.txt','r') as f:
	for line in f:
		if line.startswith('SEQRES'):
			for amion in line.strip().split()[4:]:
				if aa_codes.get(amion):
					seq += aa_codes[amion]
				else:
					print(amion)
print(seq)
