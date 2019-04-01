#!/usr/bin/env python
#-*-coding:utf-8 -*-

with open('swissport.fasta','r') as f1,open('swissport_title.txt','w') as f2:
	for line in f1:
		if line.startswith('>'):
			f2.write(line)
			print(line.split('|')[1])
