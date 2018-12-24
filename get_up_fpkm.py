#!/usr/bin/env python

symbols = {}
with open('mrna_up.txt','r') as f:
	for i in f:
		symbols[i.strip()] = i
		
ou = open('heatmap_SW1116+_vs_SW1116-.mRNA.txt.copy.txt','w')
with open('heatmap_SW1116+_vs_SW1116-.mRNA.txt','r') as f1:
	for i in f1:
		lines = i.split('\t')
		if symbols.get(lines[0]):
			ou.write(i)
