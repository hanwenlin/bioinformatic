#!/usr/bin/env python
import glob
circs ={}
with open('circname2circbase.txt') as f:
	for line in f:
		lines = line.strip().split()
		circs[lines[0]] = lines[1]
		
		
def circname2circbase(filename):
	circ ={}
	filedata = ''
	with open(filename,'r') as f1:
		for line in f1:
			lines = line.strip().split('\t')
			if circs.get(lines[0]):
				lines[0] = circs[lines[0]]
			eachline = '\t'.join(lines) +'\n'
			filedata +=eachline
	with open(filename,'w') as f2:
		f2.write(filedata)


for filename in glob.glob("heatmap.*txt"):
	circname2circbase(filename)
