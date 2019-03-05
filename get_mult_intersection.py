#!/usr/bin/env python

import xlrd

def read_excle():
	global setxl
	workbook = xlrd.open_workbook('Differentially Expressed circRNAs.xls')
    # print(workbook.sheet_names())
	sheet_names = workbook.sheet_names()
	#sheetset = set()
	for i in range(1,len(sheet_names),2):
		sheetset = set()
		sheet = workbook.sheet_by_index(i)
		for circname in sheet.col_values(0)[19:]:
			sheetset.add(circname)
			if i==1:
				setxl.add(circname)
		setxl = setxl&sheetset

setxl = set()
read_excle()
print(len(setxl),setxl)
circ ={}
for i in setxl:
	circ[i] = i
	#circ[i] =''

workbook = xlrd.open_workbook('Differentially Expressed circRNAs.xls')
sheet_names = workbook.sheet_names()
for i in range(1,len(sheet_names),2):
	sheet = workbook.sheet_by_index(i)
	for j in range(19,sheet.nrows):
		lines = sheet.row_values(j)
		if circ.get(lines[0]):
			circ[lines[0]] +='\t'
			if i <= len(sheet_names)-2:	
				circ[lines[0]] += '\t'.join([str(li) for li in lines[1:7]])
			else:
				circ[lines[0]] += '\t'.join([str(li) for li in lines[1:]])
		
	#	if circ.get(lines[0])=='' or circ.get(lines[0]):
			#circ[lines[0]] += '\t'.join([str(li) for li in lines])
		#	circ[lines[0]] +='\t'
			#print(circ[lines[0]])
						
with open('result.txt','w') as f:
	f.write(circ['CircRNAID']+'\n')
	circ.pop('CircRNAID')
	for value in circ.values():
		f.write(value+'\n')







#sheet1 = workbook.sheet_by_index(1)
#print(sheet1.row_values(19))
#'''
#	第一个sheet得到一个交集的文件，然后从第二个sheet开始依次打开这个交集文件
#'''
#with open('intersect.txt','w') as f:
	#with open('result.txt','w') as f1:
	#	for k in range(19,sheet1.nrows):
		#	lines = sheet1.row_values(k)
		#	if circ.get(lines[0]):
		#		f.write('\t'.join([str(li) for li in lines])+'\n')
			#	f1.write('\t'.join([str(li) for li in lines])+'\n')



#for i in range(3,len(sheet_names),2):
	#circss ={}
	#for j in range(19,sheet.nrows):
	#	lines = sheet.row_values(j)
	#	circss[lines[0]] = '\t'.join([str(li) for li in lines])
	#with open('result.txt','r') as f:
		#for linee in f:
			#linesss = linee.split()
			#if circss.get(linesss[0]):
			#	f.write(linee.strip()+'\t'+circss[linesss[0]]+'\n')

