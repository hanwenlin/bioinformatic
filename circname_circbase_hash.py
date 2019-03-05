#!/usr/bin/env python

import xlrd
circ = {}

def read_excle():
	workbook = xlrd.open_workbook("All Comparisons.xlsx")
	sheet_names = workbook.sheet_names()
	for i in range(0,len(sheet_names)):
		sheet = workbook.sheet_by_index(i)
		#print(sheet.row_values(19,20))
		nrows = sheet.nrows
		ncols = sheet.ncols
		for j in range(19,nrows):
			circname = sheet.row_values(j)[0]
			circbase = sheet.row_values(j)[23]
			#print(circname,circbase)
			circ[circname] = circbase
	#	for circname in sheet.col_values(0)[19:]:
	#		sheetset.add(circname)
		#	if i==1:
		#		setxl.add(circname)
		#setxl = setxl&sheetset
		
read_excle()

with open('circname2circbase.txt','w') as f:
	for k,v in circ.items():
		if v:
			f.write(k+'\t'+v+'\n')
		else:
			f.write(k+'\t'+k+'\n')
