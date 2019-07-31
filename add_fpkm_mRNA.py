#!/usr/bin/env python
#-*- coding:utf-8 -*-

import xlrd
import xlwt
from datetime import date,datetime

gene_fpkms = {}

def read_excle():
    workbook = xlrd.open_workbook('mRNA Expression Profiling.xlsx')
    sheet2_name = workbook.sheet_names()
    sheet2 = workbook.sheet_by_index(0)
    
    rowNum = sheet2.nrows
    
    for i in range(25,rowNum,1):
        gene_fpkm = '\t'.join([str(colitem.value) for colitem in sheet2.row(i)[6:18]])
        if i==25:
            gene_fpkms['gene'] = gene_fpkm
	    
        gene_fpkms[sheet2.row(i)[0].value.strip()] = gene_fpkm

#read_excle()

def set_style(name, size, bold=False):
    style = xlwt.XFStyle()  # 初始化样式

    font = xlwt.Font()  # 为样式创建字体
    font.name = name  # 'Times New Roman'
   # font.bold = bold
    # f.underline= Font.UNDERLINE_DOUBLE
    font.color_index = 4
    #font.height = height
    font.size = size

    # borders= xlwt.Borders()
    # borders.left= 6
    # borders.right= 6
    # borders.top= 6
    # borders.bottom= 6

    style.font = font
    # style.borders = borders

    return style
    
    
    
 
def write_excel():
   # f = xlwt.Workbook()  # 创建工作簿
    f2 = xlrd.open_workbook(r'Differentially Expressed mRNAs.xlsx')
    read_excle()
   
    
    sheet_namess = f2.sheet_names()
    #sheet1,sheet2 = f2.sheets()
    
    
    for sheetname in sheet_namess:
        sheet1 = f2.sheet_by_name(sheetname)
        rowNum = sheet1.nrows
       # print(sheet1)
        filename = sheetname+'.txt'
        with open(filename,'w') as f1:
            for i in range(27,rowNum,1):
                #f1.write('\t'.join([str(colitem.value) for colitem in sheet1.row(i)[0:12]]))
                values = gene_fpkms[sheet1.cell(i,0).value]
                if values.startswith('0.0'):
                    if values.endswith('0.0'):
                        va1 = values.replace('0.0','0')
                    else:
                        va1 = values.replace('0.0','0',1)
                else:
                    va1 = values
                f1.write(va1+'\n')
                #f1.write('\t'.join([str(colitem.value) for colitem in sheet1.row(i)[12:]])+'\n')
				
    
    

    
write_excel()
   

    
    
  
