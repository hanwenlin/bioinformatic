# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xlrd
import time

data = xlrd.open_workbook(r'E:/software/bioinformatics/lnc_go_pathway/LncRNA Targets.xlsx')
print('rithg')
time.sleep(4)

sheetnames = data.sheet_names()
for sheetid in range(len(sheetnames)):
    table = data.sheet_by_name(sheetnames[sheetid])
    #list1 = table.readline()
    nrows = table.nrows
    fgo = open(r'E:/software/bioinformatics/lnc_go_pathway/go/go.LncRNA.'+sheetnames[sheetid]+'.txt','w')
    fpw = open(r'E:/software/bioinformatics/lnc_go_pathway/pathway/pw.LncRNA.'+sheetnames[sheetid]+'.txt','w')
    
    for i in range(2,nrows):
        gene = table.cell_value(rowx=i, colx=17)
      
        if gene != '':
            fgo.write(gene+'\n')
            
            if 'up' in sheetnames[sheetid]:
                fpw.write(gene+'\torange' +'\n')
            else:
                fpw.write(gene+'\tyellow' +'\n')
    fgo.close()
    fpw.close()        
