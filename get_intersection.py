#!/usr/bin/env python
from sys import argv

fname1 = argv[1]   #circ_down_HPRT_VS_B.txt
fname2 = argv[2]   #circ_down_HPRT-HM_VS_B-HM.txt

f1 = open(fname1,'r')
f2 =open(fname2,'r')
ona1 = fname1.replace('.txt','')
ona2 = fname2.replace('mRNA_up_','')
seqname ={}
for i in f1:
    seqname[i.split('\t')[1]] = '\t'.join(i.split('\t')[4:12])

oname = ona1+'_intersect_'+ona2
f =open(oname,'w')

for j in f2:
    part = j.strip().split('\t')
    if seqname.get(part[1]):
        f.write('\t'.join(part[0:12])+'\t'+seqname[part[1]]+'\t'+'\t'.join(part[13:]))
        f.write('\n')
