# -*- coding: utf-8 -*-
import sys
reload(sys)
import time
sys.setdefaultencoding('utf8')
"""
This is used to format mRNA HTML Report (human version 3.0)
Two files/folders are needed:
1, config.txt file -- contains all the parameters used
2, mRNA Report folder (Just copy from Z disk)
"""
import ConfigParser
import codecs
import xlsxwriter
import xlrd
import math
import matplotlib as mpl
mpl.use('Agg')
mpl.use('Agg')
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt
from scipy import stats
from pylab import *
import os, glob, numpy,xlwt,shutil,re,ctypes,collections
from pyExcelerator import *
from rpy2.robjects import r
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats_rpy2 = importr('stats')



#############<<<<<<<<<<<<<<<Read config file and return sth
global_now_location=os.getcwd()
config_files = glob.glob("*config.txt")
config_file = config_files[0]
config_file=os.path.join(global_now_location,config_file)
print global_now_location
print config_file
cf = ConfigParser.ConfigParser()
cf.readfp(codecs.open(config_file,"r","utf_8"))
print cf.sections()
try:
	all_de_folder=os.mkdir("all_de_folder")
	all_de_loc=global_now_location+"/all_de_folder"
except:
	shutil.rmtree("all_de_folder")
	#os.removedirs("all_de_folder")
	time.sleep(1)
	all_de_folder=os.mkdir("all_de_folder")
	all_de_loc=global_now_location+"/all_de_folder"

# user_information
user_information_name           = cf.get("user_infomation", "name")
user_information_project_number = cf.get("user_infomation", "projcet_number")
user_information_department     = cf.get("user_infomation", "department")
user_information_species        = cf.get("user_infomation", "species")
user_information_sample_number  = cf.get("user_infomation", "sample_number")
#  venn_comparisons = cf.options("venn_fig_excel")
#  GeneSpring

# XoutY_X                         = cf.get('GeneSpring', 'XoutY_X')
# intentity_threshold             = cf.get('GeneSpring', 'threshold')
XoutY_X=str(1)
intentity_threshold=str(1)

#######################
heatmap_loc=global_now_location+"/heatmap"
try:
	os.mkdir("heatmap")
except:
	shutil.rmtree("heatmap")
	#os.removedirs("heatmap")
	time.sleep(1)
	os.mkdir("heatmap")
heatmap_code	=	cf.get("heatmap","heatmap_operation_code").strip()
########################################
mRNA_info_addline = -1 # increase this, if not enough space:rat:-1, mouse:-1，human: -1
cell_width=4000 # 3333 = 1" (one inch).
outof=XoutY_X+'outof'+user_information_sample_number ###adjust it
intentity_threshold=int(intentity_threshold)
intensity = intentity_threshold
#genespringversion	= cf.get('GeneSpring', 'GeneSpring_ver')
genespringversion	= '12.0'
lg="log2"
genesymbl="GeneName"
if user_information_species.lower() == 'human':
	source_description='circBase, Guojunjie2014,...'
elif user_information_species.lower() == 'mouse':
	source_description='circBase, Guojunjie2014,...'
elif user_information_species.lower() == 'rat':
	source_description='circBase, Guojunjie2014,...'
else:
	source_description=''

controltype="xxxx"
flagsuffix=".txt:gFEFlags"
#first_row_name='ProbeName'#adjust
venn_margin=.2 #adjust
collist=["red","green","blue","yellow","purple","pink","yellow","orange","brown"]
cell_width=4000 # 3333 = 1" (one inch).
############# CONSTANTS
# column_index 
column_index = {} # GZ (in most(99.99%) situations, this is enough, if not, please add in the next-next line)
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
#colalphabet = list(alphabet)+ ['A'+i for i in alphabet]+['B'+i for i in alphabet]+['C'+i for i in alphabet]+['D'+i for i in alphabet]+['E'+i for i in alphabet]+['F'+i for i in alphabet]+['G'+i for i in alphabet]+['H'+i for i in alphabet]+['I'+i for i in alphabet]+['J'+i for i in alphabet]+['K'+i for i in alphabet]
colalphabet = list(alphabet)
for i in range(1,len(alphabet)+1):
	for j in range(1,len(alphabet)+1):
		temp=alphabet[i-1]+alphabet[j-1]
		colalphabet=colalphabet + [temp]
for i in range(1, len(colalphabet)+1):
    column_index[i] = colalphabet[i-1]
# COLOR CONSTANTS
############warning colors <<<<<<<<<<<<<<<<<<<<<<<<<<<<
STD_INPUT_HANDLE = -10
STD_OUTPUT_HANDLE= -11
STD_ERROR_HANDLE = -12

FOREGROUND_WHITE = 0xf
FOREGROUND_BLACK = 0x0
FOREGROUND_BLUE = 0x01 # text color contains blue.
FOREGROUND_GREEN= 0x02 # text color contains green.
FOREGROUND_RED = 0x04 # text color contains red.

FOREGROUND_INTENSITY = 0x08 # text color is intensified.
BACKGROUND_BLUE = 0x10 # background color contains blue.
BACKGROUND_GREEN= 0x20 # background color contains green.
BACKGROUND_RED = 0x40 # background color contains red.
BACKGROUND_INTENSITY = 0x80 # background color is intensified.
############## colors
#std_out_handle = ctypes.windll.kernel32.GetStdHandle(STD_OUTPUT_HANDLE)
#def set_cmd_text_color(color, handle=std_out_handle):
#    """(color) -> BOOL
#
#    #Example: set_cmd_text_color(FOREGROUND_GREEN | FOREGROUND_INTENSITY)
#    #"""
#    #bool = ctypes.windll.kernel32.SetConsoleTextAttribute(handle, color)
#    #return bool

#def resetColor():
#    set_cmd_text_color(FOREGROUND_WHITE | FOREGROUND_INTENSITY)


def printError(mess):
    "print wrong message"
    #set_cmd_text_color(FOREGROUND_RED | FOREGROUND_INTENSITY)
    print(mess)
    #resetColor()

def printWarning(mess):
    "print warning message"
    #set_cmd_text_color(FOREGROUND_GREEN | FOREGROUND_INTENSITY)
    print(mess)
    #resetColor()
##############

title_background_color  =       'yellow'
title_FC_sample_color   =       'red'
title_log2_color        =       'blue'
fc_fdr_regulation_color =       'red'
up_color                =       'red'
down_color              =       'green'
raw_color               =       'yellow'
normalized_color        =       'blue'
annotation_color        =       'purple'

pathway_up_color        =       'orange'
pathway_down_color      =       'yellow'
############## COLOR CONSTANTS
def remove_space(string):
    #p=re.compile('\s+')
    #return re.sub(p,'',string)
    plist=re.split(",",string)
    templist=[]
    for i in plist:
    	templist.append(i.strip())
    p=",".join(templist)
    return p
##############sample-group
#sample-group,samples,groups
try:sp=cf.options("sample_group")
except:sp=[]
if sp:
	for i in range(1,len(sp)+1):
		sp_ele=remove_space(cf.get("sample_group",str(i))).split("\t")
		try:sample_group[sp_ele[0]]=sp_ele[1:]
		except:sample_group={sp_ele[0]:sp_ele[1:]}


##########################################

def formats(st):
    if '.' in str(st) or 'e' in str(st):
        try:
            return float(st)
        except:
            return st
    else:
        try:
            return int(st)
        except:
            return st
############## 数学函数 开始

def log2(a):
    return math.log(a)/math.log(2)


def unpaired_ttest(groupA, groupB):
	if groupA==groupB:
		return 1.0 # stats will calculate it as nan
	else:
		if (len(set(groupA))==len(set(groupB))==1) and (set(groupA)==set(groupB)):
			return 1.0
		else:
			return stats.ttest_ind(groupA, groupB)[1]

def paired_ttest(groupA, groupB):
	if groupA==groupB:
		return 1.0 # stats will calculate it as nan
	else:
		if (len(set(groupA))==len(set(groupB))==1) and (set(groupA)==set(groupB)):
			return 1.0
		else:
			return stats.ttest_rel(groupA, groupB)[1]


def means(content_list, numlist):
    """
    从原来list中根据索引提取数据，并计算均值
    input:
        content_list:   原list
        numlist:        要提取的list（数字）
    output:
        平均值，float
    """
    all = []
    for i in numlist:
        all.append(float(content_list[i]))
    return sum(all)/len(numlist)
    

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"): 
    """
    FDR修正，与R函数correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1])一致
    input:
        pvalues         p-value值list
        correction_type 默认是Benjamini-Hochberg
    output：
        correct_pvalues 修正后的p值list，即FDR
        
    """
    pvalues1=pvalues
    pvalues = numpy.array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = numpy.empty(n)
    if correction_type == "Bonferroni":
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue
    elif correction_type == "Benjamini-Hochberg":
        new_pvalues = stats_rpy2.p_adjust(FloatVector(pvalues1), method = 'BH')
    return new_pvalues
################ 数学函数 结束
"""
# STEP 1
# Export [Raw, Normalized,Flag,Annotation] LncRNA+mRNA GeneSpring file  -- all.txt
# Seperate into LncRNA and mRNA -- all.mRNA.txt all.LncRNA.txt


"""

def get_headline(file):
    f = open(file)
    for line in f:
        if not line.startswith('#'):
            needline = line
            break
    f.close()
    return needline
def get_wanted_txt(file):
	head=get_headline(file).strip("\n")
	headlist=head.split("\t")
	newlist=[]
	for i in headlist:
		if ".txt:gProcessedSignal" in i or ".txt(" in i:
			newname="["+i.replace(".txt:gProcessedSignal","]").replace(".txt","]")
			newlist.append(newname)
		else:newlist.append(i)
	newheadline="\t".join(newlist)
	ofile=open(file,"rU")
	otmpfile=open(file+".tmp","w")
	line="ini"
	count=0
	while line:
		line=ofile.readline().rstrip("\n")
		if line=="":
			break
		if not line.startswith("#"):
			count=count+1
			if count==1:
				otmpfile.write(newheadline+"\n")
			else:
				otmpfile.write(line+"\n")
	ofile.close()
	otmpfile.close()
	os.remove("all.mRNA")
	os.rename("all.mRNA.tmp","all.mRNA")

def check_name(rawnamelist):
	alllist=parse_2sample_comparison()+parse_2group_unpaired_comparison()+parse_2group_paired_comparison()+parse_sample_group_comparison()+parse_group_sample_comparison()
	for thelist in alllist:
		if len(thelist)==3:
			for e in thelist[:2]:
				if e in rawnamelist:
					pass
				else:
					liststr="~~~~"+" ".join(thelist)+"~~~~"
					#printError(e+' in config.txt is not in your sample names of raw data, please check comparison '+liststr+'and(or) [sample_group] names in the config.txt, there may be some spelling errors')
					os._exit(0)
		else:
			if isinstance(thelist[0],list) and isinstance(thelist[2],list):
				judgelist=thelist[0]+thelist[2]
			else:
				try:
					thelist[0].append(thelist[2])
					judgelist=thelist[0]
				except:
					thelist[2].append(thelist[0])
					judgelist=thelist[2]
			for e in judgelist:
				if e in rawnamelist:
					pass
				else:
					liststr="~~~~"+thelist[1]+' '+" ".join(thelist[3:])+"~~~~"
					#printError(e+' in config.txt is not in your sample names of raw data, please check comparison '+liststr+'and(or) [sample_group] names in the config.txt, there may be some spelling errors')
					os._exit(0)


def get_matrix_excel(file):
	headlist=get_headline(file).split("\t")
	myloc=[0]
	for i in headlist:
		if "US10450393" in i:
			#printError(u'flag中含有US10450393，请在genespring中按照样本名字改好flag，然后输出！')
			sys.exit()
	for i in headlist:
		if "(normalized)" in i:
			myloc.append(headlist.index(i))
	#print myloc,headlist
	ofile=open(file,"rU")
	line="ini"
	mymatrix=[]
	while line:
		line=ofile.readline().rstrip("\n")
		if line=="":
			break
		if not line.startswith("#"):
			sline=line.split("\t")
			mysline=[sline[j] for j in myloc]
			mymatrix.append(mysline)
	ofile.close()
	#print mymatrix[:2]
	firstrowtemp=[]
	for k in mymatrix[0]:
		if not "(normalized)" in k:
			firstrowtemp.append("ID_REF")
		else:
			firstrowtemp.append(k.replace("(normalized)","")[1:-1])
	mymatrix[0]=firstrowtemp
	#print mymatrix[:2]
	wb = xlwt.Workbook()
	ws0 = wb.add_sheet("Matrix")
	rownum=len(mymatrix)
	colnum=len(mymatrix[0])
	ws0.write(0, 0, mymatrix[0][0], style8)
	for i in range(1,colnum):
		ws0.write(0, i, mymatrix[0][i], style8)
	for i in range(1,rownum):
		for j in range(colnum):
			ws0.col(j).width = cell_width
			ws0.write(i, j, formats(unicode(mymatrix[i][j],"utf-8")), style6)
	wb.save('Matrix.xls')
	check_name(firstrowtemp)

def seperate_mRNA_LncRNA():
    """
    seperate all.txt into mRNA and LncRNA
    input: all.txt
    output: all.mRNA, all.LncRNA
    """
    headlist = get_headline('all.txt').strip('\n').split('\t')
    #GeneSymbol_index = headlist.index('GeneName')
    f = open('all.txt')
    fmRNA = open('all.mRNA','w')
    #print >>fmRNA, '\t'.join(headlist)    
    for line in f:
        if not line.startswith('#'):
            linelist = line.strip('\n').split('\t')
            print >>fmRNA, line,
    f.close()
    fmRNA.close()




def flag_x_out_y(file):
    """
    according to the flag number you choose, filter input files
    intput:     all.mRNA
    output:     all.mRNA.2outof3.txt 
    and heatmap_all.mRNA.2outof3.txt
    """
    f = open(file)
    start = 1+int(user_information_sample_number)*2
    end   = 1+int(user_information_sample_number)*3
    flag_indexs = range(start, end)
    line="ini"
    while line:
    	line=f.readline().rstrip("\n")
    	if not line.startswith("#"):
    		break
    head = line
    headerlist=re.split("\t",head) 
    newhead="\t".join(headerlist)
    fflag = open(file+'.'+XoutY_X+'outof'+ user_information_sample_number+'.txt','w') # all.mRNA.txt.3outof9
    fflag.write("%s\n"%newhead)
    for line in f:
        line_flags = []
        linelist = line.strip('\n').split('\t')
        for flag_index in flag_indexs:
            line_flags.append(linelist[flag_index])
        if ((line_flags.count('Detected')+line_flags.count('Not Detected'))!=0 and (line_flags.count('Detected')+line_flags.count('Not Detected')) >= int(XoutY_X)) or ((line_flags.count('P')+line_flags.count('M'))!=0 and (line_flags.count('P')+line_flags.count('M')) >= int(XoutY_X)):
            fflag.write("%s\n"%("\t".join(linelist)))
    f.close()
    fflag.close()
    nstart=1+int(user_information_sample_number)
    nend=int(user_information_sample_number)*2
    all_heat_matrix,data_len=all_txt(file+'.'+XoutY_X+'outof'+ user_information_sample_number+'.txt')
    matrix_len=len(all_heat_matrix)
    try:groupnumber=len(sample_group.values()[0])
    except:groupnumber=1
    for number in range(groupnumber):
        if groupnumber==1:
            fflag_for_heatmap=open('heatmap_'+file+'.'+XoutY_X+'outof'+ user_information_sample_number+'.txt','w')
        else:fflag_for_heatmap=open('heatmap_'+"group_type"+str(number)+file+'.'+XoutY_X+'outof'+ user_information_sample_number+'.txt','w')
        try:
            isinstance(sample_group,dict)
            rnstart=1
            rnend=int(user_information_sample_number)*2
            newheaderlist=[]
            for i in range(len(headerlist)):
                if i in range(rnstart,rnend+1):
                    samplename=headerlist[i].split("](")[0][1:]
                    typename=headerlist[i].split("](")[1]
                    newname="["+samplename+", "+sample_group[samplename][number]+"]("+typename
                    newheaderlist.append(newname)
                else:newheaderlist.append(headerlist[i])
            h_headerlist=[newheaderlist[0]]+newheaderlist[nstart:nend+1]
            print >>fflag_for_heatmap,"\t".join(h_headerlist)
        except:
            h_headerlist=[headerlist[0]]+headerlist[nstart:nend+1]
            print >>fflag_for_heatmap,"\t".join(h_headerlist)
        for j in range(1,matrix_len):
    	    templist=[all_heat_matrix[j][0]]+all_heat_matrix[j][nstart:nend+1]
    	    print >>fflag_for_heatmap,"\t".join(templist)
        fflag_for_heatmap.close()

def intensity_x_out_y(file):
    """
    according to the flag number you choose, filter input files
    intput:     all.mRNA
    output:     all.mRNA.2outof3.txt
    """
    #print file
    f = open(file)
    #print file
    start = 1+int(user_information_sample_number)*0
    end   = 1+int(user_information_sample_number)*1
    #print start
    nstart=1+int(user_information_sample_number)*1
    nend=1+int(user_information_sample_number)*2
    #print nstart
    flag_indexs = range(start, end)
    head = f.readline()
    headerlist=re.split("\t",head) 
    fflag = open(file+'.'+XoutY_X+'outof'+ user_information_sample_number+'.txt','w') # all.mRNA.txt.3outof9
    fflags = open('hheatmap_'+file+'.'+XoutY_X+'outof'+ user_information_sample_number+'.txt','w') # all.mRNA.txt.3outof9
    h_headerlist=[headerlist[0]]+headerlist[nstart:nend]
    #print "haha"
    print >>fflags,"\t".join(h_headerlist)
    totalrows=0
    filteredrows=0
    print >>fflag, head,
    for line in f:
        totalrows+=1
        linelist = line.strip('\n').split('\t')
        line_flags_count = 0
        #print "haha"
        for flag_index in flag_indexs:
            #print "haha"
            if(float(linelist[flag_index]) >= intentity_threshold):
				line_flags_count+=1
				
        if (line_flags_count >= int(XoutY_X)):
            filteredrows+=1
            print >>fflag, line,
            h_headerlist=[linelist[0]]+linelist[nstart:nend]
            print >>fflags,"\t".join(h_headerlist)
    f.close()
    #print "haha"
    sample_number=int(user_information_sample_number)
    print "%.f"%filteredrows +" rows out of "+"%.f"%totalrows+" rows are filtered with " + "%s"%XoutY_X+" out of "+"%.f"%sample_number+" samples intensity >"+"%s"%intensity+"!\n"
    fflag.close()
    fflags.close()




##################### 解析config文件，并进行比较
##################### 两样品比较
def parse_2sample_comparison():
    comparisons = cf.options("2_sample_comparison")  # pairs
    A_vs_B_FCs = []
    if len(comparisons) % 2 == 0:        
        for i in range(len(comparisons)/2):
            AvsB            = remove_space(cf.get('2_sample_comparison','2_sample_comparison_'+str(i+1))).split(',')
            A               = AvsB[0].strip()
            B               = AvsB[1].strip()
            AvsB_FC         = cf.get('2_sample_comparison','2_sample_comparison_'+str(i+1)+'_fc').strip()
            A_vs_B_FCs.append([A, B, AvsB_FC])
    return A_vs_B_FCs


#####################两组非配对
def parse_2group_unpaired_comparison():
    """
    [[[a1,a2],a,[b1,b2],b,2,0.05],...]
    """
    comparisons = cf.options("2_group_unpaired_comparison") # sextuplet
    groupA_vs_groupB_FCs_Ps = []
    if len(comparisons)%6 == 0:
        for i in range(len(comparisons)/6):
            groupAs = remove_space(cf.get('2_group_unpaired_comparison', '2_group_unpaired_comparison_'+str(i+1)+'_groupA_samples')).split(',')      # [A1, A2, A3]
            groupAname = cf.get('2_group_unpaired_comparison', '2_group_unpaired_comparison_'+str(i+1)+'_groupA_name')                          # A
            groupBs = remove_space(cf.get('2_group_unpaired_comparison', '2_group_unpaired_comparison_'+str(i+1)+'_groupB_samples')).split(',')      # [B1, B2, B3]
            groupBname =  cf.get('2_group_unpaired_comparison', '2_group_unpaired_comparison_'+str(i+1)+'_groupB_name')                          # B
            groupA_vs_groupB_FC = cf.get('2_group_unpaired_comparison', '2_group_unpaired_comparison_'+str(i+1)+'_FC')                          # 2.0
            groupA_vs_groupB_P = cf.get('2_group_unpaired_comparison', '2_group_unpaired_comparison_'+str(i+1)+'_PValue')                       # 0.05
            groupA_vs_groupB_FCs_Ps.append([groupAs, groupAname, groupBs, groupBname, groupA_vs_groupB_FC, groupA_vs_groupB_P])
    return groupA_vs_groupB_FCs_Ps


#################### 两组配对的
def parse_2group_paired_comparison():
    comparisons = cf.options("2_group_paired_comparison") # sextuplet
    groupA_vs_groupB_FCs_Ps = []
    if len(comparisons)%6 == 0:
        for i in range(len(comparisons)/6):
            groupAs = remove_space(cf.get('2_group_paired_comparison', '2_group_paired_comparison_'+str(i+1)+'_groupA_samples')).split(',')      # [A1, A2, A3]
            groupAname = cf.get('2_group_paired_comparison', '2_group_paired_comparison_'+str(i+1)+'_groupA_name')                          # A
            groupBs = remove_space(cf.get('2_group_paired_comparison', '2_group_paired_comparison_'+str(i+1)+'_groupB_samples')).split(',')      # [B1, B2, B3]
            groupBname =  cf.get('2_group_paired_comparison', '2_group_paired_comparison_'+str(i+1)+'_groupB_name')                          # B
            groupA_vs_groupB_FC = cf.get('2_group_paired_comparison', '2_group_paired_comparison_'+str(i+1)+'_FC')                          # 2.0
            groupA_vs_groupB_P = cf.get('2_group_paired_comparison', '2_group_paired_comparison_'+str(i+1)+'_PValue')                       # 0.05
            groupA_vs_groupB_FCs_Ps.append([groupAs, groupAname, groupBs, groupBname, groupA_vs_groupB_FC, groupA_vs_groupB_P])
    return groupA_vs_groupB_FCs_Ps


#######################一个样品vs一组样品比较
def parse_sample_group_comparison():
    comparisons = cf.options("1_sample_vs_1_group_comparison") # sextuplet
    print 'parse_sample_group_comparison: 1_sample_vs_1_group_comparison'
    sampleA_vs_groupB_FCs = []
    if len(comparisons)%5 == 0:
        for i in range(len(comparisons)/5):
            sampleA = remove_space(cf.get('1_sample_vs_1_group_comparison', '1_sample_vs_1_group_comparison_'+str(i+1)+'_sample'))      # A1
            sampleAname = cf.get('1_sample_vs_1_group_comparison', '1_sample_vs_1_group_comparison_'+str(i+1)+'_sample_name')                          # A
            groupBs = remove_space(cf.get('1_sample_vs_1_group_comparison', '1_sample_vs_1_group_comparison_'+str(i+1)+'_group')).split(',')      # [B1, B2, B3]
            groupBname =  cf.get('1_sample_vs_1_group_comparison', '1_sample_vs_1_group_comparison_'+str(i+1)+'_group_name')                          # B
            sampleA_vs_groupB_FC = cf.get('1_sample_vs_1_group_comparison', '1_sample_vs_1_group_comparison_'+str(i+1)+'_FC')                          # 2.0
            sampleA_vs_groupB_FCs.append([sampleA, sampleAname, groupBs, groupBname, sampleA_vs_groupB_FC])
            print 'parse_sample_group_comparison: ',sampleA
    return sampleA_vs_groupB_FCs
#######################一组样品vs一个样品比较
def parse_group_sample_comparison():
    comparisons = cf.options("1_group_vs_1_sample_comparison") # sextuplet
    groupA_vs_sampleB_FCs = []
    if len(comparisons)%5 == 0:
        for i in range(len(comparisons)/5):
            groupAs = remove_space(cf.get('1_group_vs_1_sample_comparison', '1_group_vs_1_sample_comparison_'+str(i+1)+'_group')).split(',')      # [A1, A2, A3]
            groupAname = cf.get('1_group_vs_1_sample_comparison', '1_group_vs_1_sample_comparison_'+str(i+1)+'_group_name')                          # A
            sampleB = remove_space(cf.get('1_group_vs_1_sample_comparison', '1_group_vs_1_sample_comparison_'+str(i+1)+'_sample'))      # B1
            sampleBname =  cf.get('1_group_vs_1_sample_comparison', '1_group_vs_1_sample_comparison_'+str(i+1)+'_sample_name')                          # B
            groupAs_vs_sampleB_FC = cf.get('1_group_vs_1_sample_comparison', '1_group_vs_1_sample_comparison_'+str(i+1)+'_FC')                          # 2.0
            groupA_vs_sampleB_FCs.append([groupAs, groupAname, sampleB, sampleBname, groupAs_vs_sampleB_FC])
    return groupA_vs_sampleB_FCs


###################
def index_location(headlist, samplesname):
    """
    返回sublist中样品在list中的位置（无重复）
    input:
        headlist        原始list
        samplesname     子list
    output:
        locations       返回位置list
    """
    locations = []
    for i in samplesname:
        locations.append(headlist.index(i))
    return locations      


def add_fdr_column(file):
    """
    对genespring的x vs y结果进行fdr修正
    input:
        file    待修正文件（第一行为head，其余行为数据，第二列是p-value（1-based）
    output:
        file+'.fdr.txt' 修正后的文件（FDR排在第三列（1-based)
    """
    f = open(file)
    head = f.readline()
    headlist = head.strip('\n').split('\t')
    f2 = open(file+'.fdr.txt','w')
    print >>f2, '\t'.join(headlist[0:2])+'\tFDR\t'+'\t'.join(headlist[2:])
    dicpp={}
    dicppnew={}
    datalines = f.readlines()
    
    for line in datalines:
        linelist = line.strip('\n').split('\t')
        ps_value=float(linelist[1])
        probename=str(linelist[0])
        dicpp[probename]=ps_value
    f.close()
    probes=dicpp.keys()
    ps=dicpp.values()
    corrected_ps = correct_pvalues_for_multiple_testing(ps)
    p_len=len(probes)
    for i in range(p_len):
    	dicppnew[probes[i]]=corrected_ps[i]
    for i in range(len(datalines)):
        line2list = datalines[i].strip('\n').split('\t')
        print >>f2, '\t'.join(line2list[0:2])+'\t%.9f'%dicppnew[line2list[0]]+'\t'+'\t'.join(line2list[2:])
    f2.close()


def make_comparison2dic(comparison):
	"""
	input:
	[['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05']
	or [a,a,[b1,b2],b,2]
	or [[a1,a2],a,b,b,2]
	output:
	put "samples" in the lists into dictionaries which the value is "group",
	and the form is sample_group_dic[sample]=group
	"""
	sample_group_dic={}
	group1_name=comparison[1]
	group2_name=comparison[3]
	samples1=comparison[0]
	samples2=comparison[2]
	if isinstance(samples1,list):
		for i in samples1:
			sample_group_dic[i]=group1_name
	else:sample_group_dic[samples1]=group1_name
	if isinstance(samples2,list):
		for i in samples2:
			sample_group_dic[i]=group2_name
	else:sample_group_dic[samples2]=group2_name
	return sample_group_dic


def name_group(name,dic):
	"""
	get "groupname" from dic, and combine sample and group as 
	[samplename, groupname](raw) or [samplename, groupname](normalized)
	"""
	samplename=name.split("](")[0][1:]
	typename=name.split("](")[1][:-1]
	try:
		groupname=dic[samplename]
		newname="["+samplename+", "+groupname+"]("+typename+")"
	except:pass
	return newname


def get_two_group_unpaired_comparison(inputfile):
    """
    根据配置文件和输入文件，计算差异
    input: 
        all.mRNA.4out8.txt                             flag filtered data
        [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
    
    output:
        all.mRNA.4out8.txt.A_vs_B
        all.mRNA.4out8.txt.A_vs_B.2.0.0.05.txt
        ...
    """
    comparisonlist = parse_2group_unpaired_comparison() #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
    if comparisonlist:
        for comparison in comparisonlist:
            sample_group_dic=make_comparison2dic(comparison)
            f = open(inputfile) # all.mRNA.4out8.txt
            head = f.readline()
            headlist = head.strip('\n').split('\t')
            groupA_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[0]])
            groupB_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[2]])
            groupA_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[0]])
            groupB_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[2]])
            # seqname, P, FDR, abs(FC), regulation, annotation, raw, normalized
            # since there are 3 repeats of sample numbers.
            favsbtmp = open("unpaired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp' ,'w')
            print >>favsbtmp, headlist[0]+'\t'+'P-value\tFold Change\tRegulation\tgroup-%s(raw)\tgroup-%s(raw)\tgroup-%s(normalized)\tgroup-%s(normalized)' % (comparison[1], comparison[3], comparison[1], comparison[3])+'\t'+'\t'.join([headlist[i] for i in groupA_raw_index])+'\t'+'\t'.join([headlist[i] for i in groupB_raw_index])+'\t'+'\t'.join([headlist[i] for i in groupA_norm_index])+'\t'+'\t'.join([headlist[i] for i in groupB_norm_index])+'\t'+'\t'.join(headlist[1+2*int(user_information_sample_number):])
            for line in f:
                linelist = line.strip('\n').split('\t')
                # just select probe, P, mean_groupA(raw), mean_groupB(raw),mean_groupA(normalized), mean_groupB(normalized), groupA(raw), groupB(raw), groupA(normalized), groupB(normalized), annotation, no fold change here, used for scatter plot [just read the file, and select 4, 5 columns], and FDR calculation
                groupA_raw = [float(linelist[i]) for i in groupA_raw_index]
                groupB_raw = [float(linelist[i]) for i in groupB_raw_index]
                groupA_norm = [float(linelist[i]) for i in groupA_norm_index]
                groupB_norm = [float(linelist[i]) for i in groupB_norm_index]
                groupA_raw_mean = means(linelist, groupA_raw_index)
                groupB_raw_mean = means(linelist, groupB_raw_index)
                groupA_norm_mean= means(linelist, groupA_norm_index)
                groupB_norm_mean= means(linelist, groupB_norm_index)
                p = unpaired_ttest(groupA_norm, groupB_norm)
                if groupA_norm_mean >= groupB_norm_mean:
                    fc = (2**groupA_norm_mean)/(2**groupB_norm_mean)
                    regulation = 'up'
                else:
                    fc = (2**groupB_norm_mean)/(2**groupA_norm_mean)
                    regulation = 'down'
                print >>favsbtmp, linelist[0]+'\t'+'%.13f\t%.7f\t%s\t%.6f\t%.6f\t%.6f\t%.6f' % (p, fc, regulation, groupA_raw_mean, groupB_raw_mean, groupA_norm_mean, groupB_norm_mean)+'\t'+'\t'.join([linelist[i] for i in groupA_raw_index])+'\t'+'\t'.join([linelist[i] for i in groupB_raw_index])+'\t'+'\t'.join([linelist[i] for i in groupA_norm_index])+'\t'+'\t'.join([linelist[i] for i in groupB_norm_index])+'\t'+'\t'.join(linelist[1+2*int(user_information_sample_number):])
            f.close()
            favsbtmp.close()
            # add_fdr_column("unpaired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp')
            # 差异表达
            all_matrix,data_len=all_txt("unpaired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp')
            #### sort matrix by foldchange
            head =all_matrix.pop(0)
            for i in range(0,len(all_matrix)):
				for j in range(0,len(all_matrix[i])):
					all_matrix[i][j]=formats(all_matrix[i][j])
            all_matrix_sorted=sorted(all_matrix,key = lambda all_matrix: all_matrix[2],reverse=True)
            all_matrix_sorted.insert(0,head)
            all_matrix=all_matrix_sorted
            #print all_matrix[0]
            print 'the matrix has been sorted'
            #### sort matrix end
            matrix_len=len(all_matrix)
            f3 = open(inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.'+comparison[5]+'.unpaired.fdr.txt','w')
            fallde = open("allde_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.1.0.1.0.unpaired.fdr.txt','w')
            newheadlist=[]
            p=re.compile('^\[[\S\s]+\]\(\w{3,10}\)$')# blank is allowed for the name or group name
            for i in all_matrix[0]:
            	#if "[" in i and "]" in i and "(" in i and ")" in i:# the sentence is not very strict
            	if p.match(i):
           			newheadlist.append(name_group(i,sample_group_dic))
            	else:
            		newheadlist.append(i)
            print >>f3, "\t".join(newheadlist)
            print >>fallde, "\t".join(newheadlist)
            ########header to f4
            f4 = open('pre_heatmap_'+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.'+comparison[5]+'.unpaired.fdr.txt','w')
            head2list=all_matrix[0]
            groupA2_norm_index = index_location(head2list, ['['+i+'](normalized)' for i in comparison[0]])
            groupB2_norm_index = index_location(head2list, ['['+i+'](normalized)' for i in comparison[2]])
            norm_headlist=[head2list[i] for i in (groupA2_norm_index+groupB2_norm_index)]
            h_headlist=[]
            for i in norm_headlist:
                samplename=i.split("](")[0][1:]
                typename=i.split("](")[1]
                newname="["+samplename+", "+sample_group_dic[samplename] +"]("+typename
                h_headlist.append(newname)
            h_headlist=[head2list[0]]+h_headlist
            print >>f4, '\t'.join(h_headlist)
            for i in range(1,matrix_len):
                line2list=all_matrix[i]
                linelist_str = [str(k) for k in line2list]
                p2 = float(line2list[1])
                fc2 = float(line2list[2])
                cline="\t".join(linelist_str)
                fallde.write("%s\n"%cline)
                if p2<=float(comparison[5]) and fc2 >= float(comparison[4]):
                    #cline="\t".join(linelist_str)
                    f3.write("%s\n"%cline)
                    ################### content to f4
                    newlist1=[line2list[0]]+[str(line2list[i]) for i in (groupA2_norm_index+groupB2_norm_index)]
                    h_cline="\t".join(newlist1)
                    f4.write("%s\n"%h_cline)
            f3.close()
            f4.close()
            fallde.close()
            try:shutil.copy("allde_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.1.0.1.0.unpaired.fdr.txt',all_de_loc)
            except:pass

def get_two_group_paired_comparison(inputfile):
    """
    根据配置文件和输入文件，计算差异
    input: 
        all.mRNA.4out8.txt                             flag filtered data
        [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
    
    output:
        all.mRNA.4out8.txt.A_vs_B
        all.mRNA.4out8.txt.A_vs_B.2.0.0.05.txt
        ...
    """
    comparisonlist = parse_2group_paired_comparison() #  [[['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
    if comparisonlist:
        for comparison in comparisonlist:
            sample_group_dic=make_comparison2dic(comparison)
            f = open(inputfile) # all.mRNA.4out8.txt
            head = f.readline()
            headlist = head.strip('\n').split('\t')
            groupA_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[0]])
            groupB_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[2]])
            groupA_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[0]])
            groupB_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[2]])
            
            # seqname, P, FDR, abs(FC), regulation, annotation, raw, normalized
            # since there are 3 repeats of sample numbers.
            favsbtmp = open("paired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp' ,'w')
            print >>favsbtmp, headlist[0]+'\t'+'P-value\tFold Change\tRegulation\tgroup-%s(raw)\tgroup-%s(raw)\tgroup-%s(normalized)\tgroup-%s(normalized)' % (comparison[1], comparison[3], comparison[1], comparison[3])+'\t'+'\t'.join([headlist[i] for i in groupA_raw_index])+'\t'+'\t'.join([headlist[i] for i in groupB_raw_index])+'\t'+'\t'.join([headlist[i] for i in groupA_norm_index])+'\t'+'\t'.join([headlist[i] for i in groupB_norm_index])+'\t'+'\t'.join(headlist[1+2*int(user_information_sample_number):])
            for line in f:
                linelist = line.strip('\n').split('\t')
                # just select probe, P, mean_groupA(raw), mean_groupB(raw),mean_groupA(normalized), mean_groupB(normalized), groupA(raw), groupB(raw), groupA(normalized), groupB(normalized), annotation, no fold change here, used for scatter plot [just read the file, and select 4, 5 columns], and FDR calculation
                groupA_raw = [float(linelist[i]) for i in groupA_raw_index]
                groupB_raw = [float(linelist[i]) for i in groupB_raw_index]
                groupA_norm = [float(linelist[i]) for i in groupA_norm_index]
                groupB_norm = [float(linelist[i]) for i in groupB_norm_index]
                groupA_raw_mean = means(linelist, groupA_raw_index)
                groupB_raw_mean = means(linelist, groupB_raw_index)
                groupA_norm_mean= means(linelist, groupA_norm_index)
                groupB_norm_mean= means(linelist, groupB_norm_index)
                p = paired_ttest(groupA_norm, groupB_norm)
                if groupA_norm_mean >= groupB_norm_mean:
                    fc = (2**groupA_norm_mean)/(2**groupB_norm_mean)
                    regulation = 'up'
                else:
                    
                    fc = (2**groupB_norm_mean)/(2**groupA_norm_mean)
                    regulation = 'down'
                print >>favsbtmp, linelist[0]+'\t'+'%.13f\t%.7f\t%s\t%.6f\t%.6f\t%.6f\t%.6f' % (p, fc, regulation, groupA_raw_mean, groupB_raw_mean, groupA_norm_mean, groupB_norm_mean)+'\t'+'\t'.join([linelist[i] for i in groupA_raw_index])+'\t'+'\t'.join([linelist[i] for i in groupB_raw_index])+'\t'+'\t'.join([linelist[i] for i in groupA_norm_index])+'\t'+'\t'.join([linelist[i] for i in groupB_norm_index])+'\t'+'\t'.join(linelist[1+2*int(user_information_sample_number):])
            f.close()
            favsbtmp.close()
            #print "add_fdr_column"
            #add_fdr_column("paired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp')
            # 差异表达
            all_matrix,datal_len=all_txt("paired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp')
            #### sort matrix by foldchange
            head =all_matrix.pop(0)
            for i in range(0,len(all_matrix)):
				for j in range(0,len(all_matrix[i])):
					all_matrix[i][j]=formats(all_matrix[i][j])
            all_matrix_sorted=sorted(all_matrix,key = lambda all_matrix: all_matrix[2],reverse=True)
            all_matrix_sorted.insert(0,head)
            all_matrix=all_matrix_sorted
            #print all_matrix[0]
            print 'the matrix has been sorted'
            #### sort matrix end
            matrix_len=len(all_matrix)
            headlist=all_matrix[0]
            f3 = open(inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.'+comparison[5]+'.paired.fdr.txt','w')
            fallde = open("allde_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.1.0.1.0.paired.fdr.txt','w')
            newheadlist=[]
            p=re.compile('^\[[\S\s]+\]\(\w{3,10}\)$')# blank is allowed for the name or group name
            for i in all_matrix[0]:
            	#if "[" in i and "]" in i and "(" in i and ")" in i: #it is not strict
            	if p.match(i):
           			newheadlist.append(name_group(i,sample_group_dic))
            	else:
            		newheadlist.append(i)
            print >>f3, "\t".join(newheadlist)
            print >>fallde, "\t".join(newheadlist)
            ########header to f4
            f4 = open('pre_heatmap_'+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.'+comparison[5]+'.paired.fdr.txt','w')
            head2list=all_matrix[0]
            groupA2_norm_index = index_location(head2list, ['['+i+'](normalized)' for i in comparison[0]])
            groupB2_norm_index = index_location(head2list, ['['+i+'](normalized)' for i in comparison[2]])
            norm_headlist=[head2list[i] for i in (groupA2_norm_index+groupB2_norm_index)]
            h_headlist=[]
            for i in norm_headlist:
                samplename=i.split("](")[0][1:]
                typename=i.split("](")[1]
                newname="["+samplename+", "+sample_group_dic[samplename] +"]("+typename
                h_headlist.append(newname)
            h_headlist=[head2list[0]]+h_headlist
            print >>f4, '\t'.join(h_headlist)
            #print "add_fdr_column"
            all_matrix,datal_len=all_txt("paired_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp')
            matrix_len=len(all_matrix)
            #print matrix_len
            for i in range(1,matrix_len):
                line2list=all_matrix[i]
                linelist_str = [str(k) for k in line2list]
                p2 = float(line2list[1])
                fc2 = float(line2list[2])
                cline="\t".join(linelist_str)
                fallde.write("%s\n"%cline)
                if p2<=float(comparison[5]) and fc2 >= float(comparison[4]):
                    #cline="\t".join(line2list)
                    f3.write("%s\n"%cline)
                    ################### content to f4
                    newlist1=[line2list[0]]+[str(line2list[i]) for i in (groupA2_norm_index+groupB2_norm_index)]
                    h_cline="\t".join(newlist1)
                    f4.write("%s\n"%h_cline)
            f3.close()
            f4.close()
            fallde.close()
            try:shutil.copy("allde_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.1.0.1.0.paired.fdr.txt',all_de_loc)
            except:pass
def two_sample_comparison(inputfile):
    """
    两个样品比较
    input: 
        all.mRNA.2out3.txt                             flag filtered data
        [['A', 'B', '2.0'], ...]
    
    output:
        all.mRNA.2out3.txt.A_vs_B
        all.mRNA.2out3.txt.A_vs_B.2.0.txt
        
        ...
    """    
    comparisonlist = parse_2sample_comparison()
    if comparisonlist:
        for comparison in comparisonlist:
            print "two_sample_comparison"
            f = open(inputfile) # all.mRNA.4out8.txt            
            head = f.readline()
            headlist = head.strip('\n').split('\t')
            A_raw_index = headlist.index('['+comparison[0]+'](raw)')
            B_raw_index = headlist.index('['+comparison[1]+'](raw)')
            A_norm_index = headlist.index('['+comparison[0]+'](normalized)')
            B_norm_index = headlist.index('['+comparison[1]+'](normalized)')
            favsbtmp = open(inputfile+'.'+comparison[0]+'_vs_'+comparison[1]+'.tmp' ,'w')
            print >>favsbtmp, headlist[0]+'\t'+ 'Fold Change\tRegulation\t%s\t%s\t%s\t%s\t%s' %  (headlist[A_raw_index], headlist[B_raw_index], headlist[A_norm_index], headlist[B_norm_index],'\t'.join(headlist[1+2*int(user_information_sample_number):]))
            for line in f:
                linelist = line.strip('\n').split('\t')
                # just select probe, P, mean_groupA(raw), mean_groupB(raw),mean_groupA(normalized), mean_groupB(normalized), groupA(raw), groupB(raw), groupA(normalized), groupB(normalized), annotation, no fold change here, used for scatter plot [just read the file, and select 4, 5 columns], and FDR calculation
                A_norm_value = float(linelist[A_norm_index])
                B_norm_value = float(linelist[B_norm_index])
                if A_norm_value >= B_norm_value:
                    fc = (2**A_norm_value)/(2**B_norm_value)
                    abs_fc = fc
                    regulation = 'up'
                else:
                    fc = - (2**B_norm_value)/(2**A_norm_value)
                    abs_fc = - fc
                    regulation = 'down'
                print >>favsbtmp, linelist[0]+'\t%.7f\t%s\t%s\t%s\t%s\t%s\t%s' % (abs_fc,regulation, linelist[A_raw_index], linelist[B_raw_index], linelist[A_norm_index], linelist[B_norm_index],'\t'.join(map(str, linelist[1+2*int(user_information_sample_number):])))
            f.close()
            favsbtmp.close()
            all_matrix,data_len=all_txt(inputfile+'.'+comparison[0]+'_vs_'+comparison[1]+'.tmp')
            #### sort matrix by foldchange
            head =all_matrix.pop(0)
            for i in range(0,len(all_matrix)):
				for j in range(0,len(all_matrix[i])):
					all_matrix[i][j]=formats(all_matrix[i][j])
            all_matrix_sorted=sorted(all_matrix,key = lambda all_matrix: all_matrix[1],reverse=True)
            all_matrix_sorted.insert(0,head)
            all_matrix=all_matrix_sorted
            #print all_matrix[0]
            print 'the matrix has been sorted'
            #### sort matrix end
            matrix_len=len(all_matrix)
            headlist=all_matrix[0]
            A_raw_index = headlist.index('['+comparison[0]+'](raw)')
            B_raw_index = headlist.index('['+comparison[1]+'](raw)')
            A_norm_index = headlist.index('['+comparison[0]+'](normalized)')
            B_norm_index = headlist.index('['+comparison[1]+'](normalized)')
            # seqname, fc, abs(fc), log(fc), regulation, annotation, raw, normalized
            # since there are 3 repeats of sample numbers.
            favsb = open(inputfile+'.'+comparison[0]+'_vs_'+comparison[1] ,'w')
            print >>favsb, '\t'.join(headlist)
            fupdown  = open(inputfile+'.'+comparison[0]+'_vs_'+comparison[1]+'.'+comparison[2]+'.txt','w')
            fallde	=	open("allde_"+inputfile+'.'+comparison[0]+'_vs_'+comparison[1]+'.1.0.txt','w')
            print >>fupdown, '\t'.join(headlist)
            print >>fallde, '\t'.join(headlist)
            f4 = open('pre_heatmap_'+inputfile+'.'+comparison[0]+'_vs_'+comparison[1]+'.'+comparison[2]+'.txt','w')
            print >>f4, headlist[0]+'\t'+headlist[A_norm_index]+'\t'+headlist[B_norm_index]
            for i in range(1,matrix_len):
                linelist=all_matrix[i]
                # just select probe, a, b, annotation, no fold change here, used for scatter plot [just read the file, and select 3, 4 columns]
                #print linelist
                linelist_str = [str(k) for k in linelist]
                print >>favsb, '\t'.join(linelist_str)
                A_norm_value = float(linelist[A_norm_index])
                B_norm_value = float(linelist[B_norm_index])
                # up
                if A_norm_value >= B_norm_value:
                    fc = (2**A_norm_value)/(2**B_norm_value)
                    abs_fc = fc
                    #logfc = log2(fc)
                    if abs_fc >= float(comparison[2]):
                        print >>fupdown, '\t'.join(linelist_str)
                        print >>f4, linelist[0]+"\t"+linelist_str[A_norm_index]+"\t"+linelist_str[B_norm_index]
                    if abs_fc >= 1.0:
                        print >>fallde, '\t'.join(linelist_str)

                # down
                else:
                    fc = - (2**B_norm_value)/(2**A_norm_value)
                    abs_fc = - fc
                    #logfc = -log2(-fc)
                    regulation = 'down'
                    if fc < -float(comparison[2]):
                        print >>fupdown, '\t'.join(linelist_str)
                        print >>f4, linelist[0]+"\t"+linelist_str[A_norm_index]+"\t"+linelist_str[B_norm_index]
                    if fc < -1.0:
                        print >>fallde, '\t'.join(linelist_str)
            favsb.close()
            fupdown.close()
            f4.close()
            fallde.close()
            print 'two_sample_comparison end'
            try:
				shutil.copy("allde_"+inputfile+'.'+comparison[0]+'_vs_'+comparison[1]+'.1.0.txt',all_de_loc)
				print 'copy allde'
            except:
				print "error copy allde";
def one_sg_one_sg_comparison(inputfile,comparisonlist):
    """
    两个样品-组/组样品比较#one_sg_one_sg_comparison
    input: 
        all.mRNA.2out3.txt                             flag filtered data
        parse_sample_group_comparison():[[a,a,[b1,b2],b,2],...]; parse_group_sample_comparison():[[[a1,a2],a,b,b,2],...]
    
    output:
        all.mRNA.2out3.txt.A_vs_B
        all.mRNA.2out3.txt.A_vs_B.2.0.txt
        
        ...
    """
    print 'one_sg_one_sg_comparison:',' ',comparisonlist
    if comparisonlist:
        for comparison in comparisonlist:
            print "one_sg_one_sg_comparison"
            sample_group_dic=make_comparison2dic(comparison)
            f = open(inputfile) # all.mRNA.4out8.txt
            print "inputfile",inputfile;
            print "comparison",comparison
            print "user_information_sample_number",user_information_sample_number
            head = f.readline()
            headlist = head.strip('\n').split('\t')
            if isinstance(comparison[0], list):
            	A_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[0]])
            	A_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[0]])
            else:
                A_raw_index = index_location(headlist, ['['+comparison[0]+'](raw)'])
                A_norm_index = index_location(headlist, ['['+comparison[0]+'](normalized)'])
            if isinstance(comparison[2], list):
                B_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[2]])
                B_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[2]])
            else:            
                B_raw_index = index_location(headlist, ['['+comparison[2]+'](raw)'])
                B_norm_index = index_location(headlist, ['['+comparison[2]+'](normalized)'])
            # seqname, P, FDR, abs(FC), regulation, annotation, raw, normalized
            # since there are 3 repeats of sample numbers.
            print "one_sg_one_sg_comparison lable001"
            favsbtmp = open(inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp' ,'w')
            print >>favsbtmp, headlist[0]+'\t'+'Fold Change\tRegulation'+'\t'+'\t'.join([headlist[i] for i in A_raw_index])+'\t'+'\t'.join([headlist[i] for i in B_raw_index])+'\t'+'\t'.join([headlist[i] for i in A_norm_index])+'\t'+'\t'.join([headlist[i] for i in B_norm_index])+'\t'+'\t'.join(headlist[1+2*int(user_information_sample_number):])
            print "one_sg_one_sg_comparison lable002"
            for line in f:
                linelist = line.strip('\n').split('\t')
                # just select probe, P, mean_groupA(raw), mean_groupB(raw),mean_groupA(normalized), mean_groupB(normalized), groupA(raw), groupB(raw), groupA(normalized), groupB(normalized), annotation, no fold change here, used for scatter plot [just read the file, and select 4, 5 columns], and FDR calculation
                groupA_raw = [float(linelist[i]) for i in A_raw_index]
                groupB_raw = [float(linelist[i]) for i in B_raw_index]
                groupA_norm = [float(linelist[i]) for i in A_norm_index]
                groupB_norm = [float(linelist[i]) for i in B_norm_index]
                groupA_raw_mean = means(linelist, A_raw_index)
                groupB_raw_mean = means(linelist, B_raw_index)
                groupA_norm_mean= means(linelist, A_norm_index)
                groupB_norm_mean= means(linelist, B_norm_index)
                #print "one_sg_one_sg_comparison lable004"
                if groupA_norm_mean >= groupB_norm_mean:
                    fc = (2**groupA_norm_mean)/(2**groupB_norm_mean)
                    regulation = 'up'
                else:
                    fc = (2**groupB_norm_mean)/(2**groupA_norm_mean)
                    regulation = 'down'
                print >>favsbtmp, linelist[0]+'\t'+'%.7f\t%s' % (fc, regulation)+'\t'+'\t'.join([linelist[i] for i in A_raw_index])+'\t'+'\t'.join([linelist[i] for i in B_raw_index])+'\t'+'\t'.join([linelist[i] for i in A_norm_index])+'\t'+'\t'.join([linelist[i] for i in B_norm_index])+'\t'+'\t'.join(linelist[1+2*int(user_information_sample_number):])
            f.close()
            favsbtmp.close()
            #print "one_sg_one_sg_comparison lable003"
            all_matrix,data_len=all_txt(inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.tmp')
            #### sort matrix by foldchange
            head =all_matrix.pop(0)
            for i in range(0,len(all_matrix)):
				for j in range(0,len(all_matrix[i])):
					all_matrix[i][j]=formats(all_matrix[i][j])
            all_matrix_sorted=sorted(all_matrix,key = lambda all_matrix: all_matrix[1],reverse=True)
            all_matrix_sorted.insert(0,head)
            all_matrix=all_matrix_sorted
            #print all_matrix[0]
            print 'the matrix has been sorted'
            #### sort matrix end
            matrix_len=len(all_matrix)
            headlist=all_matrix[0]
            if isinstance(comparison[0], list):
            	A_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[0]])
            	A_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[0]])
            else:
                A_raw_index = index_location(headlist, ['['+comparison[0]+'](raw)'])
                A_norm_index = index_location(headlist, ['['+comparison[0]+'](normalized)'])
            if isinstance(comparison[2], list):
                B_raw_index = index_location(headlist, ['['+i+'](raw)' for i in comparison[2]])
                B_norm_index = index_location(headlist, ['['+i+'](normalized)' for i in comparison[2]])
            else:            
                B_raw_index = index_location(headlist, ['['+comparison[2]+'](raw)'])
                B_norm_index = index_location(headlist, ['['+comparison[2]+'](normalized)'])
            ########header to f4
            sample_number=len(A_raw_index)+len(B_raw_index)
            print 'sample_number',sample_number;
            f4 = open('pre_heatmap_'+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.txt','w')
            norm_headlist=[headlist[i] for i in (A_norm_index+B_norm_index)]
            h_headlist=[]
            for i in norm_headlist:
            	samplename=i.split("](")[0][1:]
            	typename=i.split("](")[1]
            	newname="["+samplename+", "+sample_group_dic[samplename] +"]("+typename
            	h_headlist.append(newname)
            h_headlist=[headlist[0]]+h_headlist
            print >>f4, '\t'.join(h_headlist)
            ###########
            favsb = open(inputfile+'.'+comparison[1]+'_vs_'+comparison[3] ,'w')
            print >>favsb, "\t".join(headlist)
            fupdown  = open(inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.txt','w')
            #try:
            print >>fupdown, "\t".join(headlist)
            fallde  = open("allde_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.1.0.txt','w')
            print >>fallde, "\t".join(headlist)
            fscatter=open("scatter_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3] ,'w')
            print >>fscatter, headlist[0]+'\t'+'['+comparison[1]+'](raw)'+'\t'+'['+comparison[3]+'](raw)'+'\t'+'['+comparison[1]+'](normalized)'+'\t'+'['+comparison[3]+'](normalized)'+'\t'+'\t'.join(headlist[3+2*int(sample_number):])
            for i in range(1,matrix_len):
            	linelist=all_matrix[i]
            	linelist_str = [str(k) for k in linelist]
                A_raw = [str(linelist[i]) for i in A_raw_index]
                B_raw = [str(linelist[i]) for i in B_raw_index]
                A_norm = [str(linelist[i]) for i in A_norm_index]
                B_norm = [str(linelist[i]) for i in B_norm_index]
                A_norm_mean= means(linelist, A_norm_index)
                A_raw_mean=means(linelist, A_raw_index)
                B_norm_mean= means(linelist, B_norm_index)
                B_raw_mean= means(linelist, B_raw_index)
                print >>favsb, "\t".join(linelist_str)
                print >>fscatter,linelist[0]+'\t'+str(A_raw_mean)+'\t'+str(B_raw_mean)+'\t'+str(A_norm_mean)+'\t'+str(B_norm_mean)+'\t'+'\t'.join(linelist_str[3+2*int(sample_number):])
                # up
                #print 'one_sg_one_sg_comparison lablel001'
                if A_norm_mean >= B_norm_mean:
                    fc = (2**A_norm_mean)/(2**B_norm_mean)
                    abs_fc = fc
                    logfc = log2(fc)
                    if abs_fc >= float(comparison[4]):
                        #print 'one_sg_one_sg_comparison lablel002'
                        print >>fupdown, "\t".join(linelist_str)
                        ################### content to f4
                        newlist1=[linelist[0]]+[str(linelist[i]) for i in (A_norm_index+B_norm_index)]
                        print >>f4,"\t".join(newlist1)
                    if abs_fc >= 1.0:
                        print >>fallde, "\t".join(linelist_str)

                        #########################
                # down
                else:
                    fc = - (2**B_norm_mean)/(2**A_norm_mean)
                    abs_fc = - fc
                    logfc = -log2(-fc)
                    regulation = 'down'
                    if fc < -float(comparison[4]):
                        print >>fupdown, "\t".join(linelist_str)
                        ################### content to f4
                        newlist1=[linelist[0]]+[str(linelist[i]) for i in (A_norm_index+B_norm_index)]
                        print >>f4,"\t".join(newlist1)
                        #########################
                    if fc < -1.0:
                        print >>fallde, "\t".join(linelist_str)
            print 'one_sg_one_sg_comparison lablel100'
            favsb.close()
            fupdown.close()
            fscatter.close()
            f4.close()
            fallde.close()
            try:shutil.copy("allde_"+inputfile+'.'+comparison[1]+'_vs_'+comparison[3]+'.1.0.txt',all_de_loc)
            except:pass

###### 写入到excel
def all_txt(file):
	"""
	it can make a file into a list matrix, and the form
	[[rowele1,rowele2,...],...]
	all_txt has the remove control type function
	
	"""
	#print "all_txt1"	
	print 'all_txt',file
	ofile=open(file,"rU")
	print "all_txt2"
	head=ofile.readline().rstrip("\n")
	headlist=re.split("\t",head)
	ct_loc=""
	temp_rawhead_list=[]
	all_matrix=[]
	row_len=len(headlist)
	for i in headlist:
		if "(raw)" in i:
			temp_rawhead_list.append(i)
		j=i.replace(flagsuffix,"")
		if ct_loc!="" and j!=headlist[ct_loc]:
			try:all_matrix[0].append(j)
			except:all_matrix.append([j])
		if ct_loc=="":
			try:all_matrix[0].append(j)
			except:all_matrix.append([j])
	datalen=len(temp_rawhead_list)
	temp_data_list=[]
	line="ini"
	print line
	while line:
		line=ofile.readline().rstrip("\n")
		if line=="":
			break
		linelist=line.split("\t")
		linelist=[t.strip(" ") for t in linelist]
		if ct_loc!="":
			if (linelist[ct_loc]).lower()=="false":
				if ct_loc+1<row_len:
					linelist=linelist[:ct_loc]+linelist[ct_loc+1:]
					temp_data_list.append(linelist)
				elif ct_loc+1==row_len:
					linelist=linelist[:ct_loc]
					temp_data_list.append(linelist)
		else:
			temp_data_list.append(linelist)
			#print '\n\nline',linelist,'\n\n'
	all_matrix=[all_matrix[0]]+temp_data_list
	ofile.close()
	print 'all_txt end'
	return all_matrix,datalen
	
def all_mRNA_excel(all_mRNA_matrix,datalen):
	print "haha"	
	rownum=len(all_mRNA_matrix)
	colnum=len(all_mRNA_matrix[0])
	probename=all_mRNA_matrix[0][0]
	annotation=all_mRNA_matrix[0][2*datalen+1:]
	print annotation
	workbook=xlsxwriter.Workbook('CircRNA Expression Profiling.xlsx')
	workbook.use_zip64()
	worksheet=workbook.add_worksheet('circRNA Expression Profiling')	
	if "P" in all_mRNA_matrix[1] or "M" in all_mRNA_matrix[1] or "A" in all_mRNA_matrix[1]:
		a=str("All Targets Value (Entities where at least ")+outof.replace('outof', ' out of ')+" samples have flags in Present or Marginal)"
	else:a="circRNAs identified by DCC"
	b="Column A: CircRNAID, the ID of the identified circRNA by DCC."
	c="Column B ~ "+ column_index[datalen+1]+": junction reads, the junction read number of each sample."
	d="Column "+column_index[datalen+2]+" ~ "+column_index[2*datalen+1]+": Normalized Intensity, Normalized Intensity of each sample (log2 transformed)."
	e="Column "+str(column_index[datalen*2+2])+" ~ "+str(column_index[datalen*2+5])+": the coordinates of circRNA."
	f="Column "+ str(column_index[datalen*2+6])+": circBaseID, the identifier of circBase (http://www.circbase.org)."
	g="Column "+ str(column_index[datalen*2+7])+": source, the source of the circRNA, including circBase, Guojunjie2014,..."
	h="Column "+ str(column_index[datalen*2+8])+": best_transcript, the best transcript of the circRNA."
	i="Column "+ str(column_index[datalen*2+9])+": GeneName, the name of the circRNA-associated gene."
	j="Column "+ str(column_index[datalen*2+10])+": Catalog, the catalog of the circRNA, including exonic, intronic, ..."
	k="Column "+ str(column_index[datalen*2+11])+": predicted_sequence_length, the length of predicted circRNA sequence."
	if user_information_species.lower() == 'human':
				l="Column "+ str(column_index[datalen*2+12])+": circRNA-associated diseases indicated by (http://gyanxet-beta.com/circdb/)."
	else:
				l=""
	annotationformat =  workbook.add_format({'align':'left', 'fg_color':'#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_align('bottom')
	annotationformat.set_text_wrap()  # auto wrapping
	font0 = workbook.add_format({'font_name':'Times New Roman'})
	font0.set_font_size(10)
	font1 = workbook.add_format({'font_name':'Times New Roman'})
	font1.set_bold('true')
	font1.set_font_size(10)
	annotationformat1 =  workbook.add_format({'align':'center', 'fg_color':'#3366ff','font_name':'Times New Roman','bold':'ture'})
	annotationformat2 =  workbook.add_format({'align':'center', 'fg_color':'#ff9900','font_name':'Times New Roman','bold':'ture'})
	annotationformat3 =  workbook.add_format({'align':'center', 'fg_color':'#008080','font_name':'Times New Roman','bold':'ture'})
	print "merge1..............."
	worksheet.merge_range(0,0,13+mRNA_info_addline,colnum-1,'',font0)
	print "merge2..............."
	worksheet.write_rich_string(0,0,font1,a+'\n'+'\n',font0,b+'\n',font0,c+'\n',font0,d+'\n',font0,e+'\n',font0,f+'\n',font0,g+'\n',font0,h+'\n',font0,i+'\n',font0,j+'\n',font0,k+'\n',font0,l,annotationformat)	
	worksheet.merge_range(15+mRNA_info_addline,1,15+mRNA_info_addline,datalen,'junction reads',annotationformat1)
	worksheet.merge_range(15+mRNA_info_addline, datalen+1,15+mRNA_info_addline,   2*datalen,'Normalized Intensity',annotationformat2)
	worksheet.merge_range(15+mRNA_info_addline, 2*datalen+1,15+mRNA_info_addline,   colnum-1,'Annotations',annotationformat3)
	print "merge3..............."
	worksheet.set_column(0, 0, 30)
	for rows in range(rownum):
		for cols in range(colnum):
			if rows==0:
				worksheet.write(rows+10+6+mRNA_info_addline, cols, all_mRNA_matrix[rows][cols], font1)
			else:
				worksheet.write(rows+10+6+mRNA_info_addline, cols, formats(all_mRNA_matrix[rows][cols]), font0)
	print "merge4..............."
	workbook.close()

def rearrange_list(your_list,newindex):
	newlist=[]
	for i in newindex:
		newlist.append(your_list[i])
	return newlist
def rearrange_de_matrix(de_matrix):
	datalocationlist=[]
	textlocationlist=[]
	new_de_martrix=[]
	for i in de_matrix[0]:
		if "(raw)" in i or "(normalized)" in i:
			datalocationlist.append(de_matrix[0].index(i))
		else:textlocationlist.append(de_matrix[0].index(i))
	locationlist=textlocationlist+datalocationlist
#	print locationlist,textlocationlist,datalocationlist
	for i in de_matrix:
		j=rearrange_list(i,locationlist)
		new_de_martrix.append(j)
	return new_de_martrix
def sep_regulation_matrix(your_arranged_matrix):
	up=[]
	down=[]
	up.append(your_arranged_matrix[0])
	down.append(your_arranged_matrix[0])
	for i in your_arranged_matrix[0]:
		regulation_index=your_arranged_matrix[0].index("Regulation")
	for i in your_arranged_matrix:
		if i[regulation_index]=="up":
			up.append(i)
		if i[regulation_index]=="down":
			down.append(i)
	return up,down
	
def conclude_comparison_info():#conclude each comparison like[A_vs_B,fc,p,sta,repeatorder]
	ss01=parse_2sample_comparison() #[A, B, AvsB_FC]
	ss11=[]
	ss11_dic={}
	temp={}
	if ss01:
		order=1
		for i in ss01:
			try:temp[i[0]+"_vs_"+i[1]].append([order,ss01.index(i)])
			except:temp[i[0]+"_vs_"+i[1]]=[[order,ss01.index(i)]]
			order=order+1
		for k,v in temp.items():
			if len(v)>=2:
				for j in v:
					ss11_dic[j[0]]=[k,ss01[j[1]][2],None,None,ss01[j[1]][2]+"_"]
			else:ss11_dic[v[0][0]]=[k,ss01[v[0][1]][2],None,None,None]  #{sampleorder:[A_vs_B,fc,p,sta,repeatorder],...}, in the dict,sta :paired(pd), unpaired(upd),None; repeatorder:fc_,fc_p_sta_,None
		ss11=[value for (key,value) in sorted(ss11_dic.items())]
	ggu01=parse_2group_unpaired_comparison() #groupA_vs_groupB_FCs_Ps.append([groupAs, groupAname, groupBs, groupBname, groupA_vs_groupB_FC, groupA_vs_groupB_P])
	ggu11=[]
	ggu11_dic={}
	temp={}
	if ggu01:
		order=1
		for i in ggu01:
			try:temp[i[1]+"_vs_"+i[3]].append([order,ggu01.index(i)])
			except:temp[i[1]+"_vs_"+i[3]]=[[order,ggu01.index(i)]]
			order=order+1
		for k,v in temp.items():
			if len(v)>=2:
				for j in v:
					ggu11_dic[j[0]]=[k,ggu01[j[1]][4],ggu01[j[1]][5],"unpaired",ggu01[j[1]][4]+"_"+ggu01[j[1]][5]+"_"+"upd"]
			else:ggu11_dic[v[0][0]]=[k,ggu01[v[0][1]][4],ggu01[v[0][1]][5],"unpaired",None]  #{sampleorder:[A_vs_B,fc,p,sta,repeatorder],...}, in the dict,sta :paired, unpaired,None; repeatorder:fc_,fc_p_sta_,None
		ggu11=[value for (key,value) in sorted(ggu11_dic.items())]
	ggp01=parse_2group_paired_comparison() #groupA_vs_groupB_FCs_Ps.append([groupAs, groupAname, groupBs, groupBname, groupA_vs_groupB_FC, groupA_vs_groupB_P])
	ggp11=[]
	ggp11_dic={}
	temp={}
	if ggp01:
		order=1
		for i in ggp01:
			try:temp[i[1]+"_vs_"+i[3]].append([order,ggp01.index(i)])
			except:temp[i[1]+"_vs_"+i[3]]=[[order,ggp01.index(i)]]
			order=order+1
		for k,v in temp.items():
			if len(v)>=2:
				for j in v:
					ggp11_dic[j[0]]=[k,ggp01[j[1]][4],ggp01[j[1]][5],"paired",ggp01[j[1]][4]+"_"+ggp01[j[1]][5]+"_"+"pd"]
			else:ggp11_dic[v[0][0]]=[k,ggp01[v[0][1]][4],ggp01[v[0][1]][5],"paired",None]  #{sampleorder:[A_vs_B,fc,p,sta,repeatorder],...}, in the dict,sta :paired, unpaired,None; repeatorder:fc_,fc_p_sta_,None
		ggp11=[value for (key,value) in sorted(ggp11_dic.items())]
	sg01=parse_sample_group_comparison() #sampleA_vs_groupB_FCs.append([sampleA, sampleAname, groupBs, groupBname, sampleA_vs_groupB_FC])
	sg11=[]
	sg11_dic={}
	temp={}
	if sg01:
		order=1
		for i in sg01:
			try:temp[i[1]+"_vs_"+i[3]].append([order,sg01.index(i)])
			except:temp[i[1]+"_vs_"+i[3]]=[[order,sg01.index(i)]]
			order=order+1
		for k,v in temp.items():
			if len(v)>=2:
				for j in v:
					sg11_dic[j[0]]=[k,sg01[j[1]][4],None,None,sg01[j[1]][4]+"_"]
			else:sg11_dic[v[0][0]]=[k,sg01[v[0][1]][4],None,None,None]  #{sampleorder:[A_vs_B,fc,p,sta,repeatorder],...}, in the dict,sta :paired, unpaired,None; repeatorder:fc_,fc_p_sta_,None
		sg11=[value for (key,value) in sorted(sg11_dic.items())]
	gs01=parse_group_sample_comparison() #groupA_vs_sampleB_FCs.append([groupAs, groupAname, sampleB, sampleBname, groupAs_vs_sampleB_FC])
	gs11=[]
	gs11_dic={}
	temp={}
	if gs01:
		order=1
		for i in gs01:
			try:temp[i[1]+"_vs_"+i[3]].append([order,gs01.index(i)])
			except:temp[i[1]+"_vs_"+i[3]]=[[order,gs01.index(i)]]
			order=order+1
		for k,v in temp.items():
			if len(v)>=2:
				for j in v:
					gs11_dic[j[0]]=[k,gs01[j[1]][4],None,None,gs01[j[1]][4]+"_"]
			else:gs11_dic[v[0][0]]=[k,gs01[v[0][1]][4],None,None,None]  #{sampleorder:[A_vs_B,fc,p,sta,repeatorder],...}, in the dict,sta :paired, unpaired,None; repeatorder:fc_,fc_p_sta_,None
		gs11=[value for (key,value) in sorted(gs11_dic.items())]
	#####GET SORT INFO
	possible_dic={"[2_sample_comparison]":"ss","[2_group_unpaired_comparison]":"unpaired","[2_group_paired_comparison]":"paired","[1_sample_vs_1_group_comparison]":"sg","[1_group_vs_1_sample_comparison]":"gs"}
	file_config=open(config_file,"rU")
	line="ini"
	temp_comparison_list=[]
	all_lines=file_config.readlines()
	for line in all_lines:
		line=line.rstrip("\n")
		if line.startswith("["):
			if line in possible_dic.keys():
				temp_comparison_list.append(possible_dic[line])
	file_config.close()
	u_p_sta=[]
	info_list=[]
	for i in temp_comparison_list:
		dic={}
		if i =="paired" or i=="unpaired":
			u_p_sta.append(i)
		if i=="ss":
			dic[i]=ss11
		elif i=="sg":
			dic[i]=sg11
		elif i=="gs":
			dic[i]=gs11
		elif i=="paired":
			dic[i]=ggp11
		elif i=="unpaired":
			dic[i]=ggu11
		info_list.append(dic)
	for i in info_list:
		if i.keys()[0]==u_p_sta[0] or i.keys()[0]==u_p_sta[1]:
			try:gg01=gg01+i.values()[0]
			except:gg01=i.values()[0]
	gg11=[]
	gg11_dic={}
	temp={}
	if ggu11 and ggp11:
		order=1
		for i in gg01:
			try:temp[i[0]].append([order,gg01.index(i)])
			except:temp[i[0]]=[[order,gg01.index(i)]]
			order=order+1
		for k,v in temp.items():
			if len(v)>=2:
				for j in v:
					if gg01[j[1]][3]=="paired":
						gg11_dic[j[0]]=[k,gg01[j[1]][1],gg01[j[1]][2],gg01[j[1]][3],gg01[j[1]][1]+"_"+gg01[j[1]][2]+"_"+"pd"]
					if gg01[j[1]][3]=="unpaired":
						gg11_dic[j[0]]=[k,gg01[j[1]][1],gg01[j[1]][2],gg01[j[1]][3],gg01[j[1]][1]+"_"+gg01[j[1]][2]+"_"+"upd"]
			else:gg11_dic[v[0][0]]=[k,gg01[v[0][1]][1],gg01[v[0][1]][2],gg01[v[0][1]][3],gg01[v[0][1]][4]]  #[A_vs_B,fc,p,sta,repeatorder], in the list,sta :paired, unpaired,None; repeatorder:fc_,fc_p_sta_,None
		gg11=[value for (key,value) in sorted(gg11_dic.items())]
	if gg11:
		ggu11=[];ggp11=[];u_p_dic={}
		for j in gg11:
			if j[3]=="paired":
				try:u_p_dic["paired"].append(j)
				except:u_p_dic["paired"]=[j]
			elif j[3]=="unpaired":
				try:u_p_dic["unpaired"].append(j)
				except:u_p_dic["unpaired"]=[j]
		for i in u_p_sta:
			if i=="paired":
				ggp11=u_p_dic[i]
			elif i=="unpaired":
				ggu11=u_p_dic[i]
		for i in info_list:
			if i.keys()[0]=="paired":
				i["paired"]=ggp11
			if i.keys()[0]=="unpaired":
				i["unpaired"]=ggu11
	for i in info_list:
		try:conclusion=conclusion+i.values()[0]
		except:conclusion=i.values()[0]
	return conclusion

def judge_cc(cc):
	cpdic={}
	for i in cc:
		cpdic[i[0]]=i[0][:26]
	a=cpdic.values()
	t=[x for x, y in collections.Counter(a).items() if y > 1]
	b=[]
	for k,v in cpdic.items():
		if len(t)>0:
			for ti in t:
				if ti==v:
					b.append(k)
	b=sorted(set(b))
	if len(b)>0:
		bstring="\n".join(b)
		#printError(u'发现分组\n'+bstring+u'\n中前26个字母有重复，这会因DE excel sheet的名字重复导致生成的DE excel失败！程序跳出，请修正！')
		sys.exit()

def seperate_cc():
	cc= conclude_comparison_info()
	ssn=[]
	sgn=[]
	gsn=[]
	ggpn=[]
	ggupn=[]
	ss=parse_2sample_comparison()
	sg=parse_sample_group_comparison()
	gs=parse_group_sample_comparison()
	ggp = parse_2group_paired_comparison() #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
	ggup = parse_2group_unpaired_comparison() #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
	for i in ss:
		if len(i)==0:
			break
		for j in cc:
			if (i[0]+"_vs_"+i[1]) == j[0] and i[2]==j[1]:
				ssn.append(i+[j[-1]])
	for i in sg:
		if len(i)==0:
			break
		for j in cc:
			if (i[1]+"_vs_"+i[3])==j[0] and i[4]==j[1]:
				sgn.append(i+[j[-1]])
	for i in gs:
		if len(i)==0:
			break
		for j in cc:
			if (i[1]+"_vs_"+i[3])==j[0] and i[4]==j[1]:
				gsn.append(i+[j[-1]])
	for i in ggp:
		if len(i)==0:
			break
		for j in cc:
			if (i[1]+"_vs_"+i[3])==j[0] and i[4]==j[1] and i[5] == j[2] and j[3]=="paired":
				ggpn.append(i+[j[-1]])
	for i in ggup:
		if len(i)==0:
			break
		for j in cc:
			if (i[1]+"_vs_"+i[3])==j[0] and i[4]==j[1] and i[5] == j[2] and j[3]=="unpaired":
				ggupn.append(i+[j[-1]])
	return ssn,sgn,gsn,ggpn,ggupn

def demRNA_excel(tag,cc,location):
	os.chdir(location)
	filenote=open("note.txt","w")
	filenote.write("%s\t%s\t%s\t%s\t%s\n"%("comparison_name","DE_probe_number","fc","p_value","statistic method"))
	#cc=conclude_comparison_info()
	print cc
	judge_cc(cc)
	plist=[]
	for i in cc:
		plist.append(i[2])	
	if None in plist:
		if tag=="allde_":
			wb=xlsxwriter.Workbook('All Comparisons.xlsx')
		else:wb=xlsxwriter.Workbook('Differentially Expressed circRNAs.xlsx')
	else:
		if tag=="allde_":
			wb=xlsxwriter.Workbook('All Comparisons.xlsx')
		else:
			wb=xlsxwriter.Workbook('Differentially Expressed circRNAs.xlsx')
	print 'de1...........................'
	wb.use_zip64()
	for i in cc:
		comparison_name=i[0];fc=i[1];p=i[2];sta=i[3];repeatorder=i[4]
		if p:
			filename=tag+'all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison_name+"."+str(fc)+"."+str(p)+"."+sta+".fdr.txt"
			print filename,'haha'		
		else:
			filename=tag+'all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison_name+"."+str(fc)+".txt"
		if repeatorder:
			if p:
				go_upfile="go_up_"+comparison_name+"."+str(fc)+"."+str(p)+"."+sta+".txt"
				go_dnfile="go_down_"+comparison_name+"."+str(fc)+"."+str(p)+"."+sta+".txt"
				pw_upfile="pathway_up_"+comparison_name+"."+str(fc)+"."+str(p)+"."+sta+".txt"
				pw_dnfile="pathway_down_"+comparison_name+"."+str(fc)+"."+str(p)+"."+sta+".txt"
			else:
				go_upfile="go_up_"+comparison_name+"."+str(fc)+".txt"
				go_dnfile="go_down_"+comparison_name+"."+str(fc)+".txt"
				pw_upfile="pathway_up_"+comparison_name+"."+str(fc)+".txt"
				pw_dnfile="pathway_down_"+comparison_name+"."+str(fc)+".txt"
		else:
			go_upfile="go_up_"+comparison_name+".txt"
			go_dnfile="go_down_"+comparison_name+".txt"
			pw_upfile="pathway_up_"+comparison_name+".txt"
			pw_dnfile="pathway_down_"+comparison_name+".txt"
		ogo_upfile=open(go_upfile,"w")
		ogo_dnfile=open(go_dnfile,"w")
		opw_upfile=open(pw_upfile,"w")
		opw_dnfile=open(pw_dnfile,"w")
		#############################
		#print 'de3...........................'
		demRNA_matrix,datalen=all_txt(filename)		
		#print 'de3...........................'
		new_demRNA_matrix=demRNA_matrix
		#print 'de3...........................'
		filenote.write("%s\t%s\t%s\t%s\t%s\n"%(comparison_name,len(demRNA_matrix)-1,fc,p,sta))
		#############################33
		colnum=len(new_demRNA_matrix[0])
		probename=new_demRNA_matrix[0][0]
		mRNAup,mRNAdown=sep_regulation_matrix(new_demRNA_matrix)
		#print 'de4...........................'
		#print 'de4...........................',new_demRNA_matrix[0]
		geneloc=new_demRNA_matrix[0].index(genesymbl)
		upsymbel=[];downsymbel=[]
		#print 'de3...........................'
		for i in mRNAup[1:]:
			upsymbel.append(i[geneloc])
		upsymbel=sorted(set([x for x in upsymbel if (x!="" and x!="N/A")]))
		ogo_upfile.write("\n".join(upsymbel))
		for i in mRNAdown[1:]:
			downsymbel.append(i[geneloc])
		downsymbel=sorted(set([x for x in downsymbel if (x!="" and x!="N/A")]))
		ogo_dnfile.write("\n".join(downsymbel))
		for i in upsymbel:
			opw_upfile.write("%s\t%s\n"%(i,"orange"))
		for i in downsymbel:
			opw_dnfile.write("%s\t%s\n"%(i,"yellow"))
		ogo_upfile.close()
		ogo_dnfile.close()
		opw_upfile.close()
		opw_dnfile.close()
		if repeatorder:
			if p:
				headname=comparison_name+"_fc_"+str(fc)+"_p_"+str(p)+"_"+sta
				sheetname=str(repeatorder)+"_"+comparison_name
			else:
				headname=comparison_name+"_fc_"+str(fc)
				sheetname=str(repeatorder)+comparison_name
		else:sheetname=headname=comparison_name
		#print new_demRNA_matrix[0]
		if p is not None:
			datalen=datalen-2
			annotation=new_demRNA_matrix[0][datalen*2+8:colnum]
			relationshipend=colnum-datalen*2-4
		else:
			annotation=new_demRNA_matrix[0][datalen*2+7:colnum]
			relationshipend=colnum-datalen*2
		###############################up
		print 'de2...........................'
		ws0 = wb.add_worksheet("up_"+sheetname[:28])
		ws0.set_column(0, 0, 28)
		rownum=len(mRNAup)
		#print rownum
		a=str("Fold Change cut-off: ") 
		b=str("P-value cut-off: ")
		c=str("Differentially expressed circRNAs for: ")
		d="Column A: "+str(probename)+', the ID of the identified circRNAs by DCC.'
		if p is not None:
			e="Column B: pvalue, pvalue between two groups of samples."
			#e1="Column C: FDR, FDR is calculated from Benjamini Hochberg FDR."
			f="Column C: fold change, the fold change between two groups of samples."
			g='Column D: regulation, "up" indicates up-regulation, and "down" means down-regulation.'			
			h="Column E ~ F: average junction reads of each group."
			i="Column G ~ H: averaged normalized Intensity of each group (log2 transformed)."
			j="Column I ~ "+column_index[datalen+8]+": junction reads of each sample."
			k="Column "+str(column_index[datalen+9])+" ~ "+str(column_index[datalen*2+8])+": Normalized Intensity of each sample (log2 transformed)."
			l="Column "+str(column_index[datalen*2+9])+" ~ "+str(column_index[datalen*2+12])+": the coordinates of circRNA."
			m="Column "+ str(column_index[datalen*2+13])+": circBaseID, the identifier of circBase (http://www.circbase.org)."
			n="Column "+ str(column_index[datalen*2+14])+": source, the source of the circRNA, including circBase, Guojunjie2014,..."
			o="Column "+ str(column_index[datalen*2+15])+": best_transcript, the best transcript of the circRNA."
			p1="Column "+ str(column_index[datalen*2+16])+": GeneName, the name of the circRNA-associated gene."
			q="Column "+ str(column_index[datalen*2+17])+": Catalog, the catalog of the circRNA, including exonic, intronic, ..."
			r="Column "+ str(column_index[datalen*2+18])+": predicted_sequence_length, the length of predicted circRNA sequence."
			if user_information_species.lower() == 'human':
				s="Column "+ str(column_index[datalen*2+19])+": circRNA-associated diseases sourced from (http://gyanxet-beta.com/circdb/)."
			else:
				s=""
		else:
			f="Column B: Fold change, Absolute Fold change between two samples. "
			g="Column C: Regulation, it depicts which sample has greater or lower intensity values wrt other sample."
			h="Column D ~ "+str(column_index[datalen*1+3])+": junction reads of each sample."
			i="Column "+str(column_index[datalen*1+4])+" ~ "+str(column_index[datalen*2+3])+": Normalized Intensity of each sample (log2 transformed)."
			l="Column "+str(column_index[datalen*2+4])+" ~ "+str(column_index[datalen*2+7])+": the coordinates of circRNA."
			m="Column "+ str(column_index[datalen*2+8])+": circBaseID, the identifier of circBase (http://www.circbase.org)."
			n="Column "+ str(column_index[datalen*2+9])+": source, the source of the circRNA, including circBase, Guojunjie2014,..."
			o="Column "+ str(column_index[datalen*2+10])+": best_transcript, the best transcript of the circRNA."
			p1="Column "+ str(column_index[datalen*2+11])+": GeneName, the name of the circRNA-associated gene."
			q="Column "+ str(column_index[datalen*2+12])+": Catalog, the catalog of the circRNA, including exonic, intronic, ..."
			r="Column "+ str(column_index[datalen*2+13])+": predicted_sequence_length, the length of predicted circRNA sequence."
			if user_information_species.lower() == 'human':
				s="Column "+ str(column_index[datalen*2+14])+": circRNA-associated diseases sourced from (http://gyanxet-beta.com/circdb/)."
			else:
				s=""
		#print 'de21...........................'
		font0 = wb.add_format({'font_name':'Times New Roman'})
		font0.set_font_size(10)
		font1 = wb.add_format({'font_name':'Times New Roman','bold':'true'})
		font1.set_font_size(10)
		red = wb.add_format({'font_name':'Times New Roman','bold':'true','font_color':'red'})
		red.set_font_size(10)
		#print 'de21a...........................'
		annotationformat =  wb.add_format({'align':'left', 'fg_color':'#FFFF99','font_name':'Times New Roman'})
		annotationformat.set_align('bottom')
		annotationformat.set_text_wrap()  # auto wrapping
		annotationformat1 =  wb.add_format({'align':'center', 'fg_color':'red','font_name':'Times New Roman','bold':'true'})
		annotationformat2 =  wb.add_format({'align':'center', 'fg_color':'#3366ff','font_name':'Times New Roman','bold':'true'})
		annotationformat3 =  wb.add_format({'align':'center', 'fg_color':'#008080','font_name':'Times New Roman','bold':'true'})		
		annotationformat4 =  wb.add_format({'align':'center', 'fg_color':'#800080','font_name':'Times New Roman','bold':'true'})
		annotationformat5 =  wb.add_format({'align':'center', 'fg_color':'#ff9900','font_name':'Times New Roman','bold':'true'})
		annotationformat6 =  wb.add_format({'align':'center', 'fg_color':'#00ccff','font_name':'Times New Roman','bold':'true'})
		annotationformat7 =  wb.add_format({'align':'center', 'fg_color':'green','font_name':'Times New Roman','bold':'true'})
		print 'de21b...........................'
		seg1=(a, font0);seg2=(b, font0);seg3=(c, font0);seg4=(d, font0);
		seg6=(f, font0);seg7=(g, font0);
		seg8=(h, font0);seg9=(i, font0);
		zhushi=16
		#print 'de22...........................'
		if p is not None:
			seg5=(e, font0);seg10=(j, font0);
			#seg51=(e1,font0);
			seg11=(k,font0);seg12=(l,font0);
			seg13=(m,font0);seg14=(n,font0);
			seg15=(o,font0);seg16=(p1,font0);
			seg17=(q,font0);seg18=(r,font0);
			seg19=(s,font0);
			zhushi=21
			ws0.merge_range(0,0,zhushi+mRNA_info_addline,  colnum-1,'',font0)
			ws0.write_rich_string(0,0,font0,c,red,'up_'+sheetname+'\n',font0,a,red,str(fc)+'\n',font0,b,red,str(p)+'\n'+'\n',font0,d+'\n',font0,e+'\n',font0,f+'\n',font0,g+'\n',font0,h+'\n',font0,i+'\n',font0,j+'\n',font0,k+'\n',font0,l+'\n',font0,m+'\n',font0,n+'\n',font0,o+'\n',font0,p1+'\n',font0,q+'\n',font0,r+'\n',font0,s+'\n',annotationformat)	
			#print "merge3..............."
		else:
			seg10=(l,font0);seg11=(m,font0);
			seg12=(n,font0);seg13=(o,font0);
			seg14=(p1,font0);seg15=(q,font0);
			seg16=(r,font0);seg17=(s,font0);
			ws0.merge_range(0,0,zhushi+mRNA_info_addline,  colnum-1,'',font0)
			#print 'de21c...........................'
			ws0.write_rich_string(0,0,font0,c,red,'up_'+sheetname+'\n',font0,a,red,str(fc)+'\n'+'\n',font0,d+'\n',font0,f+'\n',font0,g+'\n',font0,h+'\n',font0,i+'\n',font0,l+'\n',font0,m+'\n',font0,n+'\n',font0,o+'\n',font0,p1+'\n',font0,q+'\n',font0,r+'\n',font0,s+'\n',annotationformat)	
			#print 'de21c...........................'
		#############
		if repeatorder:
			if p is not None:
				zhushi=21
				#print "merge4..............."
				ws0.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' up regulated circRNAs',annotationformat1)
			else:
				#print "merge4..............."
				ws0.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' up regulated circRNAs',annotationformat1)
		else:
			if p is not None:
				#print "merge4..............."
				ws0.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' '+str(fc)+' up regulated circRNAs',annotationformat1)
			else:
				#print "merge4..............."
				ws0.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' '+str(fc)+' up regulated circRNAs',annotationformat1)
		if p is not None:
			zhushi=21
			#print "merge5..............."
			ws0.merge_range(zhushi+3+mRNA_info_addline,1,zhushi+3+mRNA_info_addline,  3,'P-value, Fold change and Regulation',annotationformat2)
			ws0.merge_range(zhushi+3+mRNA_info_addline,4,zhushi+3+mRNA_info_addline,  5,'Group--junction reads',annotationformat3)
			ws0.merge_range(zhushi+3+mRNA_info_addline,6,zhushi+3+mRNA_info_addline,  7,'Group--Normalized Intensity',annotationformat4)
			ws0.merge_range(zhushi+3+mRNA_info_addline,8,zhushi+3+mRNA_info_addline, 8+datalen-1,'Junction reads',annotationformat3)
			ws0.merge_range(zhushi+3+mRNA_info_addline,8+datalen,zhushi+3+mRNA_info_addline, 8+datalen*2-1,'Normalized Intensity',annotationformat5)
			ws0.merge_range(zhushi+3+mRNA_info_addline,8+datalen*2,zhushi+3+mRNA_info_addline,  colnum-1,'Annotations',annotationformat6)	
			#print "merge6..............."		
			for rows in range(rownum):
				for cols in range(colnum):
					if rows==0:
						#print "merge7..............."
						ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font1)
					else:
						#ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font0)
						ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font0)
			#print "merge7..............."
		else:
			ws0.merge_range(zhushi+3+mRNA_info_addline,1,zhushi+3+mRNA_info_addline,  2,'Fold change and Regulation',annotationformat2)
			ws0.merge_range(zhushi+3+mRNA_info_addline,3,zhushi+3+mRNA_info_addline,  datalen+2,'Junction reads',annotationformat3)
			ws0.merge_range(zhushi+3+mRNA_info_addline,datalen+3,zhushi+3+mRNA_info_addline,  datalen*2+2,'Normalized Intensity',annotationformat5)
			ws0.merge_range(zhushi+3+mRNA_info_addline,datalen*2+3,zhushi+3+mRNA_info_addline,  colnum-1,'Annotations',annotationformat6)	
			for rows in range(rownum):
				for cols in range(colnum):
					#ws0.col(cols).width = cell_width
					if rows==0:
						#ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font1)
						ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font1)
					else:
						#ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font0)
						ws0.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAup[rows][cols]), font0)
	##################down
		print 'de3...........................'
		ws1 = wb.add_worksheet("down_"+sheetname[:26])
		ws1.set_column(0, 0, 28)
		rownum=len(mRNAdown)
		if p is not None:
			seg5=(e, font0);
			#seg51=(e1,font0);
			seg11=(k,font0);seg12=(l,font0);
			seg13=(m,font0);seg14=(n,font0);
			seg15=(o,font0);seg16=(p1,font0);
			seg17=(q,font0);seg18=(r,font0);
			seg19=(s,font0);
			zhushi=21
			#print 'de4...........................'
			ws1.merge_range(0,0,zhushi+mRNA_info_addline,  colnum-1,'',font0)
			ws1.write_rich_string(0,0,font0,c,red,'down_'+sheetname+'\n',font0,a,red,str(fc)+'\n',font0,b,red,str(p)+'\n'+'\n',font0,d+'\n',font0,e+'\n',font0,f+'\n',font0,g+'\n',font0,h+'\n',font0,i+'\n',font0,j+'\n',font0,k+'\n',font0,l+'\n',font0,m+'\n',font0,n+'\n',font0,o+'\n',font0,p1+'\n',font0,q+'\n',font0,r+'\n',font0,s+'\n',annotationformat)	
		else:
			seg10=(l,font0);seg11=(m,font0);
			seg12=(n,font0);seg13=(o,font0);
			seg14=(p1,font0);seg15=(q,font0);
			seg16=(r,font0);seg17=(s,font0);
			ws1.merge_range(0,0,zhushi+mRNA_info_addline,  colnum-1,'',font0)
			ws1.write_rich_string(0,0,font0,c,red,'down_'+sheetname+'\n',font0,a,red,str(fc)+'\n'+'\n',font0,d+'\n',font0,f+'\n',font0,g+'\n',font0,h+'\n',font0,i+'\n',font0,l+'\n',font0,m+'\n',font0,n+'\n',font0,o+'\n',font0,p1+'\n',font0,q+'\n',font0,r+'\n',font0,s+'\n',annotationformat)	
		print 'de5...........................'
		if repeatorder:
			#print 'de6...........................'
			if p is not None:
				zhushi=21
				#print 'de6...........................'
				ws1.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' down regulated circRNAs',annotationformat7)
			else:
				ws1.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' down regulated circRNAs',annotationformat7)
		else:
			#print 'de7...........................'
			if p is not None:
				zhushi=21
				ws1.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' '+str(fc)+' fold down regulated circRNAs',annotationformat7)
			else:
				ws1.merge_range(zhushi+2+mRNA_info_addline,0,zhushi+2+mRNA_info_addline,  colnum-1,headname+' '+str(fc)+' fold down regulated circRNAs',annotationformat7)
		if p is not None:
			zhushi=21
			#print 'de8...........................'
			ws1.merge_range(zhushi+3+mRNA_info_addline,1,zhushi+3+mRNA_info_addline,  3,'P-value, Fold change and Regulation',annotationformat2)
			ws1.merge_range(zhushi+3+mRNA_info_addline,4,zhushi+3+mRNA_info_addline,  5,'Group--junction reads',annotationformat3)
			ws1.merge_range(zhushi+3+mRNA_info_addline,6,zhushi+3+mRNA_info_addline,  7,'Group--Normalized Intensity',annotationformat4)
			ws1.merge_range(zhushi+3+mRNA_info_addline,8,zhushi+3+mRNA_info_addline, 8+datalen-1,'junction reads',annotationformat3)
			ws1.merge_range(zhushi+3+mRNA_info_addline,8+datalen,zhushi+3+mRNA_info_addline, 8+datalen*2-1,'Normalized Intensity',annotationformat5)
			ws1.merge_range(zhushi+3+mRNA_info_addline,8+datalen*2,zhushi+3+mRNA_info_addline,  colnum-1,'Annotations',annotationformat6)			
			for rows in range(rownum):
				for cols in range(colnum):
					if rows==0:
						#print 'de9...........................'
						ws1.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAdown[rows][cols]), font1)
					else:
						#ws1.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAdown[rows][cols]), font0)
						ws1.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAdown[rows][cols]), font0)
		else:
			ws1.merge_range(zhushi+3+mRNA_info_addline,1,zhushi+3+mRNA_info_addline,  2,'Fold change and Regulation',annotationformat2)
			ws1.merge_range(zhushi+3+mRNA_info_addline,3,zhushi+3+mRNA_info_addline,  datalen+2,'junction reads',annotationformat3)
			ws1.merge_range(zhushi+3+mRNA_info_addline,datalen+3,zhushi+3+mRNA_info_addline,  datalen*2+2,'Normalized Intensity',annotationformat5)
			ws1.merge_range(zhushi+3+mRNA_info_addline,datalen*2+3,zhushi+3+mRNA_info_addline,  colnum-1,'Annotations',annotationformat6)			
			for rows in range(rownum):
				for cols in range(colnum):
					if rows==0:
						ws1.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAdown[rows][cols]), font1)
					else:
						#ws1.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAdown[rows][cols]), font0)
						ws1.write(rows+zhushi+4+mRNA_info_addline, cols, formats(mRNAdown[rows][cols]), font0)
		print 'de10...........................'
	#print 'de11...........................'
	filenote.close()
	#print 'de12...........................'
	wb.close()
	#print 'de13...........................'
	os.chdir(global_now_location)

def get_allde_excel():
	a= conclude_comparison_info()
	lena=len(a)
	newa=[]
	for i in range(lena):
		try:newa.append([])
		except:newa=[[]]
	count=0
	for i in a:
		for j in i:
			if "." in str(j) and not "_vs_" in str(j):
				j0=j.split("_")
				temp=[]
				for k in j0:
					if "." in str(k):
						temp.append("1.0")
					else:temp.append(str(k))
				j1="_".join(temp)
				newa[count].append(j1)
			else:newa[count].append(j)
		count=count+1
	print 'newa: ', newa
	demRNA_excel("allde_",newa,all_de_loc)
	try:
		shutil.move(all_de_loc+'/All Comparisons.xlsx',global_now_location)
	except:
		pass
	shutil.rmtree("all_de_folder")
###############################################################################################################venn figure and Excel
def parse_venn_comparison():
	venn_comparisons = cf.options("venn_fig_excel") 	
	print 'venn_fig_excel'
	print venn_comparisons
	RNA_A_vs_B_ud_FC_P_S = []
	if len(venn_comparisons)%3 == 0:
		print 'venn_fig_excel0'
		for i in range(len(venn_comparisons)/3):
			print 'venn_fig_excel1'
			var_fileform = remove_space(cf.get('venn_fig_excel', 'var_fileform_'+str(i+1))).split(',')      # [a1:xx, a2:xx, a3:xx]
			logicform = remove_space(cf.get('venn_fig_excel', 'logicform_'+str(i+1)))     # such as a1&a2&a3
			setname =  cf.get('venn_fig_excel', 'setname_'+str(i+1))
			RNA_A_vs_B_ud_FC_P_S.append([var_fileform,logicform,setname])
	return RNA_A_vs_B_ud_FC_P_S
	
def get_input_data(fileform,ele_new_venn_compare): 
	"""为了得到一些组合名字得到的文件名，tags（a1,a2）,以及与probes之间的关系矩阵，方便后续使用
	"""
#fileform:[{'a1':['all.LncRNA.4outof8_NMrelationship.txt.PD_vs_Normal.1.5.0.05.paired.fdr.txt','up']}，{'a2':['all.LncRNA.4outof8_NMrelationship.txt.PD_vs_Normal.2.0.0.05.unpaired.fdr.txt','up']}]
#ele_new_venn_compare:[{'a1':['lncRNA','PD_vs_Normal','up','1.5','0.05','paired']},{'a2':['lncRNA','PD_vs_Normal','up','2.0','0.05','unpaired']},a1&a2,'comparison11']
	dic_ele_new_venn_compare={}
	matrix_list=[]
	filetags=[]
	filenames_ud=[]
	filenames=[]
	dic_list_filename_matrix={}
	filelist1_dic={}
	dic_list_tagname_matrix={}
	for i in ele_new_venn_compare[:-2]:
		dic_ele_new_venn_compare[i.keys()[0]]="_".join([str(ele) for ele in i.values()[0]])##{'a1':lncRNA_PD_vs_Normal_up_1.5_0.05_paired','a2':'lncRNA_PD_vs_Normal_up_2.0_0.05_unpaired'}
	for i in fileform:
		filetags.append(i.keys()[0])#[a1,a2,]
		filenames_ud.append(i.values()[0])#[['all.LncRNA.4outof8_NMrelationship.txt.PD_vs_Normal.1.5.0.05.paired.fdr.txt','up'],['all.LncRNA.4outof8_NMrelationship.txt.PD_vs_Normal.2.0.0.05.unpaired.fdr.txt','up']]
	count=0
	for i in filenames_ud:
		de_matrix,datalen=all_txt(i[0])
		new_de_matrix=rearrange_de_matrix(de_matrix)
		if i[1]=="up":
			x_matrix=sep_regulation_matrix(new_de_matrix)[0]
		elif i[1]=="down":
			x_matrix=sep_regulation_matrix(new_de_matrix)[1]
		else:x_matrix=new_de_matrix
		matrix_list.append(x_matrix)
		global first_row_name
		first_row_name=x_matrix[0][0]
		x={}
		for j in x_matrix:
			if j[0]==x_matrix[0][0]:
				pass
			else:
				x[j[0]]=j[1:]
		dic_list_filename_matrix[dic_ele_new_venn_compare[filetags[count]]]=x.keys()#{'mRNA_PD_vs_Normal_up_2.0_0.05_unpaired':['probename1','probename2','probename3',]
		dic_list_tagname_matrix[filetags[count]]=x.keys()#{'a1':['probename1','probename2','probename3',]
		filenames.append(dic_ele_new_venn_compare[filetags[count]])
		count=count+1
	return dic_list_filename_matrix,dic_list_tagname_matrix,matrix_list,filetags,filenames
##########################R venn figure

def list_Rform(var,LIST):
	"""
	write_R_file函数的工具函数
	"""
	temp=[]
	temp.append(var)
	temp.append(" <- c(")
	count=0
	for i in LIST:
		count=count+1
		if count<len(LIST):
			a="'"+i+"',"
			temp.append(a)
		else:
			a="'"+i+"')"
			temp.append(a)
	return temp
def get_var_list(dic_list_filename_matrix):
	"""
	write_R_file函数的工具函数
	"""
	varlist=[]
	for i in range(len(dic_list_filename_matrix.keys())):
		varlist.append("k"+str(i))
	return varlist

def write_R_file(dic_list_filename_matrix):
	"""
	将集合信息写入R文件
	"""
	ofileR=open(fileR,"w")
	ofileR.write("library(VennDiagram)\n")
	count=0
	for k in dic_list_filename_matrix.values():
		L=list_Rform(varlist[count],k)
		count=count+1
		for j in L:
			ofileR.write("%s"%j)
		ofileR.write("\n")
	ofileR.write("venn.diagram(margin =%s,list("%venn_margin)
	count=0
	for i in dic_list_filename_matrix.keys():
		if count+1<len(dic_list_filename_matrix.keys()):
			ofileR.write('''"%s"=%s,'''%(i[5:],varlist[count]))
		else: ofileR.write('''"%s"=%s),fill=c('''%(i[5:],varlist[count]))
		count=count+1
	count=0
	for i in range(len(dic_list_filename_matrix.keys())):
		count=count+1
		if count<len(dic_list_filename_matrix.keys()):
			ofileR.write('''"%s",'''%collist[i])
		else:ofileR.write('''"%s"),'''%collist[i])
	ofileR.write('''"%s.tiff")'''%(fileR[:-2]))
	ofileR.close()
def data_len(List):
	count=0
	for i in List:
		if "raw" in i:
			count=count+1
	return count

def get_anno_rn(alldefilenamelist):
	#alldefilenamelist([alldefilename,ud])
	tempalldematrixlist=[]
	alldematrixlist=[]
	for i in alldefilenamelist:
		demRNA_matrix,datalen=all_txt(i[0])
		new_demRNA_matrix=rearrange_de_matrix(demRNA_matrix)
		tempalldematrixlist.append(new_demRNA_matrix)
	alldematrixlist=tempalldematrixlist
	return alldematrixlist

def segmatrix(matrix):
	reg_loc=matrix[0].index("Regulation")
	raw_all_indexes=[]
	norm_all_indexes=[]
	raw_indexes=[]
	norm_indexes=[]
	for x in matrix[0]:
		if x[-5:]=="(raw)":
			raw_all_indexes.append(matrix[0].index(x))
		if x[-12:]=="(normalized)":
			norm_all_indexes.append(matrix[0].index(x))
		if x[-5:]=="(raw)" and x[:6]!="group-":
			raw_indexes.append(matrix[0].index(x))
		if x[-12:]=="(normalized)" and x[:6]!="group-":
			norm_indexes.append(matrix[0].index(x))
	raw1_loc=raw_all_indexes[0]
	newmatrix=[]
	for i in matrix:
		temp=[]
		seqname=i[0]
		pfr=i[1:reg_loc+1]
		annotation=i[reg_loc+1:raw1_loc]
		raw=[i[m] for m in raw_indexes]
		norm=[i[m] for m in norm_indexes]
		temp.append(seqname)
		temp.append(pfr)
		temp.append(annotation)
		temp.append(raw)
		temp.append(norm)
		newmatrix.append(temp)
	return newmatrix

def matrix2dic(matrix):
	dic={}
	for i in matrix:
		try:dic[i[0]].append(i)
		except:
			try:dic[i[0]]=[i]
			except:dic={i[0]:[i]}
	return dic

def getAB_matrix(matrix,allmatrix,AB):#get matrix from AB and probename
	"""
	得到集合运算得到的AB(probes list)中每一个探针对应的信息以及header
	"""
	pfr_len=len(get_pfr(matrix[0]))
	newmatrix=segmatrix(matrix)
	newallmatrix=segmatrix(allmatrix)
	dicnewmatrix=matrix2dic(newmatrix)
	dicallmatrix=matrix2dic(newallmatrix)
	dekeys=dicnewmatrix.keys()
	alldekeys=dicallmatrix.keys()
	newAB=[first_row_name]+AB
	ele_matrix=[]
	count=0
	for i in newAB:
		if i in dekeys:
			count=count+1
			for m in dicnewmatrix[i]:
				ele_matrix.append(m)
		elif i in alldekeys:
			for m in dicallmatrix[i]:
				x=[m[0]]+["\t"*pfr_len]+m[2:]
				ele_matrix.append(x)
	return ele_matrix,count

def get_pfr(list):
	pfr=["P-value","FDR","Fold Change","Regulation"]
	pfr=["P-value","Fold Change","Regulation"]
	pfr_loc=[]
	for i in list:
		if i in pfr:
			pfr_loc.append(list.index(i))
	return pfr_loc

def get_wanted_index(theheader):
	temp=[]
	dic={}
	for i in theheader:
		for j in i:
			temp.append(j)
	count=0
	for x in temp:
		try:dic[x].append(count)
		except:
			try:dic[x]=[count]
			except:dic={x:[count]}
		count=count+1
	newtempindex=[]
	for k,v in dic.items():
		newtempindex.append(dic[k][0])
	return newtempindex

def uniq_data(list,index):
	list1=[]
	for i in list:
		for j in i:
			list1.append(j)
	list2=[list1[i] for i in index]
	return list2


def combine_matrix(matrix_list,alldefilenamelist):
	fileslen=len(matrix_list)
	alldematrixlist=get_anno_rn(alldefilenamelist)
	dematrix2=[]
	new_filenames_list=[]
	for i in range(fileslen):
		matrix=matrix_list[i]
		allmatrix=alldematrixlist[i]
		ele_matrix,count=getAB_matrix(matrix,allmatrix,AB)
		#print count
		if count!=1:
			dematrix2.append(ele_matrix)
			new_filenames_list.append(filenames[i])
	fileslen=len(new_filenames_list)
	#print len(dematrix2)
	datalen=len(dematrix2[0])
	dematrix3=[]
	for j in range(datalen):
		for i in range(fileslen):
			k=dematrix2[i][j]
			if i==0:
				seqname=k[0]
				pfr=[k[1]]
				anno=[k[2]]
				raw=[k[3]]
				norm=[k[4]]
			else:
				pfr=pfr+[k[1]]
				#anno=anno+[k[2]]
				raw=raw+[k[3]]
				norm=norm+[k[4]]
		temp=[seqname,pfr,anno,raw,norm]
		dematrix3.append(temp)
	header=dematrix3[0]
	raw_header=header[3]
	norm_header=header[4]
	raw_index=get_wanted_index(raw_header)
	norm_index=get_wanted_index(norm_header)
	dematrix4=[]
	for i in dematrix3:
		seqname=i[0]
		pfr=i[1]
		anno=i[2][0]
		raw=i[3]
		norm=i[4]
		uniqraw=uniq_data(raw,raw_index)
		uniqnorm=uniq_data(norm,norm_index)
		temp=[seqname,pfr,anno,uniqraw,uniqnorm]
		dematrix4.append(temp)
	return dematrix4,new_filenames_list

def set_to_excel():
	fileslen=len(new_filenames_list)
	print 'set to excel'
	last_matrix=[]
	for i in all_matrix1:
		seqname=i[0]
		pfr=i[1]
		temp=[]
		for j1 in pfr:
			for j2 in j1:
				temp.append(j2)
		anno=i[2]
		raw=i[3]
		norm=i[4]
		last_matrix.append([seqname]+temp+anno+raw+norm)
	reg_loc=[]
	raw_loc=[]
	norm_loc=[]
	count=0
	for j in last_matrix[0]:
		if j=="Regulation":
			reg_loc.append(count)
		elif j[-6:]=="](raw)":
			raw_loc.append(count)
		elif j[-13:]=="](normalized)":
			norm_loc.append(count)
		count=count+1
	rownum=len(last_matrix)
	colnum=len(last_matrix[0])
	print 'set to excel0'
	workbook=xlsxwriter.Workbook(fileR[:-2]+'.xlsx')
	workbook.use_zip64()
	print 'set to excel1'
	ws0 = workbook.add_worksheet(fileR[:-2])
	print 'set to excel2'
	font0 = workbook.add_format({'font_name':'Times New Roman'})
	font0.set_font_size(10)
	font1 = workbook.add_format({'font_name':'Times New Roman'})
	font1.set_bold('true')
	font1.set_font_size(10)
	annotationformat =  workbook.add_format({'align':'left', 'fg_color':'#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_align('bottom')
	annotationformat.set_text_wrap()  # auto wrapping
	annotationformat1 =  workbook.add_format({'align':'center', 'fg_color':'red','font_name':'Times New Roman','bold':'true'})
	annotationformat2 =  workbook.add_format({'align':'center', 'fg_color':'#3366ff','font_name':'Times New Roman','bold':'true'})
	annotationformat3 =  workbook.add_format({'align':'center', 'fg_color':'#008080','font_name':'Times New Roman','bold':'true'})		
	annotationformat4 =  workbook.add_format({'align':'center', 'fg_color':'#800080','font_name':'Times New Roman','bold':'true'})
	annotationformat5 =  workbook.add_format({'align':'center', 'fg_color':'#ff9900','font_name':'Times New Roman','bold':'true'})
	annotationformat6 =  workbook.add_format({'align':'center', 'fg_color':'#00ccff','font_name':'Times New Roman','bold':'true'})
	annotationformat7 =  workbook.add_format({'align':'center', 'fg_color':'green','font_name':'Times New Roman','bold':'true'})
	ws0.merge_range(0, 0, 0,colnum-1,newheaderlogic,annotationformat1)
	count=0
	print 'set to excel3'
	for i in new_filenames_list:
		if i !="":
			if count==0:
				ws0.merge_range(1, 1, 1,reg_loc[count], i[5:], annotationformat)
			elif count%2:
				ws0.merge_range(1,  reg_loc[count-1]+1,  1, reg_loc[count], i[5:], annotationformat4)
			else:ws0.merge_range(1,  reg_loc[count-1]+1, 1, reg_loc[count], i[5:], annotationformat5)
		count=count+1
	ws0.merge_range(1,  reg_loc[-1]+1, 1,   raw_loc[0]-1, "Annotations", annotationformat3)
	ws0.merge_range(1,  raw_loc[0],    1,   raw_loc[-1], "Raw Intensity", annotationformat1)
	ws0.merge_range(1,  norm_loc[0],   1,   norm_loc[-1], "Normalized Intensity", annotationformat2)
	print 'set to excel30'
	for rows in range(rownum):
		for cols in range(colnum):
			#ws0.col(cols).width = cell_width
			ws0.set_column(cols,cols, 25)
			if rows==0:
				ws0.write(rows+2, cols, formats(unicode(last_matrix[rows][cols],"utf-8")), font1)
			else:
				ws0.write(rows+2, cols, formats(unicode(last_matrix[rows][cols],"utf-8")), font0)
	print 'set to excel4'
	workbook.close()

def get_venn_compare():
	"""
	将parse_venn_comparison()得到的list进一步处理，fileform存入字典，
	然后跟logic form和setname一起放入列表new_venn_compare
	"""
	venn_compare = parse_venn_comparison()
	new_venn_compare=[]
	for i in venn_compare:
		count=0
		b=[]
		for j in i:
			if count==0:
				for k in j:
					a={}
					k1=k.split(":")[0]
					k2=k.split(":")[1].strip(";").split(";")
					a[k1]=k2
					b.append(a)		
			else:
				b.append(j)
			count=count+1
		new_venn_compare.append(b)
	return new_venn_compare
	#1.lncRNA/mRNA 2.compare name 3.up/down 4.fc 5.p 6.sta
def headerlogic(logicform,filetags,filenames):
	fileslen=len(filetags)
	temp=logicform
	for i in range(fileslen):
		temp=temp.replace(filetags[i]," "+filenames[i][5:]+" ")
	return temp
def main_venn():
	global fileR,varlist,AB,filenames,all_matrix1,matrix_list,new_filenames_list,newheaderlogic
	new_venn_compare=get_venn_compare()
	for i in new_venn_compare:
		filenamelist=[]
		filetaglist=[]
		alldefilenamelist=[]
		for j in i[:-2]:
			filetag=j.keys()[0]
			temp=j.values()[0]
			if "up" in temp or "down" in temp:
				ud=temp[2]
				temp=temp[:2]+temp[3:]
			else:ud=None
			cc=conclude_comparison_info()
			#print temp ##为了检查命名而打印
			for k in cc:
				#print k ##为了检查命名而打印
				if set(temp[1:]).issubset(set(k[:-1])):
					comparison_name=k[0]
					fc=k[1]
					p=k[2]
					sta=k[3]
					if temp[0]=="mRNA":
						if p:
							filename='all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison_name+"."+str(fc)+"."+str(p)+"."+sta+".fdr.txt"
						else:
							filename='all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison_name+"."+str(fc)+".txt"
						if p:
							alldefilename='allde_all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison_name+".1.0.1.0."+sta+".fdr.txt"
						else:
							alldefilename='allde_all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison_name+".1.0.txt"
					filenamelist.append([filename,ud])
					filetaglist.append(filetag)
					alldefilenamelist.append([alldefilename,ud])
		logicform=i[-2]
		fileR=i[-1]+".R"
		fileform=[]
		count=0
		for m in filetaglist:
			temp_fileform={}
			temp_fileform[m]=filenamelist[count]
			fileform.append(temp_fileform)
			count=count+1
		dic_list_filename_matrix,dic_list_tagname_matrix,matrix_list,filetags,filenames=get_input_data(fileform,i)
		#print logicform,filetags,filenames
		newheaderlogic=headerlogic(logicform,filetags,filenames)
		temp_tag=dic_list_tagname_matrix.keys()
		temp_name=dic_list_tagname_matrix.values()
		p=re.compile(r'[\W]*')
		m=re.split(p,logicform)
		#print m
		for n in m:
			if n!="":
				newname="tempdata["+str(temp_tag.index(n))+"]"
				logicform=logicform.replace(n,"set("+newname+")")
		order="sorted("+str(logicform)+")"
		#print order
		tempdata=[]
		for p in range(len(temp_tag)):
			t=dic_list_tagname_matrix.values()[p]
			tempdata.append(t)
		AB=eval(order)
		varlist=get_var_list(dic_list_filename_matrix)
		write_R_file(dic_list_filename_matrix)
		try:r.source(fileR)
		except:pass
		if len(AB)!=0:
			all_matrix1,new_filenames_list=combine_matrix(matrix_list,alldefilenamelist)
			set_to_excel()
		else:
			printWarning(newheaderlogic+" waring!!"+u' 集合运算结果是0，可以生成图，但不能生成venn excel，请调整运算表达式')
			#print 
#############################all header	
algn00 = xlwt.Alignment()
algn00.wrap = 1
style00 = xlwt.XFStyle()
style00.alignment = algn00
font00=Font()
font00.name='Times New Roman'
font00.bold=True
font00.colour_index='blue'
style00.font=font00

#style0.font.name = 'Times New Roman'
# style0.font.bold = True
#style0.font.height = 10*0x14 # height
pat00 = Pattern()
pat00.pattern = Pattern.SOLID_PATTERN
pat00.pattern_fore_colour = 0x2b # background color
style00.pattern = pat00

############## de header
algn0 = xlwt.Alignment()
algn0.wrap = 1
algn0.vert=xlwt.Alignment.VERT_TOP
style0 = xlwt.XFStyle()
style0.alignment = algn0
font0=Font()
font0.name='Times New Roman'
font0.bold=True
font0.colour_index='red'
style0.font=font0

#style0.font.name = 'Times New Roman'
# style0.font.bold = True
#style0.font.height = 10*0x14 # height
pat0 = Pattern()
pat0.pattern = Pattern.SOLID_PATTERN
pat0.pattern_fore_colour = 0x2b # background color
style0.pattern = pat0
############### up style
style1 = XFStyle()
style1.font.name = 'Times New Roman'
style1.font.bold = True
pat1 = Pattern()
pat1.pattern = Pattern.SOLID_PATTERN
pat1.pattern_fore_colour = 'red'
style1.pattern = pat1
# center 
al = Alignment()
al.horz = Alignment.HORZ_CENTER
style1.alignment = al

################# down title style
style2 = XFStyle()
style2.font.name = 'Times New Roman'
style2.font.bold = True
pat2 = Pattern()
pat2.pattern = Pattern.SOLID_PATTERN
pat2.pattern_fore_colour = 'green'
style2.pattern = pat2
a2 = Alignment()
a2.horz = Alignment.HORZ_CENTER
style2.alignment = a2
##########################
style20 = XFStyle()
style20.font.name = 'Times New Roman'
style20.font.bold = True
pat20 = Pattern()
pat20.pattern = Pattern.SOLID_PATTERN
pat20.pattern_fore_colour = 'orange'
style20.pattern = pat20
a20 = Alignment()
a20.horz = Alignment.HORZ_CENTER
style20.alignment = a20

##############title merge line
style3 = XFStyle()
style3.font.name = 'Times New Roman'
style3.font.bold = True
pat3 = Pattern()
pat3.pattern = Pattern.SOLID_PATTERN
pat3.pattern_fore_colour = 0x34
style3.pattern = pat3
a3 = Alignment()
a3.horz = Alignment.HORZ_CENTER
style3.alignment = a3

style4 = XFStyle()
style4.font.name = 'Times New Roman'
style4.font.bold = True
pat4 = Pattern()
pat4.pattern = Pattern.SOLID_PATTERN
pat4.pattern_fore_colour = 0x28
style4.pattern = pat4
a4 = Alignment()
a4.horz = Alignment.HORZ_CENTER
style4.alignment = a4

style5 = XFStyle()
style5.font.name = 'Times New Roman'
style5.font.bold = True
pat5 = Pattern()
pat5.pattern = Pattern.SOLID_PATTERN
pat5.pattern_fore_colour = 0x26
style5.pattern = pat5
a5 = Alignment()
a5.horz = Alignment.HORZ_CENTER
style5.alignment = a5

#content style
style6 = XFStyle()
style6.font.name = 'Times New Roman'
left = Alignment()
left.horz = Alignment.HORZ_LEFT
style6.alignment = left

style7 = XFStyle()
style7.font.name = 'Times New Roman'
style7.font.bold = True
pat7 = Pattern()
pat7.pattern = Pattern.SOLID_PATTERN
pat7.pattern_fore_colour = 0x24
style7.pattern = pat7
a7 = Alignment()
a7.horz = Alignment.HORZ_CENTER
style7.alignment = a7

style8 = XFStyle()
style8.font.name = 'Times New Roman'
style8.font.bold = True
#pat8 = Pattern()
#pat8.pattern = Pattern.SOLID_PATTERN
#pat8.pattern_fore_colour = 0x22
#style8.pattern = pat8
# center 
a8 = Alignment()
a8.horz = Alignment.HORZ_CENTER
style8.alignment = a8

style9 = XFStyle()
style9.font.name = 'Times New Roman'
style9.font.bold = True
pat9 = Pattern()
pat9.pattern = Pattern.SOLID_PATTERN
pat9.pattern_fore_colour = 0x30
style9.pattern = pat9
# center 
a9 = Alignment()
a9.horz = Alignment.HORZ_CENTER
style9.alignment = a9
##############################

####### 绘图

# box plot
def plot_box(file,angle):
    # box plot
    # parameter of box
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['patch.facecolor'] = '#00FFFF'
    for i in range(int(user_information_sample_number)):
        try:treatments.append([])
        except:treatments=[[]]
    box_matrix,datalen=all_txt(file)
    matrix_len=len(box_matrix)
    names = [i.replace('](normalized)','').replace('[','') for i in box_matrix[0][1+int(user_information_sample_number):int(user_information_sample_number)*2+1]]
    for i in range(1,matrix_len):
    	temp=box_matrix[i][(1+int(user_information_sample_number)):(int(user_information_sample_number)*2+1)]
    	for j in range(datalen):
    		treatments[j].append(float(temp[j]))
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    bp = ax.boxplot(treatments, sym='r+',patch_artist=True, bootstrap=5000)
    text_transform= mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    setp(bp['whiskers'], color='#3030FF',  linestyle='-' )
    setp(bp['fliers'], markersize=5.0)
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9,bottom=0.3)##0.25 to 0.3
    xticks(range(1,int(user_information_sample_number)+1),names,rotation=float(angle))
    xlabel('sample',fontsize=16)
    ylabel('Normalized Intensity Values',fontsize=16)
    grid(linestyle='solid',lw=0.5, color='#C3C3C3')
    savefig(str(angle)+'_Box.png',dpi=301)
    plt.close('all')
# scatter plot
def get_limit(ylist,xlist):
	ymin=int(math.floor(min(ylist)))
	ymax=int(math.floor(max(ylist)+1))
	xmin=int(math.floor(min(xlist)))
	xmax=int(math.floor(max(xlist)+1))
	xymin=min(ymin,xmin)
	xymax=max(ymax,xmax)
	if xymax-float(max(xlist))<=0.2 or xymax-float(max(ylist))<=0.2:
		xymax=xymax+1
	elif xymax-float(max(xlist))>=1.2 and xymax-float(max(ylist))>=1.2:
		xymax=xymax-1
	else:pass
	if xymin-float(min(xlist))>=-0.2 or xymin-float(min(ylist))>=-0.2:
		xymin=xymin-1
	else:pass
	return xymin,xymax
def two_sample_scatter_plot(file,your_comparisons):
    """
    scatter plot for two sample comparison, A - col 3, B - col 4
    all.mRNA.2outof3.txt.B_vs_A
    
    """
    
    comparisons = your_comparisons   # [A, B]
    for comparison in comparisons:
        if file.startswith('all.mRNA') or file.startswith('scatter_all.mRNA'): # mRNA
            if file=='all.mRNA':
                f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison[0]+'_vs_'+comparison[1]
                distance = log2(float(comparison[2]))
            if file=='scatter_all.mRNA':
                f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison[1]+'_vs_'+comparison[3]
                distance = log2(float(comparison[4]))
            scatter_matrix,data_len=all_txt(f)
            matrix_len=len(scatter_matrix)
            a_name=scatter_matrix[0][3]
            b_name=scatter_matrix[0][4]

            filetemp=open("tempscatter.txt","w")
            for i in range(matrix_len):
    	        filetemp.write("%s\t%s\t%s\n"%(scatter_matrix[i][0],scatter_matrix[i][3],scatter_matrix[i][4]))
            filetemp.close()
            f1="tempscatter.txt"
            #########
            fig = plt.figure(figsize=(8,8))
            #t=a
            rcode = 's <- read.table("%s", header = TRUE, sep = "\t",row.names=1,fill = TRUE)' %f1
            robjects.r(rcode)
            robjects.r('''dat <-s[,1:2]''')
            robjects.r('''rbPal <- colorRampPalette(c("blue","blue","yellow","orange","darkorange","orangered1","red"))''')
            #robjects.r('''rbPal <- colorRampPalette(c("blue","blue","yellow","yellow","orange","darkorange","orangered1","red"))''')
            robjects.r('''dat$Col <- rbPal(100)[as.numeric(cut(dat[,1],breaks = 100))]''')
            t= list(r('dat[,3]'))
            a= list(r('dat[,1]'))
            b= list(r('dat[,2]'))
            os.remove(f1)
            plt.scatter(b, a, marker='s',c=t,edgecolor='gray',linewidth='0.25')
            
            ylabel(a_name)
            xlabel(b_name)
            left,right=get_limit(b,a)
            plot([left,right],[left,right], c='#63C363', lw='1.5') # mid
            plot([left,right],[left+distance,right+distance], c='#63C363', lw='1.5') # up
            plot([left,right],[left-distance,right-distance], c='#63C363', lw='1.5') # up
            xlim(left,right)
            ylim(left,right)
            grid(linestyle='solid',lw=1.5, color='#C3C3C3')
            if file=='all.mRNA' and comparison[-1]==None:
                savefig('scatter_'+comparison[0]+'_vs_'+comparison[1]+'.png', dpi=301)
            elif file=='all.mRNA' and (comparison[-1] is not None):
                savefig('scatter_'+comparison[0]+'_vs_'+comparison[1]+'.'+comparison[2]+'.png', dpi=301)
            elif file=='scatter_all.mRNA'and comparison[-1] == None:
                savefig('scatter_'+comparison[1]+'_vs_'+comparison[3]+'.png', dpi=301)
            elif file=='scatter_all.mRNA'and (comparison[-1] is not None):
                savefig('scatter_'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.png', dpi=301)
            plt.close(fig)
        plt.close('all')
#############################
 

############ 两组的scatter图
def fun_2groups_scatter(file,comparison):
    if "all.mRNA" in file:
        f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison[1]+'_vs_'+comparison[3]+'.tmp'
    else:f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'_NMrelationship.txt.'+comparison[1]+'_vs_'+comparison[3]+'.tmp'
    scatter_matrix,data_len=all_txt(f)
    matrix_len=len(scatter_matrix)
    a_name=scatter_matrix[0][6]
    b_name=scatter_matrix[0][7]

    filetemp=open("tempscatter.txt","w")
    for i in range(matrix_len):
    	filetemp.write("%s\t%s\t%s\n"%(scatter_matrix[i][0],scatter_matrix[i][6],scatter_matrix[i][7]))
    filetemp.close()
    f1="tempscatter.txt"
    ########
    fig = plt.figure(figsize=(8,8))

    rcode = 's <- read.table("%s", header = TRUE, sep = "\t",row.names=1,fill = TRUE)' %f1
    robjects.r(rcode)
    robjects.r('''dat <-s[,1:2]''')
    robjects.r('''rbPal <- colorRampPalette(c("blue","blue","yellow","orange","darkorange","orangered1","red"))''')
    #robjects.r('''rbPal <- colorRampPalette(c("blue","blue","yellow","yellow","orange","darkorange","orangered1","red"))''')
    robjects.r('''dat$Col <- rbPal(100)[as.numeric(cut(dat[,1],breaks = 100))]''')
    t= list(r('dat[,3]'))
    a= list(r('dat[,1]'))
    b= list(r('dat[,2]'))
    os.remove(f1)
    plt.scatter(b, a, marker='s',c=t,edgecolor='gray',linewidth='0.25')
    ylabel(a_name)
    xlabel(b_name)
    left,right=get_limit(b,a)
    plot([left,right],[left,right], c='#63C363', lw='1.5') # mid
    distance = log2(float(comparison[4]))
    plot([left,right],[left+distance,right+distance], c='#63C363', lw='1.5') # up
    plot([left,right],[left-distance,right-distance], c='#63C363', lw='1.5') # up
    xlim(left,right)
    ylim(left,right)
    grid(linestyle='solid',lw=1.5, color='#C3C3C3')
    if "all.mRNA" in file and comparison[-1]==None:
        savefig('scatter_'+comparison[1]+'_vs_'+comparison[3]+'.png', dpi=301)
    elif "all.mRNA" in file and (comparison[-1] is not None):
        savefig('scatter_'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.png', dpi=301)
    #else:savefig(file+'.'+XoutY_X+'outof'+user_information_sample_number+'_NMrelationship.txt.'+comparison[1]+'_vs_'+comparison[3]+'.'+comparison[4]+'.'+comparison[5]+'.scatter.png', dpi=301)
    plt.close(fig)

def two_group_scatter_plot(file):
	"""
	scatter plot for two group comparison,   6, 7 (0-based)
	all.mRNA.4outof8.txt.PD_vs_Normal.tmp 
	"""
	comparison_paired = ggpn #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
	comparison_unpaired = ggupn #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
	comparisons=comparison_paired+comparison_unpaired
	tempdic={}
	for x in range(len(comparisons)):
		if comparisons[x][-1] is not None:
			tag=comparisons[x][-1].split("_")[-1]
			judgelist=[comparisons[x][1],comparisons[x][3],comparisons[x][4]]
			judgestr="__".join(judgelist)
			try:tempdic[judgestr].append(x)
			except:tempdic={judgestr:[x]}
	loclist=[]
	ggpn_temp=[]
	ggupn_temp=[]
	for k,v in tempdic.items():
		if len(v) >1:
			loclist.append(v)
	#print loclist
	for i in loclist:
		for j in range(len(i)):
			if j==0:
				comparisons[i[0]][-1]=None
			else:
				comparisons[i[j]]=[]
	for i in range(len(comparison_paired)):
		if len(comparisons[i])!=0:
			ggpn_temp.append(comparisons[i])
	for i in range(len(comparison_paired),(len(comparison_unpaired)+len(comparison_paired))):
		if len(comparisons[i])!=0:
			ggupn_temp.append(comparisons[i])
	comparison_paired=ggpn_temp
	comparison_unpaired = ggupn_temp
	if file.startswith('unpaired_all.mRNA'):
		for comparison in comparison_unpaired:
			fun_2groups_scatter(file,comparison)
			plt.close('all')
	else:
		for comparison in comparison_paired:
			fun_2groups_scatter(file,comparison)
			plt.close('all')

mpl.rcParams['axes.axisbelow'] = True
def floats(da, listss):
    a= []
    for i in listss:
        a.append(float(da[i]))
    return a
def volcano_data(inputfile, Pcutoff, FCcutoff):
    """
    volcano plot drawing.
    inputfile   ==>     contains three colums:p-value (column 1), FCAbsolute (column 2) and regulation (column 3)
    Pcutoff     ==>     p-value cutoff, such as 0.05
    FCcutoff    ==>     fold change cutoff, such as 2.0
    """
    vol_matrix,datalen=all_txt(inputfile)
    matrix_len=len(vol_matrix)
    f=vol_matrix[1:]
    pvalues = []
    fcs = []
    cols = []
    ecolor=[]
    for linelist in f:
        if linelist[1]=='nan':
        	continue
        #print linelist[-1]
        if linelist[3] == 'down':
            fc = -numpy.log2(float(linelist[2]))
        else:
            fc = numpy.log2(float(linelist[2]))
        #print fc
        p=float(linelist[1])
        if p<=1E-12:
        	p=1E-12
        fcs.append(fc)
        pvalues.append(-numpy.log10(p))
        #print -numpy.log10(p)
        if fc >= numpy.log2(FCcutoff) and p <= Pcutoff:
            cols.append("#FF0000")
            ecolor.append("#B30000")
        elif fc <= -numpy.log2(FCcutoff) and p <= Pcutoff:
            cols.append('#FF0000')
            ecolor.append("#B30000")
        else:
            cols.append('#C2C2C2')
            ecolor.append("#8E8E8E")
    #f.close()
    return fcs, pvalues, cols, ecolor
def volcano_plot(inputfile, Pcutoff, FCcutoff,comparison,tag):
    fc, pv, co, eco= volcano_data(inputfile, Pcutoff, FCcutoff)
    fig = plt.figure(figsize=(6,8))
    grid(linestyle='solid',lw=1.5, color='#C3C3C3')
    scatter(fc, pv, marker='s',linewidths=0.75,color=co, s=100, edgecolors=eco)
    m = round(max(abs(min(fc)), max(fc)))
    #print m
    ym = max(pv)
    ymin = min(pv)
    plot([-int(m)-10, int(m)+10], [-numpy.log10(Pcutoff), -numpy.log10(Pcutoff)], c='#63C363', lw='1.5')
    plot([-numpy.log2(FCcutoff), -numpy.log2(FCcutoff)], [ymin-3, ym+3],  c='#63C363', lw='1.5')
    plot([numpy.log2(FCcutoff),   numpy.log2(FCcutoff)], [ymin-3, ym+3],  c='#63C363', lw='1.5')
    xlabel('log2(Fold Change)\n%s vs %s' % (comparison[1], comparison[3]))
    ylabel('-log10(pvalue)')
    xlim(-int(m)-0.5, int(m)+0.5)
    ylim(-0.2,ym+1) 
    if comparison[-1]==None:
        savefig('volcano_plot_'+comparison[1]+"_vs_"+comparison[3]+'.png',dpi=301)
    elif comparison[-1] is not None:
        savefig('volcano_plot_'+comparison[1]+"_vs_"+comparison[3]+"_"+str(FCcutoff)+'.'+str(Pcutoff)+'.'+tag+'.png',dpi=301)
    plt.close('all')
def two_group_volcano_plot(file):
    comparison_paired = ggpn #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
    comparison_unpaired = ggupn #  [['A1', 'A2', 'A3'],'A',['B1','B2','B3'], 'B','2.0', '0.05'], ...]
    if file.startswith('unpaired_all.mRNA'): 
        for comparison in comparison_unpaired:
            if 'all.mRNA' in file:
                f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison[1]+'_vs_'+comparison[3]+'.tmp'
                volcano_plot(f, float(comparison[5]), float(comparison[4]),comparison,"unpaired")
                plt.close('all')
            #else:
                #f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'_NMrelationship.txt.'+comparison[1]+'_vs_'+comparison[3]+'.tmp'
                #volcano_plot(f, float(comparison[5]), float(comparison[4]),comparison)
                #plt.close('all')
    else:
    	for comparison in comparison_paired:
            if 'all.mRNA' in file:
                f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'.txt.'+comparison[1]+'_vs_'+comparison[3]+'.tmp'
                volcano_plot(f, float(comparison[5]), float(comparison[4]),comparison,"paried")
                plt.close('all')
            #else:
                #f = file+'.'+XoutY_X+'outof'+user_information_sample_number+'_NMrelationship.txt.'+comparison[1]+'_vs_'+comparison[3]+'.tmp'
                #volcano_plot(f, float(comparison[5]), float(comparison[4]),comparison)
                #plt.close('all')

#two_group_volcano_plot('all.mRNA')
def compress_heatmap_data():
	for i in glob.glob("*.txt"):
		if i.startswith("pre_heatmap"):
			ofile=open(i,"rU")
			newfile=open(i[4:],"w")
			line=ofile.readline()
			newfile.write("%s"%line)
			dic={}
			while line:
				line=ofile.readline()
				if line=="":
					break
				sline=re.split("\t",line)
				dic[sline[0]]=line
			for k,v in dic.items():
				newfile.write("%s"%v)
			ofile.close()
			newfile.close()
def use_heatmap():
	nowloc=os.getcwd()
	cc= conclude_comparison_info()
	heatmaplist=[]
	heatmap_temp=[]
	for i in glob.glob("*.txt"):
		if i.startswith("heatmap_"):
			shutil.copy(i,heatmap_loc)
	os.chdir(heatmap_loc)
	for i in glob.glob("*.txt"):
		if i=="heatmap_all.mRNA."+outof+".txt":
			os.rename("heatmap_all.mRNA."+outof+".txt","Heatmap_of_All_Targets_Value.txt")
		elif i.startswith('heatmap_group_type') and i.endswith(('all.mRNA.'+outof+".txt")):
			j="Heatmap"+i[7:].replace('all.mRNA.'+outof+".txt","")+"_of_All_Targets_Value.txt"
			os.rename(i,j)
		else:
			#heatmap_all.LncRNA.3outof6_NMrelationship.txt.1-1_vs_1-2.2.0.tif
			i1=i.replace("heatmap_all.mRNA.",'').replace(outof+".txt.","").replace(".txt","")
			print 'i1',i1
			i2=i1.split(".")
			if "paired.fdr" in i1 or "unpaired.fdr" in i1:
				name=".".join(i2[:-6])
				print 'name',name
				#name='heatmap.'+name
				heatmap_temp.append([name,i2[-6]+"."+i2[-5],i2[-4]+"."+i2[-3],i2[-2]])# vs , fc, p ,sta
			else:
				name=".".join(i2[:-2])
				heatmap_temp.append([name,i2[-2]+"."+i2[-1]])# vs , fc, p ,sta
			heatmaplist.append(i)
	for i in range(len(heatmap_temp)):
		for j in cc:
			thename=None
			if heatmap_temp[i]==j[:4] and (j[2] is not None) and (j[3] is not None):
				if j[-1] is not None:
					thename=".".join(j[:-1])
					thename='heatamp.'+j[0]
				else:
					thename=j[0]
					thename='heatamp.'+j[0]
			elif heatmap_temp[i]==j[:2] and j[2:4] ==[None,None]:
				if j[-1] is not None:
					thename=".".join(j[:2])
					thename='heatamp.'+j[0]
				else:
					thename='heatamp.'+j[0]
			print thename,heatmap_temp[i],heatmaplist[i]
			if thename is not None:
				try:os.rename(heatmaplist[i],thename+".txt")
				except:pass
	os.system('%s' % heatmap_code)
	os.chdir(nowloc)
	
def remove_and_classify():
	try:os.mkdir("go")
	except:pass
	try:os.mkdir("pw")
	except:pass
	for i in glob.glob("*.txt"):
		if i.startswith("go_"):
			shutil.move(i,"go")
		elif i.startswith("pathway_"):
			shutil.move(i,"pw")
		elif i==config_file or i=="all.mRNA" or i=="note.txt":
			pass
		else:os.remove(i)
	for i in glob.glob("*"):
		if i.endswith(".tmp"):
			os.remove(i)
		elif i.endswith(".py") or i.endswith(".xls") or i.endswith(".txt") or i.endswith(".png") or i.endswith(".R") or i.endswith(".tiff") or i=="all.mRNA":
			pass
		elif os.path.isdir(i):
			pass
		else:os.remove(i)
########### HTML格式化

#printError printWarning
def main():
    """Generating LncRNA Report"""
    ################## 2, all.mRNA.txt.3out9, all.LncRNA.txt.3out9
    print u'1. '
    # mRNA
    seperate_mRNA_LncRNA()
    intensity_x_out_y('all.mRNA')
    print u'2. '
    # mRNA
    try:
        two_sample_comparison('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt')
        print u'2.1 '
    except:
        printError(u'2.1 ')
    # mRNA
    try:
        print u'2.2 '
        comparisonlist = parse_sample_group_comparison()
        one_sg_one_sg_comparison('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt',comparisonlist)        
    except:printError(u'2.2 ')
    try:
        print u'2.3 '
        comparisonlist = parse_group_sample_comparison()
        one_sg_one_sg_comparison('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt',comparisonlist)        
    except:printError(u'2.3 ')
    # mRNA
    try:
        get_two_group_unpaired_comparison('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt')
        print u'2.4 '  # all.mRNA.4outof8.txt.PD_vs_Normal.2.0.0.05.unpaired.fdr.txt        
    except:
        printError(u'2.4 ')
    # mRNA
    try:
        get_two_group_paired_comparison('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt')
        print u'2.5 '
    except:
        printError(u'2.5 ')
    
    
    ####################### 3, box plot of all.mRNA.3out9, all.LncRNA.3out9
    print u'3. '
     #mRNA
    try:
        print u'3.1 mRNA box'
        plot_box('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt',0)
        plot_box('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt',45)
        plot_box('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt',90)
    except:
        printError(u'3.1 mRNA')
    print u'4. '
    global ssn,sgn,gsn,ggpn,ggupn
    ssn,sgn,gsn,ggpn,ggupn=seperate_cc()
    sggsn=sgn+gsn
    try:
        two_sample_scatter_plot('all.mRNA',ssn)
        print u'4.1 '
    except:printError(u'4.1 ')
    try:
        two_sample_scatter_plot('scatter_all.mRNA',sggsn)    
        print u'4.2 '
    except:printError(u'4.2 ')

    try:
        two_group_scatter_plot('unpaired_all.mRNA')
        print u'4.3 mRNA '
    except:
        printError(u'4.3 ')
    try:
        two_group_scatter_plot('paired_all.mRNA')
        print u'4.4 mRNA '
    except:
        printError(u'4.4 ')
    ssn,sgn,gsn,ggpn,ggupn=seperate_cc()
    print u'5. '
    try:
        two_group_volcano_plot('unpaired_all.mRNA')
        print u'5.1 unpaired mRNA volcano_plot'
    except: printError(u'5.1 unpaired mRNA')
    try:
        two_group_volcano_plot('paired_all.mRNA')
        print u'5.2 paired mRNA volcano_plot'
    except:printError(u'5.2 paired mRNA volcano_plot')
    try:
        compress_heatmap_data()
        use_heatmap()
        print u'6. heatmap seems good'
    except:printError(u'6. heatmap failure')
    
    #############################################################
    #进行control type过滤
    ###############################################################
    print u'7. All Targets Value.xls Expression Profiling Data.xls'
    ##All targets value of mRNA excel.
    #try:
    print "all_txt function to all_mRNA_excel"
    print 'all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+".txt"
    os.chdir(global_now_location)
    all_mRNA_matrix,datalen=all_txt('all.mRNA.'+XoutY_X+'outof'+user_information_sample_number+'.txt')
    print "all.txt"
    all_mRNA_excel(all_mRNA_matrix,datalen)
    #except: printError(u"7.1 All targets value to mRNA Excel failure!")
    
    ######################################################DE value of mRNA and lncRNA to excel.
    print u'8. Be ready for making DE Excel!'
    try:
        cc= conclude_comparison_info()
        demRNA_excel("",cc,global_now_location)
        print "demRNA_excel completed"
        get_allde_excel()
    except: printError(u'8. DE mRNA Excel failure!')
    try:
        if parse_venn_comparison():
            print u'9. venn excel'
            main_venn()
            print u'9.1 venn figure_excel'
        else:pass
    except:printError(u'9.1 Venn figure or Venn Excel failed!')

main()


