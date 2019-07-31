# -*- coding: utf-8 -*-

"""
This is used for RNA-Seq data analysis all-in-one (DOing) with reference genomes to get the expression profiling of all genes (genes, lncRNAs, ...)

following the hisat2 - cuffquant -- cuffdiff pipeline

some files to prepare
1, UCSC genome  -- downloaded from genome.ucsc.edu
2, bowtie index of [1]
3, gtf file
4, config.txt file

"""

import ConfigParser
import hashlib
import threading
import re
import glob
import numpy
import os
import sys
import shutil
import os.path
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from xlsxwriter import *
from scipy import stats
import operator
from operator import itemgetter
import random
import math
from matplotlib.path import Path
import matplotlib.patches as patches
from PIL import Image

def png2tiff(file):
	pngfile = Image.open(file)
	# Save as TIFF
	pngfile.save(file.replace('.png','.tiff'), dpi=(600,600))
	os.remove(file)

mpl.rcParams['axes.axisbelow'] = True
from matplotlib.transforms import Affine2D


reload(sys)
sys.setdefaultencoding('utf-8')


#### excel column index
info_addline = 3
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
colalphabet = list(alphabet)+ ['A'+i for i in alphabet]+['B'+i for i in alphabet]+['C'+i for i in alphabet]+['D'+i for i in alphabet]+['E'+i for i in alphabet]+['F'+i for i in alphabet]+['G'+i for i in alphabet]+['H'+i for i in alphabet]+['I'+i for i in alphabet]+['J'+i for i in alphabet]+['K'+i for i in alphabet]+['L'+i for i in alphabet]+['M'+i for i in alphabet]+['N'+i for i in alphabet]+['O'+i for i in alphabet]+['P'+i for i in alphabet]+['Q'+i for i in alphabet]+['R'+i for i in alphabet]+['S'+i for i in alphabet]+['T'+i for i in alphabet]+['U'+i for i in alphabet]+['V'+i for i in alphabet]+['W'+i for i in alphabet]+['X'+i for i in alphabet]+['Y'+i for i in alphabet]+['Z'+i for i in alphabet]
#### excel column index


reload(sys)
sys.setdefaultencoding('utf-8')

####################CONSTANTS, no use here
pathway_colors = {'up':'orange', 'down': 'yellow'}
lncRNAclasses = ['3prime_overlapping_ncrna','antisense','lincRNA','retained_intron','sense_intronic','processed_transcript','misc_RNA', 'long_noncoding']

####################CONSTANTS
global_now_location=os.getcwd()
config_file = 'cloud3-RNA-seq-paired-end.human.mouse.rat.20190618.txt'  # DO NOT CHANGE FILE NAME!!!
cf = ConfigParser.ConfigParser()
cf.read(config_file)

# configparse

# adaptor
three_prime_adaptor_1 = cf.get('adaptor_information', 'three_prime_adaptor_1')		#  GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
three_prime_adaptor_2 = cf.get('adaptor_information', 'three_prime_adaptor_2')		#  GATCGGAAGAGCACACGTCTGAACTCCAGTCAC


# plot 
countthreshold = int(cf.get('plot_length_cutoff', 'expresscount'))
minlength = int(cf.get('plot_length_cutoff', 'minlength'))  
maxlength = int(cf.get('plot_length_cutoff', 'maxlength'))

# factor
FACTOR = float(cf.get('FACTOR','FACTOR')) # 0.1
lnc_associated_genefile = cf.get('mapping_parameters','gtf_nearby_gene')
geneinfofile = cf.get('mapping_parameters','geneinfofile')
# configparse over

# zero-step, check all files are there

# 0-1, check bowtie2 index

############################# First Step: create results folder (if exist, delete and recreate)

def create_folder(foldername):
	"""
	create a report folder, if exists, delete and recreate
	"""
	os.chdir(global_now_location)
	os.chdir('..')
	try:
		shutil.rmtree(foldername)
		os.mkdir(foldername)
	except:
		os.mkdir(foldername)
	os.chdir(global_now_location)

def create_all_folders():
	"""
	create all report folders, if exists, delete and recreate
	"""
	for folder in ['Alignment','Expression_Profiling','GEO_Uploading', 'Differentially_Expressed_RNAs', 'Sample_QC', 'GO_Analysis_Results', 'Pathway_Analysis_Results']:
		create_folder(folder)

def create_readme():
	"""
	create readme.txt
	"""
	os.chdir(global_now_location)
	os.chdir('..')
	try:
		os.remove('readme.txt')
	except:
		f=open('readme.txt','w')
		print >>f, """#COMPUTER SYSTEM REQUIREMENT
Because of the large raw data files, we recommend that your computer has at least 4GB RAM.

For Windows user
You can unzip the *.gz files using WINRAR (http://www.rarlab.com/), and open unzipped files using UltraEdit (http://www.ultraedit.com/index.html) or other TXT Editor software.

For Linux user
You can unzip this *.gz files using command line "gunzip *.gz", and open unzipped files using Emacs or other TXT Editor software.

----------------------------------Folders----------------------------------

Alignment	 			--	Mapped reads in BAM format
express				--	Expressed RNAs (excel format)
express_diff			--	Differentially Expressed RNA (excel format)
GO Analysis Report			--	Gene Ontology Analysis Results for Differentially Expressed protein-coding RNAs
Pathway Analysis Report		--	KEGG Pathway Analysis Results for Differentially Expressed protein-coding RNAs
Raw Data	 			--	Base calling data after Solexa Chastity quality filtering (FASTQ format)
Sample_QC_Report.pdf		--	Quanlity Control Report"""
	os.chdir(global_now_location)

def folder_main():
	create_all_folders()
	create_readme()

################################# 1, over

################################## 2, read config file and get raw fastq file

def getfastq():
	samples  = cf.options("sample_names")
	sample_names = []
	for sample in samples:
		name = cf.get('sample_names', sample)
		sample_names.append(name.strip())
	print(sample_names)
	return sample_names

############################## 2, over
 
############################# 3, FASTQC report (optional)

def _fastqc(fastqfile):
	"""
	must be fastq file, output is curdir
	"""
	fastqc_cmd = '/workplace/software/FastQC/fastqc -f fastq --noextract -contaminants /workplace/software/FastQC/Contaminants/contaminant_list.txt  %s'
	os.system(fastqc_cmd % fastqfile)

def move_fastqc(files):
	# fastqc
	try:
		shutil.rmtree('fastqc')
		os.mkdir('fastqc')
	except:
		os.mkdir('fastqc')
	for file in files:
		# KD-A_sequence_1_fastqc.html
		shutil.move(file[:-6]+'_1_fastqc.html', 'fastqc')
		shutil.move(file[:-6]+'_1_fastqc.zip', 'fastqc')
		shutil.move(file[:-6]+'_2_fastqc.html', 'fastqc')
		shutil.move(file[:-6]+'_2_fastqc.zip', 'fastqc')

def pfastqc(fqfiles):
	"""
	run fastqc parallel x_sequence_1.fastq
	"""
	threads = []
	newfiles = []
	for file in fqfiles:
		newfile1 = file[:-6]+'_1.fastq'
		newfile2 = file[:-6]+'_2.fastq'
		newfiles.append(newfile1)
		newfiles.append(newfile2)
	for i in newfiles:
		t = threading.Thread(target=_fastqc, args=(i,))
		threads.append(t)
	for t in threads:
		t.start()
	for t in threads:
		t.join()
	move_fastqc(fqfiles)
	

# add Q30 calculation, just by one sample 20160822
def q30(fqfile):
	"""
	for just a fastq file
	"""
	f = open(fqfile)
	lines=f.readlines()
	Q =0.0
	Q30=0.0
	for i in xrange(3,len(lines),4):
		#print i
		line=lines[i].strip()
		Q+=len(line)
		Q30+=len([jj for jj in line if (ord(jj)-33)>=30])
	f.close()
	return Q30/Q*100

def getq30(fqfiles):
	f=open('q30.txt','w')
	for fqfile in fqfiles:
		fq1 = fqfile[:-6]+'_1.fastq'
		fq2 = fqfile[:-6]+'_2.fastq'
		q301 = q30(fq1)
		q302 = q30(fq2)
		q = (q301+q302)/2.0
		print >>f, fqfile.replace('_sequence.fastq','')+'\t%2.2f%%' % q
	f.close()


############################# 3, over


############################# 4, cut 3' adaptor, this is for paired-end results, for single-end, please rewrite


def _cutadapt(fastqfile):
	"""
	cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq
	"""
	#cutadaptcmd = './clipPairedEnd.pl -m1 tS0_sequence_1.fastq -m2 tS0_sequence_2.fastq -o1 tS0_sequence.cutadapt_1.fastq -o2 tS0_sequence.cutadapt_2.fastq -a1 GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGT -s1 S01.stat -s2 S02.stat'
	cutadaptcmd = 'cutadapt -m 20 -q 15 -a %s -A %s -o %s -p %s %s %s' % (three_prime_adaptor_1, three_prime_adaptor_2,  fastqfile[:-6]+'.cutadapt_1.fastq',fastqfile[:-6]+'.cutadapt_2.fastq', fastqfile[:-6]+'_1.fastq', fastqfile[:-6]+'_2.fastq')
	os.system(cutadaptcmd)


def pcutadapt(fastqfiles):
	threads = []
	for i in fastqfiles:
		t = threading.Thread(target=_cutadapt, args=(i,))
		threads.append(t)
	for t in threads:
		t.start()
	for t in threads:
		t.join()

############################# 4, cut 3' adaptor over

############################# 5, QC: rRNA percentage, do not run everytime

def fastq_lines_by_4(file):
	"""
	report the read counts
	"""
	f = open(file)
	n = 0
	for line in f:
		n+=1
	f.close()
	return n/4.0

############################# 5, QC: rRNA percentage, over

############################ 6, get fasta file

def leaders(xs):
	counts = defaultdict(int)
	for x in xs:
		counts[x] += 1
	return sorted(counts.items(), reverse=True, key=lambda tup: tup[1])
   
def counts(fqfile):
	"""
	count the seq and counts
	"""
	f = open(fqfile)
	n = 0
	s = []
	for line in f:
		if n%4 == 1:
			s.append(line.strip())
		n+=1
	f.close()
	f2=open(fqfile+'.fa','w')
	x=0
	for i in leaders(s):
		print >>f2, '>%s_%d_x%d' % (fqfile[:-6], x,i[1])
		print >>f2, i[0]
		x+=1
	f2.close()

def fastq2fa(FASTQFILES):
	for file in FASTQFILES:
		counts(file[0:-6]+'.cutadapt_1.fastq')
		counts(file[0:-6]+'.cutadapt_2.fastq')
		
	
################################ 6, get fasta file over



###### change all the figure into size of 8cm (by liwei, 20170405)

def cm2inch(*tupl):  # add 20170405
	"""
	figure(figsize=cm2inch(8,8))，改下figsize的参数即可
	参考：http://stackoverflow.com/questions/14708695/specify-figure-size-in-centimeter-in-matplotlib
	"""	
	inch= 2.54
	if isinstance(tupl[0], tuple):
		return tuple(i/inch for i in tupl[0])
	else:
		return tuple(i/inch for i in tupl)


################################# 7, plot of length distribution
def total_tag_num(file, countthreshold):
	"""
	sum of tag count when tag count >= countthreshold
	"""
	f = open(file)
	allnum = []
	for line in f:
		if '>' in line:
			allnum.append(float(line.strip().split('_')[-1][1:])) # >R3+R4_1_x341392
	f.close()
	return sum([s for s in allnum if s>=countthreshold])

def tag_length_compact(file, countthreshold):
	"""
	[[15,5555],[16,6666],...] count >= countthreshold
	"""
	f = open(file)
	all_id_seq = []
	allchar = f.read()
	f.close()
	id_seqs = allchar.split('>')[1:]
	length_count_dict_tmp = {}
	for id_seq in id_seqs:
		idnum, seq = id_seq.split('\n')[0:2]
		n = int(idnum.split('_')[-1][1:])
		length_count_dict_tmp.setdefault(len(seq),[]).append(n)
	#print length_count_dict_tmp # ok
	length_counts = {}
	for x in length_count_dict_tmp.keys():
		length_counts[x] = sum([i for i in length_count_dict_tmp[x] if i>=countthreshold])
	return length_counts

def tag_length_compact_unique(file, countthreshold):
	"""
	same to tag_length_compact, but tag count was unique (a tag count once)
	"""
	f = open(file)
	all_id_seq = []
	allchar = f.read()
	f.close()
	id_seqs = allchar.split('>')[1:]
	length_count_dict_tmp = {}
	for id_seq in id_seqs:
		idnum, seq = id_seq.split('\n')[0:2]
		n = 1
		length_count_dict_tmp.setdefault(len(seq),[]).append(n)
	#print length_count_dict_tmp # ok
	length_counts = {}
	for x in length_count_dict_tmp.keys():
		length_counts[x] = sum([i for i in length_count_dict_tmp[x] if i>=countthreshold])
	return length_counts
    

#print tag_length_compact('G2_trimmed_tags.fa', 1)

def insert0(datadict, lists):
	new = {}
	for i in lists:
		try:
			new[i] = datadict[i]
		except:
			new[i] = 0
	return new


def autolabel(rects, name):
	# attach some text labels
	for ii,rect in enumerate(rects):
		height = rect.get_height()
		plt.text(rect.get_x()+rect.get_width()/2., 1.02*height, '%s'% (name[ii]),ha='center', va='bottom',fontsize=6,color='black')
 

def plot_distribution(file, countthreshold, minlength, maxlength):
	data = tag_length_compact(file, countthreshold)  # {15:100000,16:120000,...}
	ind = range(minlength, maxlength+1)  #[15,...,30]
	datatmp = insert0(data, ind)
	bardata = [datatmp[i] for i in ind]     #
	#print data, ind, bardata
	fig = figure(figsize=(len(ind)/1.5,8))
	rects = bar(ind, bardata, color='red')
	xticks([x+0.4 for x in ind], [str(i) for i in ind],fontsize="x-small")
	xlim(ind[0]-1, ind[-1]+1)
	ylim(0, max(bardata)*1.1)
	names= [str(i) for i in bardata]
	autolabel(rects,names)
	ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	title('Length Distribution: %s' % file.replace('_sequence.fastq.cutadapt.fa',''))
	savefig(file+'.png', dpi=600.1)
	png2tiff(file+'.png')


def plot_distribution_unique(file, countthreshold, minlength, maxlength):
	data = tag_length_compact_unique(file, countthreshold)  # {15:100000,16:120000,...}
	ind = range(minlength, maxlength+1)  #[15,...,30]
	datatmp = insert0(data, ind)
	bardata = [datatmp[i] for i in ind]     #
	#print data, ind, bardata
	fig = figure(figsize=(len(ind)/1.5,8))
	rects = bar(ind, bardata, color='red')
	xticks([x+0.4 for x in ind], [str(i) for i in ind],fontsize="x-small")
	xlim(ind[0]-1, ind[-1]+1)
	ylim(0, max(bardata)*1.1)
	names= [str(i) for i in bardata]
	autolabel(rects,names)
	ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	title('Length Distribution: %s' % (file.replace('_sequence.cutadapt_','').replace('.fastq.fa','')))
	savefig(file+'.unique.png', dpi=301)
	png2tiff(file+'.unique.png')
 

def plot_length_distribution(FASTQFILES, countthreshold=1, minlength=15, maxlength=30):
	for file in FASTQFILES:
		plot_distribution(file[0:-6]+'.cutadapt_1.fastq.fa', countthreshold, minlength, maxlength)
		plot_distribution(file[0:-6]+'.cutadapt_2.fastq.fa', countthreshold, minlength, maxlength)
		plot_distribution_unique(file[0:-6]+'.cutadapt_1.fastq.fa', countthreshold, minlength, maxlength)
		plot_distribution_unique(file[0:-6]+'.cutadapt_2.fastq.fa', countthreshold, minlength, maxlength)
	
###################################### 7 plot length distribution over

###################################### 8 hisat2

"""
nohup /root/hisat2-master/hisat2-build  ucsc_galGal4.fa ucsc_galGal4 &


HISAT2 version 2.0.4 by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo)
Usage: hisat2-build [options]* <reference_in> <bt2_index_base>
    reference_in            comma-separated list of files with ref sequences
    hisat2_index_base          write ht2 data to files with this dir/basename
Options:
    -c                      reference sequences given on cmd line (as
                            <reference_in>)
    --large-index           force generated index to be 'large', even if ref
                            has fewer than 4 billion nucleotides
    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting
    -p                      number of threads
    --bmax <int>            max bucket sz for blockwise suffix-array builder
    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)
    --dcv <int>             diff-cover period for blockwise (default: 1024)
    --nodc                  disable diff-cover (algorithm becomes quadratic)
    -r/--noref              don't build .3/.4.bt2 (packed reference) portion
    -3/--justref            just build .3/.4.bt2 (packed reference) portion
    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)
    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)
    --localoffrate <int>    SA (local) is sampled every 2^offRate BWT chars (default: 3)
    --localftabchars <int>  # of chars consumed in initial lookup in a local index (default: 6)
    --snp <path>            SNP file name
    --haplotype <path>      haplotype file name
    --ss <path>             Splice site file name
    --exon <path>           Exon file name
    --seed <int>            seed for random number generator
    -q/--quiet              verbose output (for debugging)
    -h/--help               print detailed description of tool and its options
    --usage                 print this usage message
    --version               print version information and quit

"""
 

def _run_hisat2_rRNA(file):  # the successor of tophat2
	"""
HISAT2 version 2.0.4 by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)
Usage: 
  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

  <ht2-idx>  Index filename prefix (minus trailing .X.ht2).
             索引文件名，首先要自己建立好，略
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             双端中1端的文件，与m2中的文件配对，可以是.gz或者.bz2格式的，这里选择.gz格式作为输入（省空间），先用ln -s xxx.gz yyy.gz建立软连接
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             双端中2端的文件，与m1中的文件配对，可以是.gz或者.bz2格式的，这里选择.gz格式作为输入（省空间），先用ln -s xxx.gz yyy.gz建立软连接
  <r>        Files with unpaired reads.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             非配对的read文件。可以是.gz或者.bz2格式的，这里选择.gz格式作为输入（省空间）
  <sam>      File for SAM output (default: stdout)
             SAM输出的文件（默认是stdout）

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.  一般不用这个

Options (defaults in parentheses):

 Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
                     输入文件是FASTQ格式（默认）
  --qseq             query input files are in Illumina's qseq format
                     输入文件是qseq格式
  -f                 query input files are (multi-)FASTA .fa/.mfa
                     输入文件是.fa格式
  -r                 query input files are raw one-sequence-per-line
                     输入文件是原始的一条序列一行格式
  -c                 <m1>, <m2>, <r> are sequences themselves, not files
  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
                     跳过输入文件的前int reads （默认：不跳过，即从头读）
  -u/--upto <int>    stop after first <int> reads/pairs (no limit)
                     在输入文件的第int reads处停止 （默认：无限制，即读到文件尾）
  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
                     从5'端（左末端）去掉int个碱基（默认：不去）
  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
                     从3'端（右末端）去掉int个碱基（默认：不去）
  --phred33          qualities are Phred+33 (default)
                     质量值是Phred+33（默认）
  --phred64          qualities are Phred+64
                     质量值是Phred+66（默认）
  --int-quals        qualities encoded as space-delimited integers

 Alignment:
  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
                     seed比对重最大错配数，可以是0或1（默认0）
  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
  --ignore-quals     treat all quality values as 30 on Phred scale (off)
  --nofw             do not align forward (original) version of read (off)
  --norc             do not align reverse-complement version of read (off)

 Spliced Alignment:
  --pen-cansplice <int>              penalty for a canonical splice site (0)
  --pen-noncansplice <int>           penalty for a non-canonical splice site (12)
  --pen-canintronlen <func>          penalty for long introns (G,-8,1) with canonical splice sites
  --pen-noncanintronlen <func>       penalty for long introns (G,-8,1) with noncanonical splice sites
  --min-intronlen <int>              minimum intron length (20)
  --max-intronlen <int>              maximum intron length (500000)
  --known-splicesite-infile <path>   provide a list of known splice sites
  --novel-splicesite-outfile <path>  report a list of splice sites
  --novel-splicesite-infile <path>   provide a list of novel splice sites
  --no-temp-splicesite               disable the use of splice sites found
  --no-spliced-alignment             disable spliced alignment
  --rna-strandness <string>          Specify strand-specific information (unstranded)
  --tmo                              Reports only those alignments within known transcriptome
  --dta                              Reports alignments tailored for transcript assemblers
  --dta-cufflinks                    Reports alignments tailored specifically for cufflinks
                                     输出针对cufflinks的比对结果

 Scoring:
  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
                     匹配奖励（端-端是0，局部是2）
  --mp <int>,<int>   max and min penalties for mismatch; lower qual = lower penalty <2,6>
                     错配的罚分
  --sp <int>,<int>   max and min penalties for soft-clipping; lower qual = lower penalty <1,2>
  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
  --rdg <int>,<int>  read gap open, extend penalties (5,3)
  --rfg <int>,<int>  reference gap open, extend penalties (5,3)
  --score-min <func> min acceptable alignment score w/r/t read length
                     (L,0.0,-0.2)

 Reporting:
  (default)          look for multiple alignments, report best, with MAPQ
                     寻找多重匹配，以MAPQ为标准，报告最好的那个（默认）
   OR
  -k <int>           report up to <int> alns per read; MAPQ not meaningful
                     报告每个read高达int的比对，MAPQ无意义
   OR
  -a/--all           report all alignments; very slow, MAPQ not meaningful
                     报告全部比对，非常慢，MAPQ无意义

 Effort:
  -D <int>           give up extending after <int> failed extends in a row (15)
  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)

 Paired-end:
  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
                     -1，-2比对到forward/reverse，默认--fr（我们用的也是）
  --no-mixed         suppress unpaired alignments for paired reads
  --no-discordant    suppress discordant alignments for paired reads

 Output:
  -t/--time          print wall-clock time taken by search phases
                     打印出搜索时间
  --un <path>           write unpaired reads that didn't align to <path>
  --al <path>           write unpaired reads that aligned at least once to <path>
  --un-conc <path>      write pairs that didn't align concordantly to <path>
  --al-conc <path>      write pairs that aligned concordantly at least once to <path>
  (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
  --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
  --quiet            print nothing to stderr except serious errors
                     除了严重错误外，不打印任何标准误
  --met-file <path>  send metrics to file at <path> (off)
  --met-stderr       send metrics to stderr (off)
  --met <int>        report internal counters & metrics every <int> secs (1)
  --no-head          supppress header lines, i.e. lines starting with @
  --no-sq            supppress @SQ header lines
  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
  --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                     Note: @RG line only printed when --rg-id is set.
  --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments.

 Performance:
  -o/--offrate <int> override offrate of index; must be >= index's offrate
  -p/--threads <int> number of alignment threads to launch (1)
                     启用的比对线程，默认1，指定28
  --reorder          force SAM output order to match order of input reads
  --mm               use memory-mapped I/O for index; many 'bowtie's can share

 Other:
  --qc-filter        filter out reads that are bad according to QSEQ filter
  --seed <int>       seed for random number generator (0)
                     随机数产生器的seed，默认0
  --non-deterministic seed rand. gen. arbitrarily instead of using read attributes
  --version          print version information and quit
                     打印版本，并退出
  -h/--help          print this usage message

    =====================================================================================
	"""
	#/root/hisat2-master/hisat2 -p 28 --dta-cufflinks /workplace/database/human/hisat2_index/ucsc_hg19 -1 T1_sequence.cutadapt_1.fastq -2 T1_sequence.cutadapt_2.fastq -S T1.sam')
	# mapping
	hisat2_rRNA_index = cf.get('QC','rRNA_index') # /workplace/database/human/hisat2_index/ucsc_hg19
	cmd = "/root/hisat2-master/hisat2 -p 12 --dta-cufflinks  -x %s  -1 %s -2 %s -S %s.rRNA.sam>%s.rRNA.mapping.rate.txt 2>&1" % (hisat2_rRNA_index, file[0:-6]+'.cutadapt_1.fastq', file[0:-6]+'.cutadapt_2.fastq', file[0:-6], file[0:-6])
	print '----------------------------hisat2 rRNA parameters----------------------------------'
	print cmd
	print '--------------------------------------------------------------------------------'
	os.system(cmd)
	os.remove('%s.rRNA.sam' % file[0:-6])


def run_hisat2_rRNA(files):
	"""
	for all files, paired end
	since use 28 core, do not parallel anymore
	"""
	print '-----------------------------------rRNA mapping results start-------------------------'
	for file in files:
		_run_hisat2_rRNA(file)
	print '-----------------------------------rRNA mapping results over--------------------------'



def mapping_rate(file):
	f=open(file)
	lines=f.readlines()
	f.close()
	#print lines[3],lines[4],lines[7],lines[12],lines[13]
	al=float(lines[0].strip().split()[0])*2
	a=int(lines[3].strip().split()[0])*2
	b=int(lines[4].strip().split()[0])*2
	c=int(lines[7].strip().split()[0])*2
	d=int(lines[12].strip().split()[0])
	e=int(lines[13].strip().split()[0])
	#print float(a+b+c+d+e)
	#print lines[14].strip().split()[0]
	return '%s\t%s\t%s' % (file, str(float(a+b+c+d+e)), lines[14].strip().split()[0])
	
def rRNA_mapping_rates(files):
	f=open('rRNA.mapping.rate.txt','w')
	for file in files:
		print >>f, mapping_rate('%s.rRNA.mapping.rate.txt' % file[0:-6])
		os.remove('%s.rRNA.mapping.rate.txt' % file[0:-6])
	f.close()


def _run_hisat2(file):  # the successor of tophat2
	"""
HISAT2 version 2.0.4 by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)
Usage: 
  hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]

  <ht2-idx>  Index filename prefix (minus trailing .X.ht2).
             索引文件名，首先要自己建立好，略
  <m1>       Files with #1 mates, paired with files in <m2>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             双端中1端的文件，与m2中的文件配对，可以是.gz或者.bz2格式的，这里选择.gz格式作为输入（省空间），先用ln -s xxx.gz yyy.gz建立软连接
  <m2>       Files with #2 mates, paired with files in <m1>.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             双端中2端的文件，与m1中的文件配对，可以是.gz或者.bz2格式的，这里选择.gz格式作为输入（省空间），先用ln -s xxx.gz yyy.gz建立软连接
  <r>        Files with unpaired reads.
             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).
             非配对的read文件。可以是.gz或者.bz2格式的，这里选择.gz格式作为输入（省空间）
  <sam>      File for SAM output (default: stdout)
             SAM输出的文件（默认是stdout）

  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be
  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'.  一般不用这个

Options (defaults in parentheses):

 Input:
  -q                 query input files are FASTQ .fq/.fastq (default)
                     输入文件是FASTQ格式（默认）
  --qseq             query input files are in Illumina's qseq format
                     输入文件是qseq格式
  -f                 query input files are (multi-)FASTA .fa/.mfa
                     输入文件是.fa格式
  -r                 query input files are raw one-sequence-per-line
                     输入文件是原始的一条序列一行格式
  -c                 <m1>, <m2>, <r> are sequences themselves, not files
  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)
                     跳过输入文件的前int reads （默认：不跳过，即从头读）
  -u/--upto <int>    stop after first <int> reads/pairs (no limit)
                     在输入文件的第int reads处停止 （默认：无限制，即读到文件尾）
  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)
                     从5'端（左末端）去掉int个碱基（默认：不去）
  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)
                     从3'端（右末端）去掉int个碱基（默认：不去）
  --phred33          qualities are Phred+33 (default)
                     质量值是Phred+33（默认）
  --phred64          qualities are Phred+64
                     质量值是Phred+66（默认）
  --int-quals        qualities encoded as space-delimited integers

 Alignment:
  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
                     seed比对重最大错配数，可以是0或1（默认0）
  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)
  --ignore-quals     treat all quality values as 30 on Phred scale (off)
  --nofw             do not align forward (original) version of read (off)
  --norc             do not align reverse-complement version of read (off)

 Spliced Alignment:
  --pen-cansplice <int>              penalty for a canonical splice site (0)
  --pen-noncansplice <int>           penalty for a non-canonical splice site (12)
  --pen-canintronlen <func>          penalty for long introns (G,-8,1) with canonical splice sites
  --pen-noncanintronlen <func>       penalty for long introns (G,-8,1) with noncanonical splice sites
  --min-intronlen <int>              minimum intron length (20)
  --max-intronlen <int>              maximum intron length (500000)
  --known-splicesite-infile <path>   provide a list of known splice sites
  --novel-splicesite-outfile <path>  report a list of splice sites
  --novel-splicesite-infile <path>   provide a list of novel splice sites
  --no-temp-splicesite               disable the use of splice sites found
  --no-spliced-alignment             disable spliced alignment
  --rna-strandness <string>          Specify strand-specific information (unstranded)
  --tmo                              Reports only those alignments within known transcriptome
  --dta                              Reports alignments tailored for transcript assemblers
  --dta-cufflinks                    Reports alignments tailored specifically for cufflinks
                                     输出针对cufflinks的比对结果

 Scoring:
  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) 
                     匹配奖励（端-端是0，局部是2）
  --mp <int>,<int>   max and min penalties for mismatch; lower qual = lower penalty <2,6>
                     错配的罚分
  --sp <int>,<int>   max and min penalties for soft-clipping; lower qual = lower penalty <1,2>
  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)
  --rdg <int>,<int>  read gap open, extend penalties (5,3)
  --rfg <int>,<int>  reference gap open, extend penalties (5,3)
  --score-min <func> min acceptable alignment score w/r/t read length
                     (L,0.0,-0.2)

 Reporting:
  (default)          look for multiple alignments, report best, with MAPQ
                     寻找多重匹配，以MAPQ为标准，报告最好的那个（默认）
   OR
  -k <int>           report up to <int> alns per read; MAPQ not meaningful
                     报告每个read高达int的比对，MAPQ无意义
   OR
  -a/--all           report all alignments; very slow, MAPQ not meaningful
                     报告全部比对，非常慢，MAPQ无意义

 Effort:
  -D <int>           give up extending after <int> failed extends in a row (15)
  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)

 Paired-end:
  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)
                     -1，-2比对到forward/reverse，默认--fr（我们用的也是）
  --no-mixed         suppress unpaired alignments for paired reads
  --no-discordant    suppress discordant alignments for paired reads

 Output:
  -t/--time          print wall-clock time taken by search phases
                     打印出搜索时间
  --un <path>           write unpaired reads that didn't align to <path>
  --al <path>           write unpaired reads that aligned at least once to <path>
  --un-conc <path>      write pairs that didn't align concordantly to <path>
  --al-conc <path>      write pairs that aligned concordantly at least once to <path>
  (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g.
  --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)
  --quiet            print nothing to stderr except serious errors
                     除了严重错误外，不打印任何标准误
  --met-file <path>  send metrics to file at <path> (off)
  --met-stderr       send metrics to stderr (off)
  --met <int>        report internal counters & metrics every <int> secs (1)
  --no-head          supppress header lines, i.e. lines starting with @
  --no-sq            supppress @SQ header lines
  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
  --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                     Note: @RG line only printed when --rg-id is set.
  --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments.

 Performance:
  -o/--offrate <int> override offrate of index; must be >= index's offrate
  -p/--threads <int> number of alignment threads to launch (1)
                     启用的比对线程，默认1，指定28
  --reorder          force SAM output order to match order of input reads
  --mm               use memory-mapped I/O for index; many 'bowtie's can share

 Other:
  --qc-filter        filter out reads that are bad according to QSEQ filter
  --seed <int>       seed for random number generator (0)
                     随机数产生器的seed，默认0
  --non-deterministic seed rand. gen. arbitrarily instead of using read attributes
  --version          print version information and quit
                     打印版本，并退出
  -h/--help          print this usage message

    =====================================================================================
	"""
	#/root/hisat2-master/hisat2 -p 28 --dta-cufflinks /workplace/database/human/hisat2_index/ucsc_hg19 -1 T1_sequence.cutadapt_1.fastq -2 T1_sequence.cutadapt_2.fastq -S T1.sam')
	# mapping >%s.rRNA.mapping.rate.txt 2>&1
	hisat2_index = cf.get('mapping_parameters','hisat2_index') # /workplace/database/human/hisat2_index/ucsc_hg19
	#hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')     # Homo_sapiens_HG19_for_RNAseq.gtf
	cmd = "/root/hisat2-master/hisat2 -p 12 --dta-cufflinks  -x %s  -1 %s -2 %s -S %s.sam>%s.mapping.rate.txt 2>&1" % (hisat2_index, file[0:-6]+'.cutadapt_1.fastq', file[0:-6]+'.cutadapt_2.fastq', file[0:-6], file[0:-6])
	print '----------------------------hisat2 parameters----------------------------------'
	print cmd
	print '--------------------------------------------------------------------------------'
	os.system(cmd)


def run_hisat2(files):
	"""
	for all files, paired end
	since use 28 core, do not parallel anymore
	"""
	for file in files:
		_run_hisat2(file)

def overall_mapping_rates(files):
	d ={}
	for fqfile in files:
		d[fqfile] = mapping_rate('%s.mapping.rate.txt' % fqfile[0:-6])
		os.remove('%s.mapping.rate.txt' % fqfile[0:-6])
	return d

###################################### 8 hisat2, over


###################################### 9 run sam to sorted bam

def sam22bam(file):
	os.system('/workplace/software/samtools-1.2/samtools view -bS %s.sam > %s.bam'  % (file[:-6], file[:-6]))
	os.remove('%s.sam' % file[:-6])
	os.system('/workplace/software/samtools-1.2/samtools sort %s.bam %s.sorted' % (file[:-6], file[:-6]))
	os.system('/workplace/software/samtools-1.2/samtools index %s.sorted.bam' % (file[:-6]))
	os.remove('%s.bam' % file[:-6])

def psam2bam(fastqfiles):
	threads = []
	for i in fastqfiles:
		t = threading.Thread(target=sam22bam, args=(i,))
		threads.append(t)
	for t in threads:
		t.start()
	for t in threads:
		t.join()


###################################### 9 run sam to sorted bam over

###################################### 10 run cuffquant

def _run_cuffquant(file):
	"""
	cuffquant v2.2.1 (4237)
-----------------------------
Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]
   Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam
General Options:
  -o/--output-dir              write all output files to this directory              [ default:     ./ ]
  -M/--mask-file               ignore all alignment within transcripts in this file  [ default:   NULL ]
  -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]
  -u/--multi-read-correct      use 'rescue method' for multi-reads                   [ default:  FALSE ]
  -p/--num-threads             number of threads used during quantification          [ default:      1 ]
  --library-type               Library prep used for input reads                     [ default:  below ]

Advanced Options:
  -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]
  -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]
  -c/--min-alignment-count     minimum number of alignments in a locus for testing   [ default:   10 ]
  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]
  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
  --seed                       value of random number generator seed                 [ default:      0 ]
  --no-update-check            do not contact server to check for update availability[ default:  FALSE ]
  --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]
  --max-frag-multihits         Maximum number of alignments allowed per fragment     [ default: unlim  ]
  --no-effective-length-correction   No effective length correction                  [ default:  FALSE ]
  --no-length-correction       No length correction                                  [ default:  FALSE ]

Debugging use only:
  --read-skip-fraction         Skip a random subset of reads this size               [ default:    0.0 ]
  --no-read-pairs              Break all read pairs                                  [ default:  FALSE ]
  --trim-read-length           Trim reads to be this long (keep 5' end)              [ default:   none ]
  --no-scv-correction          Disable SCV correction                                [ default:  FALSE ]

Supported library types:
        ff-firststrand
        ff-secondstrand
        ff-unstranded
        fr-firststrand
        fr-secondstrand
        fr-unstranded (default)
        transfrags
	
	
	here use 28 core, do not parallel anymore
	"""
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')     # Homo_sapiens_HG19.gtf
	os.system('/root/cufflinks-2.2.1.Linux_x86_64/cuffquant --library-type fr-firststrand -p 20 %s %s.sorted.bam -o %s' % (hisat2_gtf, file[:-6], file.split('_sequence.fastq')[0])) 


def obsolete_run_cuffquant(fastqfiles):  # change into parallel
	for fastqfile in fastqfiles:
		_run_cuffquant(fastqfile)


def prun_cuffquant(fastqfiles):
	threads = []
	for i in fastqfiles:
		t = threading.Thread(target=_run_cuffquant, args=(i,))
		threads.append(t)
	for t in threads:
		t.start()
	for t in threads:
		t.join()


###################################### 10 run cuffquant over


###################################### 11 run cuffdiff get expression, and for two samples comparison

def run_cuffdiff_1vs1():
	"""
	cuffdiff v2.2.1 (4237)
-----------------------------
Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]
   Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam
General Options:
  -o/--output-dir              write all output files to this directory              [ default:     ./ ]
                               写所有输出文件到这个目录                                 默认当目录
  -L/--labels                  comma-separated list of condition labels
                               逗号分割的条件标记列表
  --FDR                        False discovery rate used in testing                  [ default:   0.05 ]
                               测试中使用的错误发现率                                   默认: 0.05
  -M/--mask-file               ignore all alignment within transcripts in this file  [ default:   NULL ]
  -C/--contrast-file           Perform the constrasts specified in this file         [ default:   NULL ]
  -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]
  -u/--multi-read-correct      use 'rescue method' for multi-reads                   [ default:  FALSE ]
  -p/--num-threads             number of threads used during quantification          [ default:      1 ]
                               定量中使用的线程数                                       默认1，改成28
  --no-diff                    Don't generate differential analysis files            [ default:  FALSE ]
                               不产生差异分析文件                                       默认：产生
  --no-js-tests                Don't perform isoform switching tests                 [ default:  FALSE ]
                               不执行isoform交换检测                                   默认：执行
  -T/--time-series             treat samples as a time-series                        [ default:  FALSE ]
                               将样品看做时序的                                        默认：不看成时序的
  --library-type               Library prep used for input reads                     [ default:  below ]
                               输入reads的文库制备方法，见下
  --dispersion-method          Method used to estimate dispersion models             [ default:  below ]
                               评估离散模型的方法，见下
  --library-norm-method        Method used to normalize library sizes                [ default:  below ]
                               标准化文库大小的方法，见下

Advanced Options:
  -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]
                               平均片段长度（仅针对非配对reads）                         默认：200
  -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]
                               片段长度标准差（仅针对非配对reads）                       默认：80
  -c/--min-alignment-count     minimum number of alignments in a locus for testing   [ default:   10 ]
                               用于检验的一个基因座的最小比对数                          默认：10
  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]
                               mle计算允许的最大迭代次数                                默认：5000
  --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:   TRUE ]
                               仅计算与参考RNA兼容的计数击中                             默认：TRUE
  --total-hits-norm            count all hits for normalization                      [ default:  FALSE ]
                               计算所有的hits用于标准化                                 默认：FALSE
  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
                               友好日志输出（内容丰富）                                 默认：FALSE
  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
                               安静地日志输出                                          默认：FALSE
  --seed                       value of random number generator seed                 [ default:      0 ]
                               随机数产生器的seed值                                    默认：0
  --no-update-check            do not contact server to check for update availability[ default:  FALSE ]
                               不连接服务器检查更新                                     默认连接，可去掉
  --emit-count-tables          print count tables used to fit overdispersion         [    DEPRECATED   ]
                                                                                     弃用
  --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]                               
  --num-frag-count-draws       Number of fragment generation samples                 [ default:    100 ]
  --num-frag-assign-draws      Number of fragment assignment samples per generation  [ default:     50 ]
  --max-frag-multihits         Maximum number of alignments allowed per fragment     [ default: unlim  ]
  --min-outlier-p              Min replicate p value to admit for testing            [    DEPRECATED   ]
                                                                                     弃用
  --min-reps-for-js-test       Replicates needed for relative isoform shift testing  [ default:      3 ]
  --no-effective-length-correction   No effective length correction                  [ default:  FALSE ]
  --no-length-correction       No length correction                                  [ default:  FALSE ]
  -N/--upper-quartile-norm     Deprecated, use --library-norm-method                 [    DEPRECATED   ]
                                                                                     弃用
  --geometric-norm             Deprecated, use --library-norm-method                 [    DEPRECATED   ]
                                                                                     弃用
  --raw-mapped-norm            Deprecated, use --library-norm-method                 [    DEPRECATED   ]
                                                                                     弃用
  --poisson-dispersion         Deprecated, use --dispersion-method                   [    DEPRECATED   ]
                                                                                     弃用

Debugging use only:
  --read-skip-fraction         Skip a random subset of reads this size               [ default:    0.0 ]
  --no-read-pairs              Break all read pairs                                  [ default:  FALSE ]
  --trim-read-length           Trim reads to be this long (keep 5' end)              [ default:   none ]
  --no-scv-correction          Disable SCV correction                                [ default:  FALSE ]

Supported library types:
        ff-firststrand
        ff-secondstrand
        ff-unstranded
        fr-firststrand         我们用的这个？待lw确认
        fr-secondstrand
        fr-unstranded (default)
        transfrags



ff-firststrand
ff-secondstrand
ff-unstranded
fr-firststrand	dUTP, NSR, NNSR	与上类似，除了片段最右端的首先测序（或者仅测序单端）。相当于仅合成的第一条链被测序
fr-secondstrand	Directional Illumina (Ligation), Standard SOLiD	与上雷同，除了片段最左端的首先测序（或仅测序单端）。相当于仅合成的第一条链被测序
fr-unstranded (default)	Standard Illumina	来自片段最左端的（转录坐标）比对到转录本链，而最右端的比对到相反链
transfrags  


Supported dispersion methods:
        blind
        per-condition
        poisson
        pooled (default)

blind	All samples are treated as replicates of a single global “condition” and used to build one model.
        所有样品认为是一个“全局条件”的重复，并被用来构建一个模型
per-condition	Each replicated condition receives its own model. Only available when all conditions have replicates.
                每个重复条件使用它自身的模型。仅所有条件都有重复时才有用
poisson	
pooled	Each replicated condition is used to build a model, then these models are averaged to provide a single global model for all conditions in the experiment. (Default)
        使用每个重复条件构建模型，然后这些模型平均以提供单个全局模型用于试验中的所有条件（默认）

Supported library normalization methods:
支持的标准化方法：
        classic-fpkm
        库大小因子设置成1。即对FPKM值或片段数不设置scaling（Cufflinks默认的）。我们使用的也是这个。
        geometric (default)
        FPKMs and fragment counts are scaled via the median of the geometric means of fragment counts across all libraries, as described in Anders and Huber (Genome Biology, 2010). This policy identical to the one used by DESeq. (default for Cuffdiff)
        quartile
        FPKMs and fragment counts are scaled via the ratio of the 75 quartile fragment counts to the average 75 quartile value across all libraries.       
	"""
	samples  = getfastq()
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')     # Homo_sapiens_HG19.gtf
	hisat2_fasta = cf.get('mapping_parameters','hisat2_fasta')
	# expression results, no differentially expression
	try:
		shutil.rmtree('cuffdiff')
		os.mkdir('cuffdiff')
	except:
		os.mkdir('cuffdiff')
	os.chdir('cuffdiff')  # for expression
	inputname = ','.join([i.split('_sequence.fastq')[0] for i in samples])
	inputcxb=  ' '.join(['../'+i.split('_sequence.fastq')[0]+'/abundances.cxb' for i in samples])
	cmd = '/root/cufflinks-2.2.1.Linux_x86_64/cuffdiff --library-norm-method classic-fpkm --library-type  fr-firststrand --no-diff --no-update-check -u -v --emit-count-tables -p 20 -L %s -b %s %s %s' % (inputname, hisat2_fasta, hisat2_gtf, inputcxb)
	print cmd
	os.system(cmd)
	os.chdir(global_now_location)
	
	# 2 samples comparisons
	comparisons = get_2_samples_comparisons()
	if comparisons:
		for comparison in comparisons:
			print comparison
			sampleA, sampleB, FC, P = comparison
			foldname = 'cuffdiff_'+ sampleA.split('_sequence.fastq')[0]+'_vs_'+sampleB.split('_sequence.fastq')[0]
			try:
				shutil.rmtree(foldname)
				os.mkdir(foldname)
			except:
				os.mkdir(foldname)
			os.chdir(foldname)
			cmd2 = '/root/cufflinks-2.2.1.Linux_x86_64/cuffdiff --library-norm-method classic-fpkm --library-type fr-firststrand --no-update-check -u -v --emit-count-tables -p 20 -L %s -b %s %s %s' % ((sampleA+','+sampleB), hisat2_fasta, hisat2_gtf, ('../'+sampleA.split('_sequence.fastq')[0]+'/abundances.cxb'+' '+'../'+sampleB.split('_sequence.fastq')[0]+'/abundances.cxb'))
			print cmd2
			os.system(cmd2)
			os.chdir(global_now_location)
	else:
		print '!!!!!!!!!! There is no 2 samples comparison'


###################################### 11 run cuffdiff get expression, and for two samples comparison, over


###################################### 15 run cuffdiff for two groups
def remove_space(string):
	#p=re.compile('\s+')
	#return re.sub(p,'',string)
	plist=re.split(",",string)
	templist=[]
	for i in plist:
		templist.append(i.strip())
	p=",".join(templist)
	return p

def get_2group_comparisons(): # change into P
	"""[two_group_comparisons]
	two_group_comparison_1_groupA = G1_sequence.fastq, G2_sequence.fastq, G3_sequence.fastq
	two_group_comparison_1_groupA_name = G123
	two_group_comparison_1_groupB = N1_sequence.fastq, N2_sequence.fastq, N3_sequence.fastq
	two_group_comparison_1_groupB_name = N123
	two_group_comparison_1_FC = 2.0
	two_group_comparison_1_P  = 0.05
	"""
	comparisons = cf.options("two_group_comparisons") # six
	groupA_vs_groupB_FCs_Ps = []
	if len(comparisons)%6 == 0:
		for i in range(len(comparisons)/6):
			groupAs = remove_space(cf.get('two_group_comparisons', 'two_group_comparison_'+str(i+1)+'_groupA')).split(',')      # [A1, A2, A3]
			groupAname = cf.get('two_group_comparisons', 'two_group_comparison_'+str(i+1)+'_groupA_name')                          # A
			groupBs = remove_space(cf.get('two_group_comparisons', 'two_group_comparison_'+str(i+1)+'_groupB')).split(',')      # [B1, B2, B3]
			groupBname =  cf.get('two_group_comparisons', 'two_group_comparison_'+str(i+1)+'_groupB_name')                          # B
			groupA_vs_groupB_FC = cf.get('two_group_comparisons', 'two_group_comparison_'+str(i+1)+'_FC')                          # 2.0
			groupA_vs_groupB_P = cf.get('two_group_comparisons', 'two_group_comparison_'+str(i+1)+'_P')                       # 0.05
			groupA_vs_groupB_FCs_Ps.append([groupAs, groupAname, groupBs, groupBname, groupA_vs_groupB_FC, groupA_vs_groupB_P])
	return groupA_vs_groupB_FCs_Ps

#print get_2group_comparisons()

def run_cuffdiff_groupvsgroup():
	"""
	cuffdiff v2.2.1 (4237)
-----------------------------
Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]
   Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam
General Options:
  -o/--output-dir              write all output files to this directory              [ default:     ./ ]
                               写所有输出文件到这个目录                                 默认当目录
  -L/--labels                  comma-separated list of condition labels
                               逗号分割的条件标记列表
  --FDR                        False discovery rate used in testing                  [ default:   0.05 ]
                               测试中使用的错误发现率                                   默认: 0.05
  -M/--mask-file               ignore all alignment within transcripts in this file  [ default:   NULL ]
  -C/--contrast-file           Perform the constrasts specified in this file         [ default:   NULL ]
  -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]
  -u/--multi-read-correct      use 'rescue method' for multi-reads                   [ default:  FALSE ]
  -p/--num-threads             number of threads used during quantification          [ default:      1 ]
                               定量中使用的线程数                                       默认1，改成28
  --no-diff                    Don't generate differential analysis files            [ default:  FALSE ]
                               不产生差异分析文件                                       默认：产生
  --no-js-tests                Don't perform isoform switching tests                 [ default:  FALSE ]
                               不执行isoform交换检测                                   默认：执行
  -T/--time-series             treat samples as a time-series                        [ default:  FALSE ]
                               将样品看做时序的                                        默认：不看成时序的
  --library-type               Library prep used for input reads                     [ default:  below ]
                               输入reads的文库制备方法，见下
  --dispersion-method          Method used to estimate dispersion models             [ default:  below ]
                               评估离散模型的方法，见下
  --library-norm-method        Method used to normalize library sizes                [ default:  below ]
                               标准化文库大小的方法，见下

Advanced Options:
  -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]
                               平均片段长度（仅针对非配对reads）                         默认：200
  -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]
                               片段长度标准差（仅针对非配对reads）                       默认：80
  -c/--min-alignment-count     minimum number of alignments in a locus for testing   [ default:   10 ]
                               用于检验的一个基因座的最小比对数                          默认：10
  --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]
                               mle计算允许的最大迭代次数                                默认：5000
  --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:   TRUE ]
                               仅计算与参考RNA兼容的计数击中                             默认：TRUE
  --total-hits-norm            count all hits for normalization                      [ default:  FALSE ]
                               计算所有的hits用于标准化                                 默认：FALSE
  -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
                               友好日志输出（内容丰富）                                 默认：FALSE
  -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
                               安静地日志输出                                          默认：FALSE
  --seed                       value of random number generator seed                 [ default:      0 ]
                               随机数产生器的seed值                                    默认：0
  --no-update-check            do not contact server to check for update availability[ default:  FALSE ]
                               不连接服务器检查更新                                     默认连接，可去掉
  --emit-count-tables          print count tables used to fit overdispersion         [    DEPRECATED   ]
                                                                                     弃用
  --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]                               
  --num-frag-count-draws       Number of fragment generation samples                 [ default:    100 ]
  --num-frag-assign-draws      Number of fragment assignment samples per generation  [ default:     50 ]
  --max-frag-multihits         Maximum number of alignments allowed per fragment     [ default: unlim  ]
  --min-outlier-p              Min replicate p value to admit for testing            [    DEPRECATED   ]
                                                                                     弃用
  --min-reps-for-js-test       Replicates needed for relative isoform shift testing  [ default:      3 ]
  --no-effective-length-correction   No effective length correction                  [ default:  FALSE ]
  --no-length-correction       No length correction                                  [ default:  FALSE ]
  -N/--upper-quartile-norm     Deprecated, use --library-norm-method                 [    DEPRECATED   ]
                                                                                     弃用
  --geometric-norm             Deprecated, use --library-norm-method                 [    DEPRECATED   ]
                                                                                     弃用
  --raw-mapped-norm            Deprecated, use --library-norm-method                 [    DEPRECATED   ]
                                                                                     弃用
  --poisson-dispersion         Deprecated, use --dispersion-method                   [    DEPRECATED   ]
                                                                                     弃用

Debugging use only:
  --read-skip-fraction         Skip a random subset of reads this size               [ default:    0.0 ]
  --no-read-pairs              Break all read pairs                                  [ default:  FALSE ]
  --trim-read-length           Trim reads to be this long (keep 5' end)              [ default:   none ]
  --no-scv-correction          Disable SCV correction                                [ default:  FALSE ]

Supported library types:
        ff-firststrand
        ff-secondstrand
        ff-unstranded
        fr-firststrand         我们用的这个？待lw确认
        fr-secondstrand
        fr-unstranded (default)
        transfrags



ff-firststrand
ff-secondstrand
ff-unstranded
fr-firststrand	dUTP, NSR, NNSR	与上类似，除了片段最右端的首先测序（或者仅测序单端）。相当于仅合成的第一条链被测序
fr-secondstrand	Directional Illumina (Ligation), Standard SOLiD	与上雷同，除了片段最左端的首先测序（或仅测序单端）。相当于仅合成的第一条链被测序
fr-unstranded (default)	Standard Illumina	来自片段最左端的（转录坐标）比对到转录本链，而最右端的比对到相反链
transfrags  


Supported dispersion methods:
        blind
        per-condition
        poisson
        pooled (default)

blind	All samples are treated as replicates of a single global “condition” and used to build one model.
        所有样品认为是一个“全局条件”的重复，并被用来构建一个模型
per-condition	Each replicated condition receives its own model. Only available when all conditions have replicates.
                每个重复条件使用它自身的模型。仅所有条件都有重复时才有用
poisson	
pooled	Each replicated condition is used to build a model, then these models are averaged to provide a single global model for all conditions in the experiment. (Default)
        使用每个重复条件构建模型，然后这些模型平均以提供单个全局模型用于试验中的所有条件（默认）

Supported library normalization methods:
支持的标准化方法：
        classic-fpkm
        库大小因子设置成1。即对FPKM值或片段数不设置scaling（Cufflinks默认的）。我们使用的也是这个。
        geometric (default)
        FPKMs and fragment counts are scaled via the median of the geometric means of fragment counts across all libraries, as described in Anders and Huber (Genome Biology, 2010). This policy identical to the one used by DESeq. (default for Cuffdiff)
        quartile
        FPKMs and fragment counts are scaled via the ratio of the 75 quartile fragment counts to the average 75 quartile value across all libraries.       
	"""
	samples  = getfastq()
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')     # Homo_sapiens_HG19.gtf
	hisat2_fasta = cf.get('mapping_parameters','hisat2_fasta')
	comparisons = get_2group_comparisons()
	if comparisons:
		for comparison in comparisons:
			groupA, groupAname, groupB, groupBname, fc, p = comparison # fc, p not used here
			foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
			try:
				shutil.rmtree(foldname)
				os.mkdir(foldname)
			except:
				os.mkdir(foldname)
			os.chdir(foldname)
			groupAcxb = ','.join(['../'+i.split('_sequence.fastq')[0]+'/abundances.cxb' for i in groupA])
			groupBcxb = ','.join(['../'+i.split('_sequence.fastq')[0]+'/abundances.cxb' for i in groupB])
			labels = groupAname+','+groupBname
			inputcxb = groupAcxb+' '+groupBcxb
			# for firststrand
			# lin's 
			#nohup'/opt/cufflinks2/cufflinks-2.1.1.Linux_x86_64/cuffdiff --library-type fr-firststrand --no-update-check -u -v --emit-count-tables -p 24 -L 6W-KO-FM-2,6W-WT-FM-1,6W-KO-M-2,6W-WT-M-1 -b /public/INDEX_GENOME/BOWTIE2_UCSC_MM10/BOWTIE2_UCSC_MM10.fa /public/Annotation/ENSEMBL_MM10/Mus_musculus_MM10.gtf /data/ArrayStar-Seq/Analysis_Results_2015_05_14_10_21_17/ozm_mRNA_Mus_musculus/tmp/ozm_mRNA_Mus_musculus_6W-KO-FM-2/accepted_hits.bam /data/ArrayStar-Seq/Analysis_Results_2015_05_14_10_21_17/ozm_mRNA_Mus_musculus/tmp/ozm_mRNA_Mus_musculus_6W-WT-FM-1/accepted_hits.bam /data/ArrayStar-Seq/Analysis_Results_2015_05_14_10_21_17/ozm_mRNA_Mus_musculus/tmp/ozm_mRNA_Mus_musculus_6W-KO-M-2/accepted_hits.bam /data/ArrayStar-Seq/Analysis_Results_2015_05_14_10_21_17/ozm_mRNA_Mus_musculus/tmp/ozm_mRNA_Mus_musculus_6W-WT-M-1/accepted_hits.bam &
			#cmd = '/opt/cufflinks2/cufflinks-2.1.1.Linux_x86_64/cuffdiff --library-norm-method classic-fpkm --library-type fr-firststrand --no-update-check -u -v --emit-count-tables -p 16 -L %s -b %s %s %s' % (inputname, bowtie2_fasta, tophat2_gtf, inputbam)
			
			# here use classic-fpkm to keep the expression and diff results the same
			# change library-type if needed
			cmd = '/root/cufflinks-2.2.1.Linux_x86_64/cuffdiff --library-norm-method classic-fpkm --library-type fr-firststrand --no-update-check -u -v --emit-count-tables -p 20 -L %s -b %s %s %s' % (labels, hisat2_fasta, hisat2_gtf, inputcxb)
			print cmd
			os.system(cmd)
			os.chdir(global_now_location)
	else:
		print '!!!!!!!!!! There is no group comparison'

###################################### 15 run cuffdiff for two groups, over

 
###################################### 16 parse gtf file, extract geneid, transcriptid, biotype and strand,....

def index1(lists, sample):
	"""
	output the index of sample in lists
	"""
	idx =[]
	for i in range(len(lists)):
		if lists[i] == sample:
			idx.append(i)
	return idx

def index2(lists, samplelist):
	"""
	output the indexs of sample list in lists
	"""
	idx = []
	for j in range(len(lists)):
		for s in samplelist:
			if lists[j] == s:
				idx.append(j)
	return idx


def xbiotype(lists):
	idx =[]
	for i in range(len(lists)):
		if 'gene_biotype' in lists[i]:
			idx.append(i)
	if idx:
		biotype = lists[idx[0]].split('"')[1]
	else:
		biotype= 'long_noncoding'
	return biotype



def gene_id_index(list):
	return [i for i in range(len(list)) if 'gene_id' in list[i]]


def parse_gtf_gene(file):  # to test
	"""
	some annotation information, indexed by ensembl id
	['3prime_overlapping_ncrna','antisense','lincRNA','retained_intron','sense_intronic','processed_transcript','misc_RNA', 'long_noncoding']:
	"""
	f = open(file)
	# 1		2				3		4		5		6	7	8	9
	#chr1	protein_coding	gene	171669300	171711387	.	-	.	gene_id "ENSG00000117533"; gene_name "VAMP4"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
	#chr1	protein_coding	transcript	171669300	171711387	.	-	.	gene_id "ENSG00000117533"; transcript_id "ENST00000236192"; gene_name "VAMP4"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "VAMP4-004"; transcript_source "ensembl_havana"; tag "CCDS"; ccds_id "CCDS1298";
	#chr1	long_noncoding	exon	171669296	171671223	.	-	.	gene_id "VAMP4"; transcript_id "uc001ghv.2"; exon_number "4"; gene_name "VAMP4";

	genes = {}
	for line in f:
		if not line.startswith('#'):
			if 'gene_id' in line:
				linelist = line.strip('\n').split('\t')
				annotation_parts = linelist[8].split(';')
				geneid = annotation_parts[gene_id_index(annotation_parts)[0]].split('"')[1]
				biotype = xbiotype(annotation_parts)
				strand = linelist[6]
				genes[geneid]=biotype+'\t'+strand
	f.close()
	return genes


def transcript_id_index(list):
	return [i for i in range(len(list)) if 'transcript_id' in list[i]]


def parse_gtf_transcript(file):  
	"""
	some annotation information, indexed by ensembl id
	
	"""
	f = open(file)
	# 1		2				3		4		5		6	7	8	9
	# chrX	protein_coding	exon	3076920	3078823	.	+	.	 gene_id "ENSMUSG00000095801"; transcript_id "ENSMUST00000178827"; exon_number "1"; gene_name "GMCL1P1"; gene_biotype "protein_coding"; transcript_name "GMCL1P1.1-001"; exon_id "ENSMUSE00001059231";
	transcripts = {}
	for line in f:
		if not line.startswith('#'):
			if 'transcript_id' in line:
				linelist = line.strip('\n').split('\t')
				annotation_parts = linelist[8].split(';')
				#print annotation_parts
				transcriptid = annotation_parts[transcript_id_index(annotation_parts)[0]].split('"')[1]
				biotype = xbiotype(annotation_parts)  #for lncRNA added
				strand = linelist[6]
				transcripts[transcriptid]= biotype+'\t'+strand
	f.close()
	return transcripts
###################################### 16 parse gtf file, extract geneid, transcriptid, biotype and strand,....

###################################### 17 gene expression results (txt) fpkm

# add lncRNA-associated genes (transcript level lncRNA)
def read_associated(file):
	"""
	transID	name	chrom	strand	txStart	txEnd	blockCount	blockSizes	qStarts	tStarts	GeneName	Source	sourceID	RNAlength	relationship	associated_gene_transID	association_gene_acc	association_gene_name	associated_gene_strand	associated_gene_start	associated_gene_end
	"""
	f = open(file)
	head = f.readline()
	d = {}
	for line in f:
		linelist = line.strip('\n').split('\t')
		# transID: Source, RNAlength	relationship	associated_gene_transID	association_gene_acc	association_gene_name	associated_gene_strand	associated_gene_start	associated_gene_end
		d.setdefault(linelist[0], [ ]).append(linelist[11]+'\t'+'\t'.join(linelist[13:]))
	f.close()
	return d

def lncRNAaddgene(file, associatedfile):
	dat = read_associated(associatedfile)
	f = open(file)
	f2 = open(file.replace('lnctmp','lnc'),'w')
	head = f.readline().strip()
	"""
	# 1           2       3        4      5        6      7             8          9       10     11         12        13            14         15                16                      17                          18                     19                     20                      21
	#transID	name	chrom	strand	txStart	txEnd	blockCount	blockSizes	qStarts	tStarts	GeneName	Source	sourceID	RNAlength	relationship	associated_gene_transID	association_gene_acc	association_gene_name	associated_gene_strand	associated_gene_start	associated_gene_end
	AK022045	AK022045	chr10	+	24536053	24544975	2	1267,925,	7,1274,	24536053,24544050,	PRINS	lncRNAdisease	6	2199	exon sense-overlapping	ENST00000430453	ENST00000430453	KIAA1217	+	24544256	24820821
	AK022045	AK022045	chr10	+	24536053	24544975	2	1267,925,	7,1274,	24536053,24544050,	PRINS	lncRNAdisease	6	2199	intron sense-overlapping	ENST00000376456	ENST00000376456	KIAA1217	+	24497719	24816871
	"""
	print >>f2, head+'\tLncRNA_Source\tLncRNA_length\trelationship\tassociated_gene_transID\tassociation_gene_acc\tassociation_gene_name\tassociated_gene_strand\tassociated_gene_start\tassociated_gene_end'
	for line in f:
		line = line.strip('\n')
		linelist=line.split('\t')
		try:
			for one in dat[linelist[0]]:  # transID
				print >>f2, line+'\t'+one
		except:
			print >>f2, line+'\tEnsembl\t\t\t\t\t\t\t\t'
	f.close()
	f2.close()
# add lncRNA (transcript level) -associated genes, over


###################################### 17 gene expression results (txt) fpkm
# add mRNA information from NCBI FTP, pathway, go infomation (see add_gene_pathway_go.py)
def read_geneinfo(file):
	"""
	note: one gene one info, although there are some replicates
	tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\ttype_of_gene\tSymbol_from_nomenclature_authority\tFull_name_from_nomenclature_authority\tNomenclature_status\tOther_designations\tModification_date\tPathway\tGO-CC\tGO-MF\tGO-BP
	"""
	f = open(file)
	head = f.readline()
	d = {}
	for line in f:
		linelist = line.strip('\n').split('\t')
		# tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\ttype_of_gene\tSymbol_from_nomenclature_authority\tFull_name_from_nomenclature_authority\tNomenclature_status\tOther_designations\tModification_date\tPathway\tGO-CC\tGO-MF\tGO-BP
		d[linelist[2].upper()] = linelist[1]+'\t'+linelist[4]+'\t'+linelist[5]+'\t'+linelist[6]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[15]+'\t'+linelist[16]+'\t'+linelist[17]+'\t'+linelist[18]
	f.close()
	return d

def addgeneinfo(file, geneinfofile):
	dat = read_geneinfo(geneinfofile)
	f = open(file)
	f2 = open(file.replace('protein_codingtmp','protein_coding'),'w')
	head = f.readline().strip('\n')
	"""
	# 1           2       3        4      5        6      7             8          9       10     11         12        13            14         15                16                      17                          18                     19                     20                      21
	#transID	name	chrom	strand	txStart	txEnd	blockCount	blockSizes	qStarts	tStarts	GeneName	Source	sourceID	RNAlength	relationship	associated_gene_transID	association_gene_acc	association_gene_name	associated_gene_strand	associated_gene_start	associated_gene_end
	AK022045	AK022045	chr10	+	24536053	24544975	2	1267,925,	7,1274,	24536053,24544050,	PRINS	lncRNAdisease	6	2199	exon sense-overlapping	ENST00000430453	ENST00000430453	KIAA1217	+	24544256	24820821
	AK022045	AK022045	chr10	+	24536053	24544975	2	1267,925,	7,1274,	24536053,24544050,	PRINS	lncRNAdisease	6	2199	intron sense-overlapping	ENST00000376456	ENST00000376456	KIAA1217	+	24497719	24816871
	"""
	print >>f2, head+'\tGeneID\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\tpathway\tGO-CC\tGO-MF\tGO-BP'
	for line in f:
		line = line.strip('\n')
		linelist=line.split('\t')
		try:
			print >>f2, line+'\t'+dat[linelist[0].upper()]
		except:
			print >>f2, line+'\t\t\t\t\t\t\t\t\t\t'
	f.close()
	f2.close()
# add gene info, over


def get_gene_expression(): # 1, all gene level, 2, only protein_coding gene level
	"""
	gene expression files, both fpkm and read count
	genes.fpkm_tracking		tracking_id	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	6W-WT-FM-1_FPKM	6W-WT-FM-1_conf_lo	6W-WT-FM-1_conf_hi	6W-WT-FM-1_status	6W-KO-FM-2_FPKM	6W-KO-FM-2_conf_lo	6W-KO-FM-2_conf_hi	6W-KO-FM-2_status	6W-KO-M-2_FPKM	6W-KO-M-2_conf_lo	6W-KO-M-2_conf_hi	6W-KO-M-2_status	6W-WT-M-1_FPKM	6W-WT-M-1_conf_lo	6W-WT-M-1_conf_hi	6W-WT-M-1_status
	genes.count_tracking	tracking_id	6W-WT-FM-1_count	6W-WT-FM-1_count_variance	6W-WT-FM-1_count_uncertainty_var	6W-WT-FM-1_count_dispersion_var	6W-WT-FM-1_status	6W-KO-FM-2_count	6W-KO-FM-2_count_variance	6W-KO-FM-2_count_uncertainty_var	6W-KO-FM-2_count_dispersion_var	6W-KO-FM-2_status	6W-KO-M-2_count	6W-KO-M-2_count_variance	6W-KO-M-2_count_uncertainty_var	6W-KO-M-2_count_dispersion_var	6W-KO-M-2_status	6W-WT-M-1_count	6W-WT-M-1_count_variance	6W-WT-M-1_count_uncertainty_var	6W-WT-M-1_count_dispersion_var	6W-WT-M-1_status
	the two files are of same line numbers
	mRNA is protein_coding, add gene info to mRNA
	"""
	samples  = [sample.split('_sequence.fastq')[0] for sample in getfastq()]
	fpkm = open('./cuffdiff/genes.fpkm_tracking') # genes
	fpkmlines = fpkm.readlines()
	fpkm.close()
	fcount = open('./cuffdiff/genes.count_tracking')
	fcountlines = fcount.readlines()
	fcount.close()
	headfpkm = fpkmlines[0].strip('\n').split('\t')
	#print headfpkm
	headfcount = fcountlines[0].strip('\n').split('\t')
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	gene_biotype_strand_dict = parse_gtf_gene(hisat2_gtf)
	f = open('./cuffdiff/Gene_Level_Profile.txt','w')
	fprotein = open('./cuffdiff/Gene_Level_Profile_protein_codingtmp.txt','w')  # only protein coding
	# GeneName (4),ensembl gene id (3), locus (6), chr17:74449432-74466199,  # name_FPKMs
	sample_index_fpkm = index2(headfpkm, [i+'_FPKM' for i in samples])
	sample_index_fcount = index2(headfcount, [i+'_count' for i in samples])
	#print sample_index_fpkm
	print >>f, headfpkm[4]+'\t'+headfpkm[3]+'\t'+'biotype'+'\t'+'strand'+'\t'+headfpkm[6]+'\t'+ 'Length'+'\t'+'\t'.join([headfpkm[i] for i in sample_index_fpkm])+'\t'+'\t'.join([headfcount[i] for i in sample_index_fcount])
	print >>fprotein, headfpkm[4]+'\t'+headfpkm[3]+'\t'+'biotype'+'\t'+'strand'+'\t'+headfpkm[6]+'\t'+ 'Length'+'\t'+'\t'.join([headfpkm[i] for i in sample_index_fpkm])+'\t'+'\t'.join([headfcount[i] for i in sample_index_fcount])
	for linenum in range(1, len(fpkmlines)):
		fpkmlist = fpkmlines[linenum].strip('\n').split('\t')
		start = int(fpkmlist[6].split(':')[1].split('-')[0])
		end = int(fpkmlist[6].split(':')[1].split('-')[1])
		length = str(end-start)
		fcountlist = fcountlines[linenum].strip('\n').split('\t')
		#print gene_biotype_strand_dict[fpkmlist[3]]
		#print >>f, fpkmlist[4]+'\t'+fpkmlist[3]+'\t'+gene_biotype_strand_dict[fpkmlist[3]]+'\t'+fpkmlist[6]+'\t'+ length+'\t'+'\t'.join([fpkmlist[i] for i in sample_index_fpkm])+'\t'+'\t'.join([fcountlist[i] for i in sample_index_fcount])
		print >>f, fpkmlist[4]+'\t'+fpkmlist[3]+'\t'+gene_biotype_strand_dict[fpkmlist[3]]+'\t'+fpkmlist[6]+'\t'+ length+'\t'+'\t'.join([fpkmlist[i] for i in sample_index_fpkm])+'\t'+'\t'.join([fcountlist[i] for i in sample_index_fcount])
		# protein_coding
		if gene_biotype_strand_dict[fpkmlist[3]].split()[0] == 'protein_coding':
			print >>fprotein, fpkmlist[4]+'\t'+fpkmlist[3]+'\t'+gene_biotype_strand_dict[fpkmlist[3]]+'\t'+fpkmlist[6]+'\t'+ length+'\t'+'\t'.join([fpkmlist[i] for i in sample_index_fpkm])+'\t'+'\t'.join([fcountlist[i] for i in sample_index_fcount])
	f.close()
	fprotein.close()
	
	# add mRNA gene info
	addgeneinfo('./cuffdiff/Gene_Level_Profile_protein_codingtmp.txt', geneinfofile)


###################################### 17 gene expression results (txt) fpkm, over


###################################### 17 isoform expression results (txt), fpkm

def get_isoform_expression(): # 1, all transcript level, 2, only lncRNA transcript level
	"""
	isoforms expression files, both fpkm and read count
	isoforms.fpkm_tracking	tracking_id	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	6W-WT-FM-1_FPKM	6W-WT-FM-1_conf_lo	6W-WT-FM-1_conf_hi	6W-WT-FM-1_status	6W-KO-FM-2_FPKM	6W-KO-FM-2_conf_lo	6W-KO-FM-2_conf_hi	6W-KO-FM-2_status	6W-KO-M-2_FPKM	6W-KO-M-2_conf_lo	6W-KO-M-2_conf_hi	6W-KO-M-2_status	6W-WT-M-1_FPKM	6W-WT-M-1_conf_lo	6W-WT-M-1_conf_hi	6W-WT-M-1_status
	isoforms.count_tracking	tracking_id	6W-WT-FM-1_count	6W-WT-FM-1_count_variance	6W-WT-FM-1_count_uncertainty_var	6W-WT-FM-1_count_dispersion_var	6W-WT-FM-1_status	6W-KO-FM-2_count	6W-KO-FM-2_count_variance	6W-KO-FM-2_count_uncertainty_var	6W-KO-FM-2_count_dispersion_var	6W-KO-FM-2_status	6W-KO-M-2_count	6W-KO-M-2_count_variance	6W-KO-M-2_count_uncertainty_var	6W-KO-M-2_count_dispersion_var	6W-KO-M-2_status	6W-WT-M-1_count	6W-WT-M-1_count_variance	6W-WT-M-1_count_uncertainty_var	6W-WT-M-1_count_dispersion_var	6W-WT-M-1_status
	the two files are of same line numbers
	"""
	samples  = [sample.split('_sequence.fastq')[0] for sample in getfastq()]
	fpkm = open('./cuffdiff/isoforms.fpkm_tracking')
	fpkmlines = fpkm.readlines()
	fpkm.close()
	fcount = open('./cuffdiff/isoforms.count_tracking')
	fcountlines = fcount.readlines()
	fcount.close()
	headfpkm = fpkmlines[0].strip('\n').split('\t')
	#print headfpkm
	headfcount = fcountlines[0].strip('\n').split('\t')
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	transcript_biotype_strand_dict = parse_gtf_transcript(hisat2_gtf)
	f = open('./cuffdiff/Transcript_Level_Profile.txt','w')
	flnc = open('./cuffdiff/Transcript_Level_Profile_lnctmp.txt','w')  # only lnc transcript
	# GeneName (4), ensemble Transcriptid (3), locus (6), length (7),  # name_FPKMs
	sample_index_fpkm = index2(headfpkm, [i+'_FPKM' for i in samples])
	sample_index_fcount = index2(headfcount, [i+'_count' for i in samples])
	#print sample_index_fpkm
	print >>f, headfpkm[0].replace('tracking_id','transcript_id')+'\t'+headfpkm[4]+'\t'+headfpkm[3]+'\t'+'biotype\tstrand'+'\t'+headfpkm[6]+'\t'+ 'Length'+'\t'+'\t'.join([headfpkm[i] for i in sample_index_fpkm])+'\t'+'\t'.join([headfcount[i] for i in sample_index_fcount])
	print >>flnc, headfpkm[0].replace('tracking_id','transcript_id')+'\t'+headfpkm[4]+'\t'+headfpkm[3]+'\t'+'biotype\tstrand'+'\t'+headfpkm[6]+'\t'+ 'Length'+'\t'+'\t'.join([headfpkm[i] for i in sample_index_fpkm])+'\t'+'\t'.join([headfcount[i] for i in sample_index_fcount])
	for linenum in range(1, len(fpkmlines)):
		fpkmlist = fpkmlines[linenum].strip('\n').split('\t')
		start = int(fpkmlist[6].split(':')[1].split('-')[0])
		end = int(fpkmlist[6].split(':')[1].split('-')[1])
		length = str(end-start)
		fcountlist = fcountlines[linenum].strip('\n').split('\t')
		print >>f, fpkmlist[0]+'\t'+fpkmlist[4]+'\t'+fpkmlist[3]+'\t'+transcript_biotype_strand_dict[fpkmlist[0]]+'\t'+fpkmlist[6]+'\t'+ length+'\t'+'\t'.join([fpkmlist[i] for i in sample_index_fpkm])+'\t'+'\t'.join([fcountlist[i] for i in sample_index_fcount])
		if transcript_biotype_strand_dict[fpkmlist[0]].split()[0] in lncRNAclasses:
			print >>flnc, fpkmlist[0]+'\t'+fpkmlist[4]+'\t'+fpkmlist[3]+'\t'+transcript_biotype_strand_dict[fpkmlist[0]]+'\t'+fpkmlist[6]+'\t'+ length+'\t'+'\t'.join([fpkmlist[i] for i in sample_index_fpkm])+'\t'+'\t'.join([fcountlist[i] for i in sample_index_fcount])
	f.close()
	flnc.close()
	lncRNAaddgene('./cuffdiff/Transcript_Level_Profile_lnctmp.txt', lnc_associated_genefile)


###################################### 17 isoform expression results (txt), fpkm, over


###################################### 18 gene expression excel (fpkm), start

def read_data_content(file):
	"""
	read txt file into list [[a...], [b...], ..., []]
	"""
	f = open(file)
	alls = []
	for line in f:
		linelist = line.strip('\n').split('\t')
		alls.append(linelist)
	f.close()
	return alls

def formats(st):
	if '.' in st or 'e' in st:
		try:
			return float(st)
		except:
			return st
	else:
		try:
			return int(st)
		except:
			return st

def read_data_content_with_float(file):
	"""
	read txt file into list [[a...], [b...], ..., []]
	if there are float numbers, return float
	"""
	f = open(file)
	alls = []
	for line in f:
		linelist = line.strip(os.linesep).split('\t')
		alls.append([formats(i) for i in linelist])
	f.close()
	return alls

def excel_expression_all_samples_gene(): # 1, all, 2, protein_coding, only gene level
	"""
	# file is Gene Expression.xlsx, Gene Expression_protein_coding.xlsx
	gene: /home/workplace/RNA-seq/test/cuffdiff/Gene_Level_Profile.txt
 
1, merge range
worksheet.merge_range(first_row, first_col, last_row, last_col, data, cell_format)

2, string with url (not display)
link_format = workbook.add_format({'color': 'blue', 'underline': 1})
worksheet.write_url('A1', 'http://www.python.org', link_format, 'Python')

3, format within a cell

italic = workbook.add_format({'italic': True})

worksheet.write_rich_string('A1','This is ',
                            bold, 'bold',
                            ' and this is ',
                            italic, 'italic')

4, insert image  (PNG, JPEG or BMP)
worksheet.insert_image('B2', 'python.png')

5, cell format
format = workbook.add_format({'bold': True})

{'bold':     True,
'vjustify': True,
'border':   6,
'italic': True,
'align':    'center',
'valign':   'vcenter',
'fg_color': '#D7E4BC',
'align':    'center',

}

6, xlsx file properties
workbook.set_properties({
    'title':    'MeDIP-Seq Signal',
    'subject':  '',
    'author':   'jimmy',
    'manager':  '',
    'company':  'KC',
    'category': 'MeDIP-Seq',
    'keywords': 'MeDIP-Seq, Sequencing',
    'comments': ''})

7, background color
format.set_pattern(1)  # This is optional when using a solid fill.
format.set_bg_color('green')
 
	"""
	samplenumber = len(getfastq())
	excel_title_info_gene = """Gene Expression for all samples

Column A: gene_short_name, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, including protein_coding, pseudogene, antisense, ....
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F: Length, the length of the locus. It equals to end-start of the locus.
Column G ~ %s: FPKM values of all samples of the gene.
Column %s ~ %s: counts of all samples of the gene."""
	
	# write into excel
	get_gene_expression()  # all & protein_coding
	################ for all
	workbook = Workbook('Gene Expression (protein+lnc).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align':    'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	data = read_data_content_with_float('cuffdiff/Gene_Level_Profile.txt')  # if not change this
	#print len(data)
	sheetname = 'Gene Expression'
	worksheet = workbook.add_worksheet(sheetname)  # sheetname
	colnum = len(data[0])
	worksheet.merge_range(0, 0, 19+info_addline, colnum-1, '', annotationformat)
	worksheet.write_rich_string(0, 0, red, "'", excel_title_info_gene % (colalphabet[5+samplenumber], colalphabet[5+samplenumber+1],colalphabet[5+samplenumber*2]), annotationformat)
	data_title = data[0]
	rownum = len(data)
	for i in range(colnum):
		worksheet.write(0+22+info_addline,i,data_title[i], titleformat)
	for row in range(1, rownum):
		for col in range(colnum):
			worksheet.write(row+22+info_addline, col, data[row][col],zhengwenformat)
	workbook.close()
	
	############### for protein coding
	excel_title_info_gene_protein_coding = """mRNA Expression for all samples

Column A: gene_short_name, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, protein_coding.
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F: Length, the length of the locus. It equals to end-start of the locus.
Column G ~ %s: FPKM values of all samples of the gene.
Column %s ~ %s: counts of all samples of the gene
Column %s ~ %s: annotation information, including GeneID, Synonyms, dbXrefs (other database ID), chromosome, map_location, description, pathway and GO annotations.

note: since the different sources, some genes without annotation information"""

	# write into excel
	workbook2 = Workbook('mRNA Expression Profiling.xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary to see
	annotationformat2 =  workbook2.add_format({'align': 'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	pkkmformat2 = workbook2.add_format({'align':'center', 'bold':True, 'fg_color': '#FABF8F', 'font_name':'Times New Roman'}) # add 20170325
	countformat2= workbook2.add_format({'align':'center', 'bold':True, 'fg_color': '#92CDDC', 'font_name':'Times New Roman'}) # add 20170325
	annoformat2 = workbook2.add_format({'align':'center', 'bold':True, 'fg_color': '#B1A0C7', 'font_name':'Times New Roman'}) # add 20170325
	data2 = read_data_content_with_float('cuffdiff/Gene_Level_Profile_protein_coding.txt')  
	sheetname2 = 'mRNA Expression'
	worksheet2 = workbook2.add_worksheet(sheetname2)  # sheetname
	colnum2 = len(data2[0])
	worksheet2.merge_range(0, 0, 19+info_addline, colnum2-1, '', annotationformat2)
	worksheet2.write_rich_string(0, 0, red2, "'", excel_title_info_gene_protein_coding % (colalphabet[5+samplenumber], colalphabet[5+samplenumber+1],colalphabet[5+samplenumber*2], colalphabet[6+samplenumber*2], colalphabet[15+samplenumber*2]), annotationformat2)
	data_title2 = data2[0]
	rownum2 = len(data2)
	worksheet2.merge_range(21+info_addline, 6,                21+info_addline, 6+samplenumber-1, 'FPKM Values', pkkmformat2)  # 6, 14
	worksheet2.merge_range(21+info_addline, 6+samplenumber,   21+info_addline, 6+samplenumber*2-1, 'Counts',      countformat2) # 15, 23
	worksheet2.merge_range(21+info_addline, 6+samplenumber*2, 21+info_addline, colnum2-1,        'Annotations', annoformat2) # 24~
	for i2 in range(colnum2):
		worksheet2.write(0+22+info_addline,i2,data_title2[i2], titleformat2)
	for row2 in range(1, rownum2):
		for col2 in range(colnum2):
			worksheet2.write(row2+22+info_addline, col2, data2[row2][col2],zhengwenformat2)
	workbook2.close()
 

###################################### 18 gene expression excel (fpkm), over

###################################### 19 isoform expression excel (fpkm), start


def excel_expression_all_samples_transcript(): #1, all, 2 lncRNA transcript level (add associated)
	"""
	# file is Transcript_Level_Profile.txt
	# write into excel Transcript Expression.xlsx & lncRNA transcript
	"""
	samplenumber = len(getfastq())
	excel_title_info_transcript = """Transcript Expression for all samples

Column A: Transcript_id, the Ensembl transcript identifier.
Column B: gene_short_name, the gene name.
Column C: gene_id, the Ensembl gene identifier.
Column D: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column E: strand, transcription direction.
Column F: locus, genomic locus. chromosome: start-end.
Column G: Length, the length of the gene. It equals to end-start of the locus.
Column H ~ %s: FPKM values of all samples of the transcript.
Column %s ~ %s: counts of all samples of the transcript."""

	
	get_isoform_expression()  # 2 files generated
	# write into excel 
	workbook = Workbook('Transcript Expression(protein+LncRNA).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align':    'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	data = read_data_content_with_float('cuffdiff/Transcript_Level_Profile.txt')  # if not change this
	#print len(data)
	sheetname = 'Transcript Expression'
	worksheet = workbook.add_worksheet(sheetname)  # sheetname
	colnum = len(data[0])
	worksheet.merge_range(0, 0, 19+info_addline, colnum-1, '', annotationformat)
	worksheet.write_rich_string(0, 0, red, "'",excel_title_info_transcript % (colalphabet[6+samplenumber], colalphabet[6+samplenumber+1],colalphabet[6+samplenumber*2]), annotationformat)
	data_title = data[0]
	rownum = len(data)
	for i in range(colnum):
		worksheet.write(0+22+info_addline,i,data_title[i], titleformat)
	for row in range(1, rownum):
		for col in range(colnum):
			worksheet.write(row+22+info_addline, col, data[row][col],zhengwenformat)
	workbook.close()
	# for lncRNA transcript level
	excel_title_info_lnc_transcript = """Transcript level LncRNA expression for all samples

Column A: Transcript_id, the Ensembl transcript identifier.
Column B: gene_short_name, the gene name.
Column C: gene_id, the Ensembl gene identifier.
Column D: biotype, including 3prime_overlapping_ncrna, antisense, lincRNA, retained_intron, sense_intronic, processed_transcript, misc_RNA, long_noncoding.
Column E: strand, transcription direction.
Column F: locus, genomic locus. chromosome: start-end.
Column G: Length, the length of the gene. It equals to end-start of the locus.
Column H ~ %s: FPKM values of all samples of the transcript.
Column %s ~ %s: counts of all samples of the transcript.
Column %s ~ %s: the genomic organization of LncRNAs and the information of the associated coding genes.

the classes of LncRNAs defined by genomic organization:
"exon sense-overlapping": the LncRNA's exon is overlapping a coding transcript exon on the same genomic strand;
"intron sense-overlapping": the LncRNA is overlapping the intron of a coding transcript on the same genomic strand;
"intronic antisense": the LncRNA is overlapping the intron of a coding transcript on the antisense strand;
"natural antisense": the LncRNA is transcribed from the antisense strand and overlapping with a coding transcript; 
"bidirectional": the LncRNA is oriented head to head to a coding transcript within 1000 bp; 
"intergenic": there are no overlapping or bidirectional coding transcripts nearby the LncRNA."""

	# write into excel 
	workbook2 = Workbook('LncRNA Expression Profiling.xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 = workbook2.add_format({'align':'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2=  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	pkkmformat2 = workbook2.add_format({'align':'center', 'bold':True, 'fg_color': '#FABF8F', 'font_name':'Times New Roman'}) # add 20170325
	countformat2= workbook2.add_format({'align':'center', 'bold':True, 'fg_color': '#92CDDC', 'font_name':'Times New Roman'}) # add 20170325
	annoformat2 = workbook2.add_format({'align':'center', 'bold':True, 'fg_color': '#B1A0C7', 'font_name':'Times New Roman'}) # add 20170325
	data2 = read_data_content_with_float('cuffdiff/Transcript_Level_Profile_lnc.txt')  #lnctmp --> lnc
	#print len(data)
	sheetname2 = 'LncRNA Expression'
	worksheet2 = workbook2.add_worksheet(sheetname2[:30]) # sheetname
	colnum2 = len(data2[0])
	worksheet2.merge_range(0, 0, 19+info_addline, colnum2-1, '', annotationformat2)
	worksheet2.write_rich_string(0, 0, red2, "'",excel_title_info_lnc_transcript % (colalphabet[6+samplenumber], colalphabet[6+samplenumber+1],colalphabet[6+samplenumber*2], colalphabet[6+samplenumber*2+1], colalphabet[6+samplenumber*2+1+8]), annotationformat2)
	data_title2 = data2[0]
	rownum2 = len(data2)
	worksheet2.merge_range(21+info_addline, 7,                21+info_addline, 7+samplenumber-1, 'FPKM Values', pkkmformat2)  # 7, 15
	worksheet2.merge_range(21+info_addline, 7+samplenumber,   21+info_addline, 7+samplenumber*2-1, 'Counts',      countformat2) # 16, 24
	worksheet2.merge_range(21+info_addline, 7+samplenumber*2, 21+info_addline, colnum2-1,        'Annotations', annoformat2) # 25~
	for i2 in range(colnum2):
		worksheet2.write(0+22+info_addline,i2,data_title2[i2], titleformat2)
	for row2 in range(1, rownum2):
		for col2 in range(colnum2):
			worksheet2.write(row2+22+info_addline, col2, data2[row2][col2],zhengwenformat2)
	workbook2.close()

###################################### 19 isoform expression excel (fpkm), over


###################################### 20 two sample all comparison gene (txt), start

def get_2_samples_comparisons():
	"""
	output 2 samples comparison
	"""
	
	comparisons  = cf.options("two_sample_comparisons")
	complist = []
	for comp in comparisons:
		c = [x.strip().replace('_sequence.fastq','') for x in cf.get('two_sample_comparisons',comp).split(',')]
		complist.append(c)
	return complist
	
	
def log2fc2fc(logfc):
	"""
	logfc is string and inverted, so here we must change + into -, - into +,  0 into 0
	change log2fc into fc, just 2**
	"""
	if logfc == '-inf':
		fc = 'inf'
	elif logfc == 'inf':
		fc= '-inf'
	elif float(logfc)>0:
		try:
			fc=-2**float(logfc)
		except:
			fc = '-inf'
	elif float(logfc)<0:
		try:
			fc=1.0/(2**(float(logfc)))
		except:
			fc = 'inf'
	else:
		fc=1.0
	return fc

def two_sample_comparison_gene(): # 1, protein+lnc gene level, 2, protein gene level
	"""
	test_id				gene_id				gene	locus					sample_1	sample_2	status	value_1		value_2	log2(fold_change)	test_stat	p_value	q_value	significant
	ENSMUSG00000000058	ENSMUSG00000000058	Cav2	chr6:17197750-17385604	6W-WT-FM-1	6W-KO-FM-2	OK		0.245057	14.4872	5.88552				1.05052		0.20925	0.366351	no
	sample_1 (6W-WT-FM-1) 
	sample_2 (6W-KO-FM-2)
	value_1  (0.245057)
	value_2  (14.4872)
	
	sample_2 vs sample_1 (6W-KO-FM-2 vs 6W-WT-FM-1)
	14.4872/0.245057 = 59.11767
	2**5.88552 = 59.11767
	
	The name order (A, B) was the same to value order (A, B), but the logfc is opposite
	so if you want to keep the logfc, the sample name order should be reversed,
	or just add minus in front of the logfc, we use this
	
	gene_exp.diff.6W-KO-M-2_vs_6W-KO-FM-2.txt
	"""
	
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	gene_biotype_strand_dict = parse_gtf_gene(hisat2_gtf)
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		f =open(foldname+'/gene_exp.diff') # all comparisons
		f2 = open(foldname+'/gene_exp.diff.%s_vs_%s.txt' % (samplea, sampleb),'w') # all for A vs B, no need to add associated genes
		fprotein = open(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.txt' % (samplea, sampleb),'w') # A vs B protein, no need to change
		fproteinup = open(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.up.txt' % (samplea, sampleb),'w') # protein, up
		fproteindn = open(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.dn.txt' % (samplea, sampleb),'w') # protein, up
		lines = f.readlines()
		f.close()
		headlist = lines[0].strip('\n').split('\t')
		print >>f2, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fprotein, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fproteinup, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>fproteindn, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == samplea and linelist[4] == sampleb:  # 样品顺序不一样，fpkm要反过来， logfc不变 # 貌似这里不用
				print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
				if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
					print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if float(linelist[9])>=0 or linelist[9]=='inf':  # >=0或者是inf的归为上调
						newfc = log2fc2fc(linelist[9])
						regulation = 'up'
						print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:									  # 其他的归为下调（即<0或者-inf）
						newfc = log2fc2fc(linelist[9])
						regulation = 'down'
						print >>fproteindn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
			elif linelist[4] == samplea and linelist[5] == sampleb: #样品顺序一样，则fpkm不变，logfc取反
				if linelist[9] == 'inf': # inf
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						#sampleAfpkm = float(linelist[7])+1.0
						#sampleBfpkm = float(linelist[8])+1.0
						newfc = log2fc2fc(linelist[9])
						regulation='down'
						print >>fproteindn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] =='-inf': # -inf
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						#sampleAfpkm = float(linelist[7])+1.0
						#sampleBfpkm = float(linelist[8])+1.0
						newfc = log2fc2fc(linelist[9])
						regulation='up'
						print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] == '0':
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						newfc='1.0'
						regulation='up'
						print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						
						
						if float(linelist[9])<0:
							newfc = log2fc2fc(linelist[9])
							regulation = 'up'
							print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						else:
							newfc = log2fc2fc(linelist[9])
							regulation = 'down'
							print >>fproteindn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
			else:
				continue
		f2.close()
		fprotein.close()
		fproteinup.close()
		fproteindn.close()
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.txt' % (samplea, sampleb), geneinfofile) # added gene info
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.up.txt' % (samplea, sampleb), geneinfofile) # added gene info
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.dn.txt' % (samplea, sampleb), geneinfofile) # added gene info


###################################### 20 two sample all comparison gene (txt), over

def at_least_one_expressed(valuelist):
	if len([i for i in valuelist if float(i) >= FACTOR]) >=1:
		return 1
	else:
		return 0

###################################### 21 two sample de gene (txt), start


def uniquegopw(file):
	"""
	remove '-' symbols in go, pathway, for rat, zebrafish .....
	"""
	f=open(file)
	lines = []
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[0]!= '-':
			lines.append(line.strip())
	f.close()
	f2=open('u.'+file,'w')
	for i in set(lines):
		print >>f2, i
	f2.close()
	os.remove(file)
	os.rename('u.'+file, file)

def two_sample_comparison_de_gene(): # ok with protein_coding results, and up, down go, pathway results
	"""
	test_id				gene_id				gene	locus					sample_1	sample_2	status	value_1		value_2	log2(fold_change)	test_stat	p_value	q_value	significant
	ENSMUSG00000000058	ENSMUSG00000000058	Cav2	chr6:17197750-17385604	6W-WT-FM-1	6W-KO-FM-2	OK		0.245057	14.4872	5.88552				1.05052		0.20925	0.366351	no
	sample_1 (6W-WT-FM-1) 
	sample_2 (6W-KO-FM-2)
	value_1  (0.245057)
	value_2  (14.4872)
	
	sample_2 vs sample_1 (6W-KO-FM-2 vs 6W-WT-FM-1)
	14.4872/0.245057 = 59.11767
	2**5.88552 = 59.11767
	
	The name order (A, B) was the same to value order (A, B), but the logfc is opposite
	so if you want to keep the logfc, the sample name order should be reversed,
	or just add minus in front of the logfc
	
	gene_exp.diff.6W-KO-M-2_vs_6W-KO-FM-2.txt
	"""
	
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0, 0.05
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	gene_biotype_strand_dict = parse_gtf_gene(hisat2_gtf)
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		abfc = float(comparison[2])
		abp = float(comparison[3])
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		f =open(foldname+'/gene_exp.diff')  # all comparisons both protein and lnc
		fup = open(foldname+'/gene_exp.diff.%s_vs_%s.up.txt' % (samplea, sampleb),'w') #     A vs B.up protein
		fdn = open(foldname+'/gene_exp.diff.%s_vs_%s.dn.txt' % (samplea, sampleb),'w') #     A vs B.dn protein
		# gene_exp.diff.33_vs_S0.up.protein_coding.txt
		fupprotein = open(foldname+'/gene_exp.diff.%s_vs_%s.up.protein_codingtmp.txt' % (samplea, sampleb),'w')
		fupprotein_go = open(foldname+'/go.mRNA.%s_vs_%s.up.txt' % (samplea, sampleb),'w')
		fupprotein_pw = open(foldname+'/pathway.mRNA.%s_vs_%s.up.txt' % (samplea, sampleb),'w')
		fdnprotein = open(foldname+'/gene_exp.diff.%s_vs_%s.dn.protein_codingtmp.txt' % (samplea, sampleb),'w')
		fdnprotein_go = open(foldname+'/go.mRNA.%s_vs_%s.dn.txt' % (samplea, sampleb),'w')
		fdnprotein_pw = open(foldname+'/pathway.mRNA.%s_vs_%s.dn.txt' % (samplea, sampleb),'w')
		
		lines = f.readlines()
		headlist = lines[0].strip('\n').split('\t')
		print >>fup, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fdn, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fupprotein, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>fdnprotein, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == samplea and linelist[4] == sampleb:
				if linelist[9] == 'inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
				elif linelist[9] == '-inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
					
				else:
					if float(linelist[9]) >= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
						
					elif float(linelist[9]) <= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] =='protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
						
					else:
						continue
			# newto
			elif linelist[4] == samplea and linelist[5] == sampleb:
				if linelist[9] == 'inf':
					if float(linelist[11]) <=abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
						
				elif linelist[9] =='-inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] =='protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
						
				else:
					if float(linelist[9])>= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
						
					elif float(linelist[9])<= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] =='protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
					else:
						continue
			else:
				continue
		f.close()
		fup.close()
		fdn.close()
		fupprotein.close()
		fdnprotein.close()
		fupprotein_go.close()
		fupprotein_pw.close()
		fdnprotein_go.close()
		fdnprotein_pw.close()
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.up.protein_codingtmp.txt' % (samplea, sampleb), geneinfofile) # added gene info
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.dn.protein_codingtmp.txt' % (samplea, sampleb), geneinfofile) # added gene info
		# 		fupprotein_go = open(foldname+'/go.mRNA.%s_vs_%s.up.txt' % (samplea, sampleb),'w')
		os.chdir(foldname)
		uniquegopw('pathway.mRNA.%s_vs_%s.up.txt' % (samplea, sampleb))
		uniquegopw('go.mRNA.%s_vs_%s.up.txt' % (samplea, sampleb))
		uniquegopw('go.mRNA.%s_vs_%s.dn.txt' % (samplea, sampleb))
		uniquegopw('pathway.mRNA.%s_vs_%s.dn.txt' % (samplea, sampleb))
		os.chdir(global_now_location)


###################################### 21 two sample de gene (txt), over, gene level


###################################### 22 two sample de gene (excel), start
def excel_two_sample_comparison_de_gene():
	"""
	gene_exp.diff.6W-KO-M-2_vs_6W-KO-FM-2.txt
	"""
	# gene	gene_id	biotype	strand	locus	6W-KO-M-2_FPKM	6W-KO-FM-2_FPKM	log2(fold change)	p_value	q_value
	two_sample_comparison_de_gene()
	############# for all
	excel_title_info_two_sample_de_gene = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: log2(fold change), fold change between two samples, differentially genes were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two samples.
Column J: FDR, adjusted p-value."""

	# write into excel 
	workbook = Workbook('Differentially Expressed Genes (protein+lnc).xlsx')  # xlsx's name
	annotationformat =  workbook.add_format({'align': 'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	uptitleformat = workbook.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat = workbook.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0, 0.05
	for comparison in comparisons:
		#print comparison
		#print type(comparison[-1])
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		abfc = float(comparison[2])
		abp = float(comparison[3])
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		#print samplea, sampleb, abfc, abp
		# up data
		updata = read_data_content_with_float(foldname+'/gene_exp.diff.%s_vs_%s.up.txt' % (samplea, sampleb))
		#print updata[0]
		# gene_exp.diff.33_vs_S0.dn.protein_coding.txt
		upsheetname = 'up.%s_vs_%s' % (samplea, sampleb)
		upworksheet = workbook.add_worksheet(upsheetname[:30])  # upsheetname
		upcolnum = len(updata[0])
		upworksheet.merge_range(0, 0, 19+info_addline, upcolnum-1, '', annotationformat)
		upworksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb,'\nFold Change cutoff: ', red, comparison[2], '\np-value cutoff: ', red, comparison[3],  '\nFPKM >= ', red, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_gene, annotationformat)
		updata_title = updata[0]
		uprownum = len(updata)
		upworksheet.merge_range(21+info_addline, 0, 21+info_addline, upcolnum-1, '%s vs %s up-regulated genes' % (samplea, sampleb), uptitleformat)
		for i in range(upcolnum):
			upworksheet.write(0+22+info_addline,i,updata_title[i], titleformat)
		for row in range(1, uprownum):
			for col in range(upcolnum):
				upworksheet.write(row+22+info_addline, col, updata[row][col],zhengwenformat)
		# down data
		dndata = read_data_content_with_float(foldname+'/gene_exp.diff.%s_vs_%s.dn.txt' % (samplea, sampleb))
		# gene_exp.diff.33_vs_S0.dn.txt
		dnsheetname = 'dn.%s_vs_%s' % (samplea, sampleb)
		dnworksheet = workbook.add_worksheet(dnsheetname[:30])  # dnsheetname
		dncolnum = len(dndata[0])
		dnworksheet.merge_range(0, 0, 19+info_addline, dncolnum-1, '', annotationformat)
		dnworksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb,'\nFold Change cutoff: ', red, comparison[2], '\np-value cutoff: ', red, comparison[3],  '\nFPKM >= ', red, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_gene, annotationformat)
		dndata_title = dndata[0]
		dnrownum = len(dndata)
		dnworksheet.merge_range(21+info_addline, 0, 21+info_addline, dncolnum-1, '%s vs %s down-regulated genes' % (samplea, sampleb), dntitleformat)
		for ii in range(dncolnum):
			dnworksheet.write(0+22+info_addline,ii,dndata_title[ii], titleformat)
		for row2 in range(1, dnrownum):
			for col2 in range(dncolnum):
				dnworksheet.write(row2+22+info_addline, col2, dndata[row2][col2],zhengwenformat)
	workbook.close()
 
		############# for protein coding
	excel_title_info_two_sample_de_gene_protein_coding = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, protein_coding.
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: fold change, fold change between two samples (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two samples (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two samples.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ V: annotation information, including GeneID, Synonyms, dbXrefs (other database ID), chromosome, map_location, description, pathway and GO annotation.

note: since the different sources, some genes without annotation information."""

	# write into excel 
	workbook2 = Workbook('Differentially Expressed mRNAs.xlsx')  # xlsx's name
	annotationformat2 =  workbook2.add_format({'align': 'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	uptitleformat2 = workbook2.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat2 = workbook2.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0, 0.05
	for comparison2 in comparisons2:
		samplea2 = comparison2[0].replace('_sequence.fastq', '')
		sampleb2 = comparison2[1].replace('_sequence.fastq', '')
		abfc2 = float(comparison2[2])
		abp2 = float(comparison2[3])
		foldname2 = 'cuffdiff_'+ samplea2.split('_sequence.fastq')[0]+'_vs_'+sampleb2.split('_sequence.fastq')[0]
		# up2 data
		updata2 = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.up.protein_coding.txt' % (samplea2, sampleb2))
		upsheetname2 = 'up.%s_vs_%s' % (samplea2, sampleb2)
		upworksheet2 = workbook2.add_worksheet(upsheetname2[:30])  # upsheetname2
		upcolnum2 = len(updata2[0])
		upworksheet2.merge_range(0, 0, 21+info_addline, upcolnum2-1, '', annotationformat2)
		upworksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.up','\nFold Change cutoff: ', red2, comparison2[2], '\np-value cutoff: ', red2, comparison2[3],  '\nFPKM >= ', red2, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_gene_protein_coding, annotationformat2)
		updata_title2 = updata2[0]
		uprownum2 = len(updata2)
		upworksheet2.merge_range(23+info_addline, 0, 23+info_addline, upcolnum2-1, '%s vs %s up-regulated genes' % (samplea2, sampleb2), uptitleformat2)
		for i2 in range(upcolnum2):
			upworksheet2.write(0+24+info_addline,i2,updata_title2[i2].replace('q_value','FDR'), titleformat2)
		for row2 in range(1, uprownum2):
			for col2 in range(upcolnum2):
				upworksheet2.write(row2+24+info_addline, col2, updata2[row2][col2],zhengwenformat2)
		# down2 data
		dndata2 = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.dn.protein_coding.txt' % (samplea2, sampleb2))
		dnsheetname2 = 'dn.%s_vs_%s' % (samplea2, sampleb2)
		dnworksheet2 = workbook2.add_worksheet(dnsheetname2[:30])  # dnsheetname2
		dncolnum2 = len(dndata2[0])
		dnworksheet2.merge_range(0, 0, 21+info_addline, dncolnum2-1, '', annotationformat2)
		dnworksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.down','\nFold Change cutoff: ', red2, comparison2[2], '\np-value cutoff: ', red2, comparison2[3],  '\nFPKM >= ', red2, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_gene_protein_coding, annotationformat2)
		dndata_title2 = dndata2[0]
		dnrownum2 = len(dndata2)
		dnworksheet2.merge_range(23+info_addline, 0, 23+info_addline, dncolnum2-1, '%s vs %s down-regulated genes' % (samplea2, sampleb2), dntitleformat2)
		for ii2 in range(dncolnum2):
			dnworksheet2.write(0+24+info_addline,ii2,dndata_title2[ii2].replace('q_value','FDR'), titleformat2)
		for row22 in range(1, dnrownum2):
			for col22 in range(dncolnum2):
				dnworksheet2.write(row22+24+info_addline, col22, dndata2[row22][col22],zhengwenformat2)
	workbook2.close()
  

###################################### 22 two sample de gene (excel), over
#excel_two_sample_comparison_de_gene()

 
###################################### 23 two sample de transcript (txt), start

def two_sample_comparison_de_transcript(): # 1, all 2, lncRNA transcript level (add associated)
	"""
test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
ENSMUST00000000001	ENSMUSG00000000001	Gnai3	chr3:108107279-108146146	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0.175945	inf	0	1	1	no
ENSMUST00000000003	ENSMUSG00000000003	Pbsn	chrX:77837900-77853623	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0	0	0	1	1	no
	
	The name order (A, B) was the same to value order (B, A), but the logfc is opposite
	so if you want to keep logfc, the sample name order should be reversed,
	or just add minus in front of the logfc
	"""
	
	comparisons = get_2_samples_comparisons()  # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	tophat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	transcript_biotype_strand_dict = parse_gtf_transcript(tophat2_gtf)
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		abfc = float(comparison[2])
		abp = float(comparison[3])
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		f =open(foldname+'/isoform_exp.diff')
		fup = open(foldname+'/isoform_exp.diff.%s_vs_%s.up.txt' % (samplea, sampleb),'w')
		fdn = open(foldname+'/isoform_exp.diff.%s_vs_%s.dn.txt' % (samplea, sampleb),'w')
		flncup = open(foldname+'/isoform_exp.diff.%s_vs_%s.up.lnctmp.txt' % (samplea, sampleb),'w')
		flncdn = open(foldname+'/isoform_exp.diff.%s_vs_%s.dn.lnctmp.txt' % (samplea, sampleb),'w')
		lines = f.readlines()
 
		headlist = lines[0].strip('\n').split('\t')
		print >>fup, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fdn, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>flncup, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>flncdn, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == samplea and linelist[4] == sampleb:
				if linelist[9] =='inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] =='-inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					if float(linelist[9]) >= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]) :
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					elif float(linelist[9]) <= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:
						continue
			elif linelist[4] == samplea and linelist[5] == sampleb:
				if linelist[9] == 'inf':
					if float(linelist[11]) <=abp and at_least_one_expressed([linelist[8], linelist[7]]) :
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] =='-inf':
					if float(linelist[11]) <=abp and at_least_one_expressed([linelist[8], linelist[7]]) :
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					if float(linelist[9])>= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					elif float(linelist[9])<= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:
						continue
						
			else:
				continue
		f.close()
		fup.close()
		fdn.close()
		flncup.close()
		flncdn.close()
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.up.lnctmp.txt' % (samplea, sampleb), lnc_associated_genefile) # added, with header
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.dn.lnctmp.txt' % (samplea, sampleb), lnc_associated_genefile) # added, with header
		

###################################### 23 two sample de transcript (txt), over

###################################### 24 two sample de transcript (excel), start

#two_sample_comparison_de_transcript()
def excel_two_sample_comparison_de_transcript(): # 1, all, 2, lncRNA transcript level (with associated gene)
	"""
	gene_exp.diff.6W-KO-M-2_vs_6W-KO-FM-2.txt
	"""
	# gene	gene_id	biotype	strand	locus	6W-KO-M-2_FPKM	6W-KO-FM-2_FPKM	log2(fold change)	p_value	q_value
	two_sample_comparison_de_transcript()
	excel_title_info_two_sample_de_transcript = """

Column A: transcript_id, the transcript id.
Column B: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: log2(fold change), fold change between two samples, differentially expressed lncRNAs were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two samples.
Column J: FDR, adjusted p-value."""

	# write into excel 
	workbook = Workbook('Differentially Expressed (protein+lnc).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True})
	uptitleformat = workbook.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat = workbook.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0, 0.05
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		abfc = float(comparison[2])
		abp = float(comparison[3])
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		# up data
		updata = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.up.txt' % (samplea, sampleb))
		upsheetname = 'up.%s_vs_%s' % (samplea, sampleb)
		upworksheet = workbook.add_worksheet(upsheetname[:30])  # upsheetname
		upcolnum = len(updata[0])
		upworksheet.merge_range(0, 0, 19+info_addline, upcolnum-1, '', annotationformat)
		upworksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb,'\nFold Change cutoff: ', red, comparison[2], '\n','p-value cutoff: ', red, comparison[3], '\nFPKM >= ', red, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_transcript, annotationformat)
		updata_title = updata[0]
		uprownum = len(updata)
		upworksheet.merge_range(21+info_addline, 0, 21+info_addline, upcolnum-1, '%s vs %s up-regulated lncRNAs' % (samplea, sampleb), uptitleformat)
		for i in range(upcolnum):
			upworksheet.write(0+22+info_addline,i,updata_title[i], titleformat)
		for row in range(1, uprownum):
			for col in range(upcolnum):
				upworksheet.write(row+22+info_addline, col, updata[row][col],zhengwenformat)
		# down data
		dndata = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.dn.txt' % (samplea, sampleb))
		dnsheetname = 'dn.%s_vs_%s' % (samplea, sampleb)
		dnworksheet = workbook.add_worksheet(dnsheetname[:30])  # dnsheetname
		dncolnum = len(dndata[0])
		dnworksheet.merge_range(0, 0, 19+info_addline, dncolnum-1, '', annotationformat)
		dnworksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb,'\nFold Change cutoff: ', red, comparison[2], '\n','p-value cutoff: ', red, comparison[3],  '\nFPKM >= ', red, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_transcript, annotationformat)
		dndata_title = dndata[0]
		dnrownum = len(dndata)
		dnworksheet.merge_range(21+info_addline, 0, 21+info_addline, dncolnum-1, '%s vs %s down-regulated lncRNAs' % (samplea, sampleb), dntitleformat)
		for ii in range(dncolnum):
			dnworksheet.write(0+22+info_addline,ii,dndata_title[ii], titleformat)
		for row2 in range(1, dnrownum):
			for col2 in range(dncolnum):
				dnworksheet.write(row2+22+info_addline, col2, dndata[row2][col2],zhengwenformat)
	workbook.close()
	
	# for lncRNA 2 samples DE with associated
	
	excel_title_info_two_sample_de_lnc_transcript = """

Column A: transcript_id, the transcript id.
Column B: biotype, type of gene, including antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: fold change, fold change between two samples (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two samples (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two samples.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ U: the genomic organization of LncRNAs and the information of the associated coding genes.

the classes of LncRNAs defined by genomic organization:
"exon sense-overlapping": the LncRNA's exon is overlapping a coding transcript exon on the same genomic strand;
"intron sense-overlapping": the LncRNA is overlapping the intron of a coding transcript on the same genomic strand;
"intronic antisense": the LncRNA is overlapping the intron of a coding transcript on the antisense strand;
"natural antisense": the LncRNA is transcribed from the antisense strand and overlapping with a coding transcript; 
"bidirectional": the LncRNA is oriented head to head to a coding transcript within 1000 bp; 
"intergenic": there are no overlapping or bidirectional coding transcripts nearby the LncRNA."""

	# write into excel
	workbook2 = Workbook('Differentially Expressed LncRNAs.xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 =  workbook2.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True})
	uptitleformat2 = workbook2.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat2 = workbook2.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0, 0.05
	for comparison2 in comparisons2:
		samplea2 = comparison2[0].replace('_sequence.fastq', '')
		sampleb2 = comparison2[1].replace('_sequence.fastq', '')
		abfc2 = float(comparison2[2])
		abp2 = float(comparison2[3])
		foldname2 = 'cuffdiff_'+ samplea2.split('_sequence.fastq')[0]+'_vs_'+sampleb2.split('_sequence.fastq')[0]
		# up lncRNA		
		updata2 = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.up.lnc.txt' % (samplea2, sampleb2))
		upsheetname2 = 'up.%s_vs_%s' % (samplea2, sampleb2)
		upworksheet2 = workbook2.add_worksheet(upsheetname2[:30])  # upsheetname
		upcolnum2 = len(updata2[0])
		upworksheet2.merge_range(0, 0, 19+info_addline, upcolnum2-1, '', annotationformat2)
		upworksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.up','\nFold Change cutoff: ', red2, comparison2[2], '\n','p-value cutoff: ', red2, comparison2[3], '\nFPKM >= ', red2, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_lnc_transcript, annotationformat2)
		updata_title2 = updata2[0]
		uprownum2 = len(updata2)
		upworksheet2.merge_range(21+info_addline, 0, 21+info_addline, upcolnum2-1, '%s vs %s up-regulated lncRNAs' % (samplea2, sampleb2), uptitleformat2)
		for i2 in range(upcolnum2):
			upworksheet2.write(0+22+info_addline,i2,updata_title2[i2].replace('q_value','FDR'), titleformat2)
		for row2 in range(1, uprownum2):
			for col2 in range(upcolnum2):
				upworksheet2.write(row2+22+info_addline, col2, updata2[row2][col2],zhengwenformat2)
		# down lncRNA
		dndata2 = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.dn.lnc.txt' % (samplea2, sampleb2))
		dnsheetname2 = 'dn.%s_vs_%s' % (samplea2, sampleb2)
		dnworksheet2 = workbook2.add_worksheet(dnsheetname2[:30])  # dnsheetname
		dncolnum2 = len(dndata2[0])
		dnworksheet2.merge_range(0, 0, 19+info_addline, dncolnum2-1, '', annotationformat2)
		dnworksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.down','\nFold Change cutoff: ', red2, comparison2[2], '\n','p-value cutoff: ', red2, comparison2[3],  '\nFPKM >= ', red2, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_lnc_transcript, annotationformat2)
		dndata_title2 = dndata2[0]
		dnrownum2 = len(dndata2)
		dnworksheet2.merge_range(21+info_addline, 0, 21+info_addline, dncolnum2-1, '%s vs %s down-regulated lncRNAs' % (samplea2, sampleb2), dntitleformat2)
		for ii2 in range(dncolnum2):
			dnworksheet2.write(0+22+info_addline,ii2,dndata_title2[ii2].replace('q_value','FDR'), titleformat2)
		for row22 in range(1, dnrownum2):
			for col22 in range(dncolnum2):
				dnworksheet2.write(row22+22+info_addline, col22, dndata2[row22][col22],zhengwenformat2)
	workbook2.close()

###################################### 24 two sample de transcript (excel), over
#excel_two_sample_comparison_de_transcript()

###################################### 25 two sample All comparison gene (excel), start

def excel_two_sample_comparison_all_gene():
	"""
	gene_exp.diff.6W-KO-M-2_vs_6W-KO-FM-2.txt
	"""
	# gene	gene_id	biotype	strand	locus	6W-KO-M-2_FPKM	6W-KO-FM-2_FPKM	log2(fold change)	p_value	q_value
	excel_title_info_two_sample_all_comparison = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: log2(fold change), fold change between two samples, differentially expressed genes were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two samples.
Column J: FDR, adjusted p-value."""

	# write into excel 
	two_sample_comparison_gene()
	workbook = Workbook('All Comparison (2 sample protein+lnc).xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		#abfc = float(comparison[1].replace('_sequence.fastq'))
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		data = read_data_content_with_float(foldname+'/gene_exp.diff.%s_vs_%s.txt' % (samplea, sampleb))
		sheetname = '%s_vs_%s' % (samplea, sampleb)
		worksheet = workbook.add_worksheet(sheetname[:30])  # sheetname
		colnum = len(data[0])
		worksheet.merge_range(0, 0, 19+info_addline, colnum-1, '', annotationformat)
		worksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb, excel_title_info_two_sample_all_comparison, annotationformat)
		data_title = data[0]
		rownum = len(data)
		for i in range(colnum):
			worksheet.write(0+22+info_addline,i,data_title[i].replace('q_value','FDR'), titleformat)
			                                                          #q_value
		for row in range(1, rownum):
			for col in range(colnum):
				worksheet.write(row+22+info_addline, col, data[row][col])
	workbook.close()
	
	############## for protein coding
	
	excel_title_info_two_sample_all_comparison_protein_coding = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, protein_coding.
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: fold change, fold change between two samples (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two samples (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two samples.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ V: annotation information, including GeneID, Synonyms, dbXrefs (other database ID), chromosome, map_location, description, pathway and GO annotation.

note: since the different sources, some genes without annotation information"""

	workbook2 = Workbook('All Comparison (2samples, mRNA).xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 =  workbook2.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	for comparison2 in comparisons2:
		samplea2 = comparison2[0].replace('_sequence.fastq', '')
		sampleb2 = comparison2[1].replace('_sequence.fastq', '')
		#abfc = float(comparison[1].replace('_sequence.fastq')) # since do not filter
		foldname2 = 'cuffdiff_'+ samplea2.split('_sequence.fastq')[0]+'_vs_'+sampleb2.split('_sequence.fastq')[0]
		
		# up, testing
		data2up = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.protein_coding.up.txt' % (samplea2, sampleb2))
		sheetname2up = 'up.%s_vs_%s' % (samplea2, sampleb2)
		worksheet2 = workbook2.add_worksheet(sheetname2up[:30])  # sheetname2
		colnum2up = len(data2up[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2up-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.up', excel_title_info_two_sample_all_comparison_protein_coding, annotationformat2)
		data_title2up = data2up[0]
		rownum2up = len(data2up)
		print data_title2up
		for i2up in range(colnum2up):
			worksheet2.write(0+22+info_addline, i2up, data_title2up[i2up].replace('q_value','FDR'), titleformat2)
		for row2up in range(1, rownum2up):
			for col2up in range(colnum2up):
				worksheet2.write(row2up+22+info_addline, col2up, data2up[row2up][col2up],zhengwenformat2)   
		
		# dn, testing
		data2dn = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.protein_coding.dn.txt' % (samplea2, sampleb2))
		sheetname2dn = 'dn.%s_vs_%s' % (samplea2, sampleb2)
		worksheet2 = workbook2.add_worksheet(sheetname2dn[:30])  # sheetname2
		colnum2dn = len(data2dn[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2dn-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.down', excel_title_info_two_sample_all_comparison_protein_coding, annotationformat2)
		data_title2dn = data2dn[0]
		rownum2dn = len(data2dn)
		for i2dn in range(colnum2dn):
			worksheet2.write(0+22+info_addline, i2dn, data_title2dn[i2dn].replace('q_value','FDR'), titleformat2)
		for row2dn in range(1, rownum2dn):
			for col2dn in range(colnum2dn):
				worksheet2.write(row2dn+22+info_addline, col2dn, data2dn[row2dn][col2dn],zhengwenformat2)
	workbook2.close()

###################################### 25 two sample All comparison gene (excel), over

###################################### 26 DE of transcript expression (txt), start


def two_sample_comparison_transcript(): # 1, all, 2, lncRNA, associated
	"""
test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
ENSMUST00000000001	ENSMUSG00000000001	Gnai3	chr3:108107279-108146146	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0.175945	inf	0	1	1	no
ENSMUST00000000003	ENSMUSG00000000003	Pbsn	chrX:77837900-77853623	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0	0	0	1	1	no
	
	The name order (A, B) was the same to value order (B, A), but the logfc is opposite
	so if you want to keep logfc, the sample name order should be reversed,
	or just add minus in front of the logfc
	"""
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	transcript_biotype_strand_dict = parse_gtf_transcript(hisat2_gtf)
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		#abfc = float(comparison[1].replace('_sequence.fastq')), since no filtering
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		f =open(foldname+'/isoform_exp.diff')
		f2 = open(foldname+'/isoform_exp.diff.%s_vs_%s.txt' % (samplea, sampleb),'w')
		flnc2 = open(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.txt' % (samplea, sampleb),'w')
		flnc2up = open(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.up.txt' % (samplea, sampleb),'w')
		flnc2dn = open(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.dn.txt' % (samplea, sampleb),'w')
		lines = f.readlines()

		headlist = lines[0].strip('\n').split('\t')
		print >>f2, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>flnc2, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>flnc2up, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>flnc2dn, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+samplea+'_FPKM'+'\t'+sampleb+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == samplea and linelist[4] == sampleb: # 样品顺序不一样，fpkm要反过来， logfc不变 # 貌似这里不用
				print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
				if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
					print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if float(linelist[9])>=0 or linelist[9]=='inf':  # >=0或者是inf的归为上调
						newfc = log2fc2fc(linelist[9])
						regulation = 'up'
						print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:									  # 其他的归为下调（即<0或者-inf）
						newfc = log2fc2fc(linelist[9])
						regulation = 'down'
						print >>flnc2dn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
			elif linelist[4] == samplea and linelist[5] == sampleb: #样品顺序一样，则fpkm不变，logfc取反
				if linelist[9] == 'inf':
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='down'
						print >>flnc2dn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] == '-inf':
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='up'
						print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] == '0':
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						newfc='1.0'
						regulation='up'
						print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if float(linelist[9])<0:
							newfc = log2fc2fc(linelist[9])
							regulation = 'up'
							print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						else:
							newfc = log2fc2fc(linelist[9])
							regulation = 'down'
							print >>flnc2dn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
			else:
				continue
		f.close()
		f2.close()
		flnc2.close()
		flnc2up.close()
		flnc2dn.close()
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.txt' % (samplea, sampleb), lnc_associated_genefile)
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.up.txt' % (samplea, sampleb), lnc_associated_genefile)
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.dn.txt' % (samplea, sampleb), lnc_associated_genefile)
		
###################################### 26 DE of transcript expression txt, over (all , lnc)

###################################### 27 All comparison of transcript expression excel, start
def excel_two_sample_comparison_all_transcript():
	"""
	isoform_exp.diff.6W-KO-M-2_vs_6W-WT-FM-1.txt
	"""
	# test_id	biotype	strand	gene_id	locus	6W-KO-M-2_FPKM	6W-WT-FM-1_FPKM	log2(fold change)	p_value	q_value
	excel_title_info_two_sample_all_comparison = """

Column A: transcript_id, the Ensembl transcript id.
Column B: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: log2(fold change), fold change between two samples, differentially expressed lncRNAs were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two samples.
Column J: FDR, adjusted p-value."""

	# write into excel 
	two_sample_comparison_transcript()
	workbook = Workbook('All Comparison (2 sample, protein+lnc).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat = workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2_samples_comparisons()  # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		#abfc = float(comparison[1].replace('_sequence.fastq'))
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		data = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.txt' % (samplea, sampleb))
		sheetname = '%s_vs_%s' % (samplea, sampleb)
		worksheet = workbook.add_worksheet(sheetname[:30])  # sheetname
		colnum = len(data[0])
		worksheet.merge_range(0, 0, 19+info_addline, colnum-1, '', annotationformat)
		worksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb, excel_title_info_two_sample_all_comparison, annotationformat)
		data_title = data[0]
		rownum = len(data)
		for i in range(colnum):
			worksheet.write(0+22+info_addline,i,data_title[i].replace('q_value', 'FDR'), titleformat)
		for row in range(1, rownum):
			for col in range(colnum):
				worksheet.write(row+22+info_addline, col, data[row][col], zhengwenformat)
	workbook.close()
	
	# for lncRNA all comparison (2 samples)
	
	excel_title_info_two_sample_all_comparison_lnc = """

Column A: transcript_id, the Ensembl transcript id.
Column B: biotype, type of gene, including antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: fold change, fold change between two samples (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two samples (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two samples.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ U: the genomic organization of LncRNAs and the information of the associated coding genes.

the classes of LncRNAs defined by genomic organization:
"exon sense-overlapping": the LncRNA's exon is overlapping a coding transcript exon on the same genomic strand;
"intron sense-overlapping": the LncRNA is overlapping the intron of a coding transcript on the same genomic strand;
"intronic antisense": the LncRNA is overlapping the intron of a coding transcript on the antisense strand;
"natural antisense": the LncRNA is transcribed from the antisense strand and overlapping with a coding transcript; 
"bidirectional": the LncRNA is oriented head to head to a coding transcript within 1000 bp; 
"intergenic": there are no overlapping or bidirectional coding transcripts nearby the LncRNA."""

	# write into excel
	workbook2 = Workbook('All Comparison (2 sample, LncRNA).xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 = workbook2.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2_samples_comparisons()  # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	for comparison2 in comparisons2:
		samplea2 = comparison2[0].replace('_sequence.fastq', '')
		sampleb2 = comparison2[1].replace('_sequence.fastq', '')
		#abfc = float(comparison[1].replace('_sequence.fastq'))
		foldname2 = 'cuffdiff_'+ samplea2.split('_sequence.fastq')[0]+'_vs_'+sampleb2.split('_sequence.fastq')[0]
		
		# up, testing
		data2up = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.lnc.up.txt' % (samplea2, sampleb2))
		sheetname2up = 'up.%s_vs_%s' % (samplea2, sampleb2)
		worksheet2 = workbook2.add_worksheet(sheetname2up[:30])  # sheetname
		colnum2up = len(data2up[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2up-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.up', excel_title_info_two_sample_all_comparison_lnc, annotationformat2)
		data_title2up = data2up[0]
		rownum2up = len(data2up)
		for i2up in range(colnum2up):
			worksheet2.write(0+22+info_addline, i2up, data_title2up[i2up].replace('q_value', 'FDR'), titleformat2)
		for row2up in range(1, rownum2up):
			for col2up in range(colnum2up):
				worksheet2.write(row2up+22+info_addline, col2up, data2up[row2up][col2up],zhengwenformat2)
		
		# down, testing
		data2dn = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.lnc.dn.txt' % (samplea2, sampleb2))
		sheetname2dn = 'dn.%s_vs_%s' % (samplea2, sampleb2)
		worksheet2 = workbook2.add_worksheet(sheetname2dn[:30])  # sheetname
		colnum2dn = len(data2dn[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2dn-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, samplea2+' vs '+sampleb2+'.down', excel_title_info_two_sample_all_comparison_lnc, annotationformat2)
		data_title2dn = data2dn[0]
		rownum2dn = len(data2dn)
		for i2dn in range(colnum2dn):
			worksheet2.write(0+22+info_addline, i2dn, data_title2dn[i2dn].replace('q_value', 'FDR'), titleformat2)
		for row2dn in range(1, rownum2dn):
			for col2dn in range(colnum2dn):
				worksheet2.write(row2dn+22+info_addline, col2dn, data2dn[row2dn][col2dn],zhengwenformat2)
	workbook2.close()
###################################### 27 All comparison of transcript expression excel, over (all, lncRNA)


###################################### 28 group vs group all gene (txt), start
def two_group_comparison_gene():
	"""
	similar to two sample comparison gene
	gene_exp.diff.G2_vs_G1.txt into up and down
	"""
	comparisons = get_2group_comparisons()
	#print comparisons
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	gene_biotype_strand_dict = parse_gtf_gene(hisat2_gtf)
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		f =open(foldname+'/gene_exp.diff') # all comparisons, both mRNA and LncRNA
		f2 = open(foldname+'/gene_exp.diff.%s_vs_%s.txt' % (groupAname, groupBname),'w') # A vs B all
		fprotein = open(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.txt' % (groupAname, groupBname),'w') # A vs B, protein
		fproteinup = open(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.up.txt' % (groupAname, groupBname),'w') # A vs B, protein, up
		fproteindn = open(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.dn.txt' % (groupAname, groupBname),'w') # A vs B, protein, dn
		lines = f.readlines()
		f.close()
		headlist = lines[0].strip('\n').split('\t')
		print headlist
		print >>f2, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fprotein, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fproteinup, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>fproteindn, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == groupAname and linelist[4] == groupBname: # 样品顺序不一样，fpkm要反过来， logfc不变 # 貌似这里不用
				print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
				if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
					print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if float(linelist[9])>=0 or linelist[9]=='inf': # >=0或者是inf的归为上调
						newfc = log2fc2fc(linelist[9])
						regulation = 'up'
						print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:									  # 其他的归为下调（即<0或者-inf）
						newfc = log2fc2fc(linelist[9])
						regulation = 'down'
						print >>fproteindn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation

			elif linelist[4] == groupAname and linelist[5] == groupBname:  #样品顺序一样，则fpkm不变，logfc取反
				if linelist[9] == 'inf':
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='down'
						print >>fproteindn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						
				elif linelist[9] == '-inf':
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='up'
						print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				
				elif linelist[9]=='0':
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = '1.0'
						regulation='up'
						print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation

				else:
					print >>f2, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
					if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
						print >>fprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if float(linelist[9])<0:
							newfc = log2fc2fc(linelist[9])
							regulation = 'up'
							print >>fproteinup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						else:
							newfc = log2fc2fc(linelist[9])
							regulation = 'down'
							print >>fproteindn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
			else:
				continue
		f2.close()  # all
		fprotein.close() # protein
		fproteinup.close() # protein up
		fproteindn.close() # protein dn
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.txt' % (groupAname, groupBname), geneinfofile) # added, geneinfo
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.up.txt' % (groupAname, groupBname), geneinfofile) # added, geneinfo
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.dn.txt' % (groupAname, groupBname), geneinfofile) # added, geneinfo
###################################### 26 group vs group all gene (txt, all, protein), over

###################################### 27 group vs group de gene (txt), start

def two_group_comparison_de_gene(): # with protein_coding results, and up, down go, pathway results
	"""
	similar to two sample comparison, without the FPKM of each samples, or means
	gene_exp.diff.G2_vs_G1.txt
	"""
	comparisons = get_2group_comparisons()
	print comparisons
	tophat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	gene_biotype_strand_dict = parse_gtf_gene(tophat2_gtf)
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		abfc = float(fc)
		abp = float(p)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		f =open(foldname+'/gene_exp.diff')
		fup = open(foldname+'/gene_exp.diff.%s_vs_%s.up.txt' % (groupAname, groupBname),'w') # A vs B
		fdn = open(foldname+'/gene_exp.diff.%s_vs_%s.dn.txt' % (groupAname, groupBname),'w') # A vs B
		# gene_exp.diff.33_vs_S0.up.protein_coding.txt # group vs group
		fupprotein = open(foldname+'/gene_exp.diff.%s_vs_%s.up.protein_codingtmp.txt' % (groupAname, groupBname),'w')   # protein
		fupprotein_go = open(foldname+'/go.mRNA.%s_vs_%s.up.txt' % (groupAname, groupBname),'w')
		fupprotein_pw = open(foldname+'/pathway.mRNA.%s_vs_%s.up.txt' % (groupAname, groupBname),'w')
		fdnprotein = open(foldname+'/gene_exp.diff.%s_vs_%s.dn.protein_codingtmp.txt' % (groupAname, groupBname),'w')
		fdnprotein_go = open(foldname+'/go.mRNA.%s_vs_%s.dn.txt' % (groupAname, groupBname),'w')
		fdnprotein_pw = open(foldname+'/pathway.mRNA.%s_vs_%s.dn.txt' % (groupAname, groupBname),'w')
		lines = f.readlines()
 
		headlist = lines[0].strip('\n').split('\t')
		print >>fup, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fdn, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fupprotein, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>fdnprotein, headlist[2]+'\t'+headlist[1]+'\tbiotype\tstrand\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == groupAname and linelist[4] == groupBname:
				if linelist[9] == 'inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
				elif linelist[9] == '-inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
				else:
					if float(linelist[9]) >= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
					elif float(linelist[9]) <= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
					else:
						continue
			elif linelist[4] == groupAname and linelist[5] == groupBname:
				if linelist[9] == 'inf':
					if float(linelist[11]) <=abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
				elif linelist[9] =='-inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
				else:
					if float(linelist[9])>= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc =log2fc2fc(linelist[9])
							regulation='down'
							print >>fdnprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fdnprotein_go, linelist[2]
							print >>fdnprotein_pw, linelist[2]+'\t'+pathway_colors['down']
					elif float(linelist[9])<= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if gene_biotype_strand_dict[linelist[1]].split()[0] == 'protein_coding':
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>fupprotein, linelist[2]+'\t'+linelist[1]+'\t'+gene_biotype_strand_dict[linelist[1]]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
							print >>fupprotein_go, linelist[2]
							print >>fupprotein_pw, linelist[2]+'\t'+pathway_colors['up']
					else:
						continue
			else:
				continue
		f.close()
		fup.close()
		fdn.close()
		fupprotein.close()
		fdnprotein.close()
		fupprotein_go.close()
		fdnprotein_go.close()
		fupprotein_pw.close()
		fdnprotein_pw.close()
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.up.protein_codingtmp.txt' % (groupAname, groupBname), geneinfofile) # added, geneinfo
		addgeneinfo(foldname+'/gene_exp.diff.%s_vs_%s.dn.protein_codingtmp.txt' % (groupAname, groupBname), geneinfofile) # added, geneinfo
		
		os.chdir(foldname)
		uniquegopw('go.mRNA.%s_vs_%s.up.txt' % (groupAname, groupBname))
		uniquegopw('pathway.mRNA.%s_vs_%s.up.txt' % (groupAname, groupBname))
		uniquegopw('go.mRNA.%s_vs_%s.dn.txt' % (groupAname, groupBname))
		uniquegopw('pathway.mRNA.%s_vs_%s.dn.txt' % (groupAname, groupBname))
		os.chdir(global_now_location)

###################################### 27 group vs group de gene (all, protein,txt), over


###################################### 28 two group de gene (excel), start
def excel_two_group_comparison_de_gene():
	"""
	gene_exp.diff.B_vs_A.txt, up and down
	"""
	# gene	gene_id	biotype	strand	locus	6W-KO-M-2_FPKM	6W-KO-FM-2_FPKM	log2(fold change)	p_value	q_value
	two_group_comparison_de_gene() # all and protein_coding
	excel_title_info_two_group_de_gene = """
Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: log2(fold change), fold change between two samples, differentially expressed genes were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two groups.
Column J: FDR, adjusted p-value."""

	# write into excel 
	workbook = Workbook('Differentially Expressed (lnc+protein,group).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	uptitleformat = workbook.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat = workbook.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2group_comparisons()
	#print comparisons
	#tophat2_gtf = cf.get('mapping_parameters','tophat2_gtf')
	#gene_biotype_strand_dict = parse_gtf_gene(tophat2_gtf)
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		abfc = float(fc)
		abp = float(p)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		# up data
		updata = read_data_content_with_float(foldname+'/gene_exp.diff.%s_vs_%s.up.txt' % (groupAname, groupBname))
		upsheetname = 'up.%s_vs_%s' % (groupAname, groupBname)
		upworksheet = workbook.add_worksheet(upsheetname[:30])  # upsheetname
		upcolnum = len(updata[0])
		upworksheet.merge_range(0, 0, 19+info_addline, upcolnum-1, '', annotationformat)
		upworksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname+'.up','\nFold Change: ', red, fc, '\n','p-value: ', red, p, '\nFPKM >= ', red, str(FACTOR), ' in at least one group', excel_title_info_two_group_de_gene, annotationformat)
		updata_title = updata[0]
		uprownum = len(updata)
		upworksheet.merge_range(21+info_addline, 0, 21+info_addline, upcolnum-1, '%s vs %s up-regulated lncRNA' % (groupAname, groupBname), uptitleformat)
		for i in range(upcolnum):
			upworksheet.write(0+22+info_addline,i,updata_title[i].replace('q_value','FDR'), titleformat)
		for row in range(1, uprownum):
			for col in range(upcolnum):
				upworksheet.write(row+22+info_addline, col, updata[row][col],zhengwenformat)
		# down data
		dndata = read_data_content_with_float(foldname+'/gene_exp.diff.%s_vs_%s.dn.txt' % (groupAname, groupBname))
		dnsheetname = 'dn.%s_vs_%s' % (groupAname, groupBname)
		dnworksheet = workbook.add_worksheet(dnsheetname[:30])  # dnsheetname
		dncolnum = len(dndata[0])
		dnworksheet.merge_range(0, 0, 19+info_addline, dncolnum-1, '', annotationformat)
		dnworksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname+'.down','\nFold Change: ', red, fc, '\n','p-value: ', red, p,  '\nFPKM >= ', red, str(FACTOR), ' in at least one group',  excel_title_info_two_group_de_gene, annotationformat)
		dndata_title = dndata[0]
		dnrownum = len(dndata)
		dnworksheet.merge_range(21+info_addline, 0, 21+info_addline, dncolnum-1, '%s vs %s down-regulated lncRNA' % (groupAname, groupBname), dntitleformat)
		for ii in range(dncolnum):
			dnworksheet.write(0+22+info_addline,ii,dndata_title[ii].replace('q_value', 'FDR'), titleformat)
		for row2 in range(1, dnrownum):
			for col2 in range(dncolnum):
				dnworksheet.write(row2+22+info_addline, col2, dndata[row2][col2],zhengwenformat)
	workbook.close()
	
	# for protein_coding
	excel_title_info_two_group_de_gene_protein_coding = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, protein_coding.
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: fold change, fold change between two groups (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two groups (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two groups.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ V: annotation information, including GeneID, Synonyms, dbXrefs (other database ID), chromosome, map_location, description, pathway and GO annotations.

note: since the different sources, some genes without annotation information"""

	# write into excel 
	workbook2 = Workbook('Differentially Expressed mRNAs (group).xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 =  workbook2.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	uptitleformat2 = workbook2.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat2 = workbook2.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2group_comparisons()
	#print comparisons
	for comparison2 in comparisons2:
		groupA2, groupAname2, groupB2, groupBname2, fc2, p2 = comparison2
		abfc2 = float(fc2)
		abp2 = float(p2)
		foldname2 = 'cuffdiff_'+ groupAname2+'_vs_'+groupBname2
		# up2 data
		updata2 = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.up.protein_coding.txt' % (groupAname2, groupBname2))
		upsheetname2 = 'up.%s_vs_%s' % (groupAname2, groupBname2)
		upworksheet2 = workbook2.add_worksheet(upsheetname2[0:30])  # upsheetname2
		upcolnum2 = len(updata2[0])
		upworksheet2.merge_range(0, 0, 21+info_addline, upcolnum2-1, '', annotationformat2)
		upworksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.up','\nFold Change: ', red2, fc2, '\n','p-value: ', red2, p2, '\nFPKM >= ', red2, str(FACTOR), ' in at least one group', excel_title_info_two_group_de_gene_protein_coding, annotationformat2)
		updata_title2 = updata2[0]
		uprownum2 = len(updata2)
		upworksheet2.merge_range(23+info_addline, 0, 23+info_addline, upcolnum2-1, '%s vs %s up-regulated genes' % (groupAname2, groupBname2), uptitleformat2)
		for i2 in range(upcolnum2):
			upworksheet2.write(0+24+info_addline,i2,updata_title2[i2].replace('q_value', 'FDR'), titleformat2)
		for row2 in range(1, uprownum2):
			for col2 in range(upcolnum2):
				upworksheet2.write(row2+24+info_addline, col2, updata2[row2][col2])
		# down2 data
		dndata2 = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.dn.protein_coding.txt' % (groupAname2, groupBname2))
		dnsheetname2 = 'dn.%s_vs_%s' % (groupAname2, groupBname2)
		dnworksheet2 = workbook2.add_worksheet(dnsheetname2[:30])  # dnsheetname2
		dncolnum2 = len(dndata2[0])
		dnworksheet2.merge_range(0, 0, 21+info_addline, dncolnum2-1, '', annotationformat2)
		dnworksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.down','\nFold Change: ', red2, fc2, '\n','p-value: ', red2, p2, '\nFPKM >= ', red2, str(FACTOR), ' in at least one group',excel_title_info_two_group_de_gene_protein_coding, annotationformat2)
		dndata_title2 = dndata2[0]
		dnrownum2 = len(dndata2)
		dnworksheet2.merge_range(23+info_addline, 0, 23+info_addline, dncolnum2-1, '%s vs %s down-regulated genes' % (groupAname2, groupBname2), dntitleformat2)
		for ii2 in range(dncolnum2):
			dnworksheet2.write(0+24+info_addline,ii2,dndata_title2[ii2].replace('q_value', 'FDR'), titleformat2)
		for row22 in range(1, dnrownum2):
			for col22 in range(dncolnum2):
				dnworksheet2.write(row22+24+info_addline, col22, dndata2[row22][col22])
	workbook2.close()

###################################### 28 two group de gene (excel, lnc+protein, protein), end


###################################### 29 two group de transcript (txt), start

def two_group_comparison_de_transcript():
	"""
	similar to two_sample_comparison_de_transcript
	"""
	comparisons = get_2group_comparisons()
	#print comparisons
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	transcript_biotype_strand_dict = parse_gtf_transcript(hisat2_gtf)
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		abfc = float(fc)
		abp = float(p)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		f =open(foldname+'/isoform_exp.diff')
		fup = open(foldname+'/isoform_exp.diff.%s_vs_%s.up.txt' % (groupAname, groupBname),'w')
		fdn = open(foldname+'/isoform_exp.diff.%s_vs_%s.dn.txt' % (groupAname, groupBname),'w')
		flncup = open(foldname+'/isoform_exp.diff.%s_vs_%s.up.lnctmp.txt' % (groupAname, groupBname),'w')
		flncdn = open(foldname+'/isoform_exp.diff.%s_vs_%s.dn.lnctmp.txt' % (groupAname, groupBname),'w')
		lines = f.readlines()
		headlist = lines[0].strip('\n').split('\t')
		print >>fup, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>fdn, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>flncup, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>flncdn, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == groupAname and linelist[4] == groupBname:
				if linelist[9] =='inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] =='-inf':
					if float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					if float(linelist[9]) >= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					elif float(linelist[9]) <= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc= log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:
						continue
			#todo
			elif linelist[4] == groupAname and linelist[5] == groupBname:
				if linelist[9] == 'inf':
					if float(linelist[11]) <=abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				elif linelist[9] =='-inf':
					if float(linelist[11]) <=abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					if float(linelist[9])>= numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='down'
							print >>flncdn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					elif float(linelist[9])<= -numpy.log2(abfc) and float(linelist[11]) <= abp and at_least_one_expressed([linelist[8], linelist[7]]):
						print >>fup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
							newfc = log2fc2fc(linelist[9])
							regulation='up'
							print >>flncup, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:
						continue
			else:
				continue
		f.close()
		fup.close()
		fdn.close()
		flncup.close()
		flncdn.close()
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.up.lnctmp.txt' % (groupAname, groupBname), lnc_associated_genefile)
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.dn.lnctmp.txt' % (groupAname, groupBname), lnc_associated_genefile)
		
		
###################################### 29 two group de transcript (lnc+protein, lnc, txt), over

###################################### 30 two group de transcript (excel), start

def excel_two_group_comparison_de_transcript(): 
	"""
	isoform_exp.diff.G2_vs_G1.txt
	"""
	# gene	gene_id	biotype	strand	locus	6W-KO-M-2_FPKM	6W-KO-FM-2_FPKM	log2(fold change)	p_value	q_value
	two_group_comparison_de_transcript()
	excel_title_info_two_group_de_transcript = """

Column A: transcript_id, the transcript identifier.
Column B: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: log2(fold change), fold change between two groups, differentially expressed lncRNAs were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two groups.
Column J: FDR, adjusted p-value."""

	# write into excel 
	workbook = Workbook('Differentially Expressed (lnc+protein,group).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align':  'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	uptitleformat = workbook.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat = workbook.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})	
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2group_comparisons()
	#print comparisons
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		abfc = float(fc)
		abp = float(p)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		# up data
		updata = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.up.txt' % (groupAname, groupBname))
		upsheetname = 'up.%s_vs_%s' % (groupAname, groupBname)
		upworksheet = workbook.add_worksheet(upsheetname[:30])  # upsheetname
		upcolnum = len(updata[0])
		upworksheet.merge_range(0, 0, 19+info_addline, upcolnum-1, '', annotationformat)
		upworksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname,'\nFold Change: ', red, fc, '\n','p-value: ', red, p , '\nFPKM >= ', red, str(FACTOR), ' in at least one group',excel_title_info_two_group_de_transcript, annotationformat)
		updata_title = updata[0]
		uprownum = len(updata)
		upworksheet.merge_range(21+info_addline, 0, 21+info_addline, upcolnum-1, '%s vs %s up-regulated lncRNAs' % (groupAname, groupBname), uptitleformat)
		for i in range(upcolnum):
			upworksheet.write(0+22+info_addline,i,updata_title[i].replace('q_value', 'FDR'), titleformat)
		for row in range(1, uprownum):
			for col in range(upcolnum):
				upworksheet.write(row+22+info_addline, col, updata[row][col],zhengwenformat)
		# down data
		dndata = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.dn.txt' % (groupAname, groupBname))
		dnsheetname = 'dn.%s_vs_%s' % (groupAname, groupBname)
		dnworksheet = workbook.add_worksheet(dnsheetname[:30])  # dnsheetname
		dncolnum = len(dndata[0])
		dnworksheet.merge_range(0, 0, 19+info_addline, dncolnum-1, '', annotationformat)
		dnworksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname,'\nFold Change: ', red, fc, '\n','p-value: ', red, p , '\nFPKM >= ', red, str(FACTOR), ' in at least one group',excel_title_info_two_group_de_transcript, annotationformat)
		dndata_title = dndata[0]
		dnrownum = len(dndata)
		dnworksheet.merge_range(21+info_addline, 0, 21+info_addline, dncolnum-1, '%s vs %s down-regulated lncRNAs' % (groupAname, groupBname), dntitleformat)
		for ii in range(dncolnum):
			dnworksheet.write(0+22+info_addline,ii,dndata_title[ii].replace('q_value', 'FDR'), titleformat)
		for row2 in range(1, dnrownum):
			for col2 in range(dncolnum):
				dnworksheet.write(row2+22+info_addline, col2, dndata[row2][col2],zhengwenformat)
	workbook.close()
	
	######## for lncRNA
	excel_title_info_two_group_de_transcript_lnc = """

Column A: transcript_id, the transcript id.
Column B: biotype, type of gene, including antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: fold change, fold change between two groups (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two groups (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two groups.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ U: the genomic organization of LncRNAs and the information of the associated coding genes.

the classes of LncRNAs defined by genomic organization:
"exon sense-overlapping": the LncRNA's exon is overlapping a coding transcript exon on the same genomic strand;
"intron sense-overlapping": the LncRNA is overlapping the intron of a coding transcript on the same genomic strand;
"intronic antisense": the LncRNA is overlapping the intron of a coding transcript on the antisense strand;
"natural antisense": the LncRNA is transcribed from the antisense strand and overlapping with a coding transcript; 
"bidirectional": the LncRNA is oriented head to head to a coding transcript within 1000 bp; 
"intergenic": there are no overlapping or bidirectional coding transcripts nearby the LncRNA."""

	# write into excel
	workbook2 = Workbook('Differentially Expressed LncRNAs (group).xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 =  workbook2.add_format({'align':  'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	uptitleformat2 = workbook2.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	dntitleformat2 = workbook2.add_format({'bold':True,'fg_color':'green','align':'center','font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2group_comparisons()
	#print comparisons
	for comparison2 in comparisons2:
		groupA2, groupAname2, groupB2, groupBname2, fc2, p2 = comparison2
		abfc2 = float(fc2)
		abp2 = float(p2)
		foldname2 = 'cuffdiff_'+ groupAname2+'_vs_'+groupBname2
		# up data
		updata2 = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.up.lnc.txt' % (groupAname2, groupBname2))
		upsheetname2 = 'up.%s_vs_%s' % (groupAname2, groupBname2)
		upworksheet2 = workbook2.add_worksheet(upsheetname2[:30])  # upsheetname
		upcolnum2 = len(updata2[0])
		upworksheet2.merge_range(0, 0, 20+info_addline, upcolnum2-1, '', annotationformat2)
		upworksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.up','\nFold Change: ', red2, fc2, '\n','p-value: ', red2, p2, '\nFPKM >= ', red2, str(FACTOR), ' in at least one group',excel_title_info_two_group_de_transcript_lnc, annotationformat2)
		updata_title2 = updata2[0]
		uprownum2 = len(updata2)
		upworksheet2.merge_range(22+info_addline, 0, 22+info_addline, upcolnum2-1, '%s vs %s up-regulated lncRNAs' % (groupAname2, groupBname2), uptitleformat2)
		for i2 in range(upcolnum2):
			upworksheet2.write(0+23+info_addline, i2, updata_title2[i2].replace('q_value', 'FDR'), titleformat2)
		for row2 in range(1, uprownum2):
			for col2 in range(upcolnum2):
				upworksheet2.write(row2+23+info_addline, col2, updata2[row2][col2],zhengwenformat2)
		# down data
		dndata2 = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.dn.lnc.txt' % (groupAname2, groupBname2))
		dnsheetname2 = 'dn.%s_vs_%s' % (groupAname2, groupBname2)
		dnworksheet2 = workbook2.add_worksheet(dnsheetname2[:30])  # dnsheetname
		dncolnum2 = len(dndata2[0])
		dnworksheet2.merge_range(0, 0, 20+info_addline, dncolnum2-1, '', annotationformat2)
		dnworksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.down','\nFold Change: ', red2, fc2, '\n','p-value: ', red2, p2, '\nFPKM >= ', red2, str(FACTOR), ' in at least one group',excel_title_info_two_group_de_transcript_lnc, annotationformat2)
		dndata_title2 = dndata2[0]
		dnrownum2 = len(dndata2)
		dnworksheet2.merge_range(22+info_addline, 0, 22+info_addline, dncolnum2-1, '%s vs %s down-regulated lncRNAs' % (groupAname2, groupBname2), dntitleformat2)
		for ii2 in range(dncolnum2):
			dnworksheet2.write(0+23+info_addline,ii2, dndata_title2[ii2].replace('q_value', 'FDR'), titleformat2)
		for row22 in range(1, dnrownum2):
			for col22 in range(dncolnum2):
				dnworksheet2.write(row22+23+info_addline, col22, dndata2[row22][col22],zhengwenformat2)
	workbook2.close()
	
###################################### 30 two group de transcript (excel, lncRNA+protein, lnc), start
###### 31 All comparison gene expression excel (group), start

def excel_two_group_comparison_all_gene(): 
	"""
	gene_exp.diff.G2_vs_G1.txt
	"""
	# gene	gene_id	biotype	strand	locus	6W-KO-M-2_FPKM	6W-KO-FM-2_FPKM	log2(fold change)	p_value	q_value
	excel_title_info_two_group_all_comparison = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: fold change, fold change between two groups (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two groups (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two groups.
Column K: FDR, adjusted p-value."""

	# write into excel
	two_group_comparison_gene()  # all and protein_coding
	workbook = Workbook('All Comparison (lnc+protein,group).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2group_comparisons()
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		data = read_data_content_with_float(foldname+'/gene_exp.diff.%s_vs_%s.txt' % (groupAname, groupBname))
		sheetname = '%s_vs_%s' % (groupAname, groupBname)
		worksheet = workbook.add_worksheet(sheetname[:30])  # sheetname
		colnum = len(data[0])
		worksheet.merge_range(0, 0, 19+info_addline, colnum-1, '', annotationformat)
		worksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname, excel_title_info_two_group_all_comparison, annotationformat)
		data_title = data[0]
		rownum = len(data)
		for i in range(colnum):
			worksheet.write(0+22+info_addline,i,data_title[i], titleformat)
		for row in range(1, rownum):
			for col in range(colnum):
				worksheet.write(row+22+info_addline, col, data[row][col],zhengwenformat)
	workbook.close()
	
	# for protein_coding
	excel_title_info_two_group_all_comparison_protein_coding = """

Column A: gene, the gene name.
Column B: gene_id, the Ensembl gene identifier.
Column C: biotype, protein_coding.
Column D: strand, transcription direction.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: fold change, fold change between two groups (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two groups (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two groups.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ V: annotation information, including GeneID, Synonyms, dbXrefs (other database ID), chromosome, map_location and description.

note: since the different sources, some genes without annotation information"""

	# write into excel
	workbook2 = Workbook('All Comparison (mRNA, group).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat2 =  workbook2.add_format({'align': 'left',  'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'})
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	
	comparisons2 = get_2group_comparisons()
	for comparison2 in comparisons2:
		groupA2, groupAname2, groupB2, groupBname2, fc2, p2 = comparison2
		foldname2 = 'cuffdiff_'+ groupAname2+'_vs_'+groupBname2
		
		# testing
		# up
		data2up = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.protein_coding.up.txt' % (groupAname2, groupBname2))
		sheetname2up = 'up.%s_vs_%s' % (groupAname2, groupBname2)
		worksheet2 = workbook2.add_worksheet(sheetname2up[:30])  # sheetname2
		colnum2up = len(data2up[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2up-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.up', excel_title_info_two_group_all_comparison_protein_coding, annotationformat2)
		data_title2up = data2up[0]
		rownum2up = len(data2up)
		for i2up in range(colnum2up):
			worksheet2.write(0+22+info_addline,i2up,data_title2up[i2up].replace('q_value', 'FDR'), titleformat2)
		for row2up in range(1, rownum2up):
			for col2up in range(colnum2up):
				worksheet2.write(row2up+22+info_addline, col2up, data2up[row2up][col2up],zhengwenformat2)
				
		# down
		
		data2dn = read_data_content_with_float(foldname2+'/gene_exp.diff.%s_vs_%s.protein_coding.dn.txt' % (groupAname2, groupBname2))
		sheetname2dn = 'dn.%s_vs_%s' % (groupAname2, groupBname2)
		worksheet2 = workbook2.add_worksheet(sheetname2dn[:30])  # sheetname2
		colnum2dn = len(data2dn[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2dn-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.down', excel_title_info_two_group_all_comparison_protein_coding, annotationformat2)
		data_title2dn = data2dn[0]
		rownum2dn = len(data2dn)
		for i2dn in range(colnum2dn):
			worksheet2.write(0+22+info_addline,i2dn,data_title2dn[i2dn].replace('q_value', 'FDR'), titleformat2)
		for row2dn in range(1, rownum2dn):
			for col2dn in range(colnum2dn):
				worksheet2.write(row2dn+22+info_addline, col2dn, data2dn[row2dn][col2dn],zhengwenformat2)
	workbook2.close()


###################################### 31 All comparison gene expression excel (protein+lnc, protein, group), over


###################################### 32 DE of transcript expression txt (group), start

def two_group_comparison_transcript():
	"""
	only lncRNA transcript levels
	"""
	hisat2_gtf = cf.get('mapping_parameters','hisat2_gtf')
	transcript_biotype_strand_dict = parse_gtf_transcript(hisat2_gtf)
	comparisons = get_2group_comparisons()
	#print comparisons
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		abfc = float(fc)
		abp = float(p)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		f =open(foldname+'/isoform_exp.diff')
		f2 = open(foldname+'/isoform_exp.diff.%s_vs_%s.txt' % (groupAname, groupBname),'w')
		flnc2 = open(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.txt' % (groupAname, groupBname),'w') # lnc
		flnc2up = open(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.up.txt' % (groupAname, groupBname),'w') # lnc.up
		flnc2dn = open(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.dn.txt' % (groupAname, groupBname),'w') # lnc.dn
		lines = f.readlines()
		f.close()
		headlist = lines[0].strip('\n').split('\t')
		print >>f2, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>flnc2, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]
		print >>flnc2up, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		print >>flnc2dn, headlist[0].replace('test_id','transcript_id')+'\tbiotype\tstrand\t'+headlist[1]+'\t'+headlist[3]+'\t'+groupAname+'_FPKM'+'\t'+groupBname+'_FPKM'+'\t'+'Fold change'+'\t'+'log2(fold change)'+'\t'+headlist[11]+'\t'+headlist[12]+'\t'+'Regulation'
		for line in lines:
			linelist = line.strip('\n').split('\t')
			if linelist[5] == groupAname and linelist[4] == groupBname:  # 样品顺序不一样，fpkm要反过来， logfc不变 # 貌似这里不用
				print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
				if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
					print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]

					if float(linelist[9])>=0 or linelist[9]=='inf':   # >=0或者是inf的归为上调
						newfc = log2fc2fc(linelist[9])
						regulation = 'up'
						print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
					else:									  # 其他的归为下调（即<0或者-inf）
						newfc = log2fc2fc(linelist[9])
						regulation = 'down'
						print >>flnc2dn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[8]+'\t'+linelist[7]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation

			elif linelist[4] == groupAname and linelist[5] == groupBname: #样品顺序一样，则fpkm不变，logfc取反
				if linelist[9] == 'inf':
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='down'
						print >>flnc2dn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'-inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						
				elif linelist[9] == '-inf':
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='up'
						print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+'inf'+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						
				elif linelist[9] == '0':
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]
						newfc = log2fc2fc(linelist[9])
						regulation='up'
						print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+linelist[9]+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
				else:
					print >>f2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
					if transcript_biotype_strand_dict[linelist[0]].split()[0] in lncRNAclasses:
						print >>flnc2, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]
						if float(linelist[9])<0:
							newfc = log2fc2fc(linelist[9])
							regulation = 'up'
							print >>flnc2up, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
						else:
							newfc = log2fc2fc(linelist[9])
							regulation = 'down'
							print >>flnc2dn, linelist[0]+'\t'+transcript_biotype_strand_dict[linelist[0]]+'\t'+linelist[1]+'\t'+linelist[3]+'\t'+linelist[7]+'\t'+linelist[8]+'\t'+str(newfc)+'\t'+str(-float(linelist[9]))+'\t'+linelist[11]+'\t'+linelist[12]+'\t'+regulation
			else:
				continue
		f2.close()
		flnc2.close()
		flnc2up.close() # lncRNA up
		flnc2dn.close() # lncRNA down
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.txt' % (groupAname, groupBname), lnc_associated_genefile)
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.up.txt' % (groupAname, groupBname), lnc_associated_genefile)
		lncRNAaddgene(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.dn.txt' % (groupAname, groupBname), lnc_associated_genefile)


###################################### 32 DE of transcript expression txt (group), over

###################################### 33 All comparison of transcript expression excel (group), start
def excel_two_group_comparison_all_transcript():
	"""
	isoform_exp.diff.G2_vs_G1.txt
	"""
	# test_id	biotype	strand	gene_id	locus	6W-KO-M-2_FPKM	6W-WT-FM-1_FPKM	log2(fold change)	p_value	q_value
	excel_title_info_two_group_all_comparison = """

Column A: transcript_id, the Ensembl transcript id.
Column B: biotype, type of gene, including protein_coding, pseudogene, antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: log2(fold change), fold change between two groups, differentially expressed lncRNAs were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two groups.
Column J: FDR, adjusted p-value."""

	# write into excel 
	two_group_comparison_transcript()
	workbook = Workbook('All Comparison (group, mRNA+lnc).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat = workbook.add_format({'bold':True,'font_name':'Times New Roman'}) 
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2group_comparisons()
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		abfc = float(fc)
		abp = float(p)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		data = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.txt' % (groupAname, groupBname))
		sheetname = '%s_vs_%s' % (groupAname, groupBname)
		worksheet = workbook.add_worksheet(sheetname[:30])  # sheetname
		colnum = len(data[0])
		worksheet.merge_range(0, 0, 19+info_addline, colnum-1, '', annotationformat)
		worksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname, excel_title_info_two_group_all_comparison, annotationformat)
		data_title = data[0]
		rownum = len(data)
		for i in range(colnum):
			worksheet.write(0+22+info_addline,i,data_title[i], titleformat)
		for row in range(1, rownum):
			for col in range(colnum):
				worksheet.write(row+22+info_addline, col, data[row][col], zhengwenformat)
	workbook.close()
	##### for lncRNA
	excel_title_info_two_group_all_comparison_lnc = """

Column A: transcript_id, the Ensembl transcript id.
Column B: biotype, type of gene, including antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each group.
Column H: fold change, fold change between two groups (inf: infinity, -inf: negative infinity).
Column I: log2(fold change), fold change between two groups (inf: infinity, -inf: negative infinity).
Column J: p_value, p-value between two groups.
Column K: FDR, adjusted p-value.
Column L: Regulation, up indicates up-regulation, and down indicates down-regulation.
Column M ~ U: the genomic organization of LncRNAs and the information of the associated coding genes.

the classes of LncRNAs defined by genomic organization:
"exon sense-overlapping": the LncRNA's exon is overlapping a coding transcript exon on the same genomic strand;
"intron sense-overlapping": the LncRNA is overlapping the intron of a coding transcript on the same genomic strand;
"intronic antisense": the LncRNA is overlapping the intron of a coding transcript on the antisense strand;
"natural antisense": the LncRNA is transcribed from the antisense strand and overlapping with a coding transcript; 
"bidirectional": the LncRNA is oriented head to head to a coding transcript within 1000 bp; 
"intergenic": there are no overlapping or bidirectional coding transcripts nearby the LncRNA."""


	# write into excel 
	workbook2 = Workbook('All Comparison (lncRNA, group).xlsx')
	# the following is xlsx properties, press right key --> summary
	annotationformat2 =  workbook2.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat2.set_text_wrap()  # auto wrapping
	red2 = workbook2.add_format({'color':'red','bold':True,'font_name':'Times New Roman'})
	titleformat2 = workbook2.add_format({'bold':True,'font_name':'Times New Roman'}) 
	zhengwenformat2 =  workbook2.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons2 = get_2group_comparisons()
	for comparison2 in comparisons2:
		groupA2, groupAname2, groupB2, groupBname2, fc2, p2 = comparison2
		abfc2 = float(fc2)
		abp2 = float(p2)
		foldname2 = 'cuffdiff_'+ groupAname2+'_vs_'+groupBname2
		
		# up, testing
		data2up = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.lnc.up.txt' % (groupAname2, groupBname2))
		sheetname2up = 'up.%s_vs_%s' % (groupAname2, groupBname2)
		worksheet2 = workbook2.add_worksheet(sheetname2up[:30])  # sheetname
		colnum2up = len(data2up[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2up-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.up', excel_title_info_two_group_all_comparison_lnc, annotationformat2)
		data_title2up = data2up[0]
		rownum2up = len(data2up)
		for i2up in range(colnum2up):
			worksheet2.write(0+22+info_addline, i2up, data_title2up[i2up].replace('q_value', 'FDR'), titleformat2)
		for row2up in range(1, rownum2up):
			for col2up in range(colnum2up):
				worksheet2.write(row2up+22+info_addline, col2up, data2up[row2up][col2up],zhengwenformat2)
		
		# dn, testing
		data2dn = read_data_content_with_float(foldname2+'/isoform_exp.diff.%s_vs_%s.lnc.dn.txt' % (groupAname2, groupBname2))
		sheetname2dn = 'dn.%s_vs_%s' % (groupAname2, groupBname2)
		worksheet2 = workbook2.add_worksheet(sheetname2dn[:30])  # sheetname
		colnum2dn = len(data2dn[0])
		worksheet2.merge_range(0, 0, 19+info_addline, colnum2dn-1, '', annotationformat2)
		worksheet2.write_rich_string(0, 0, red2, groupAname2+' vs '+groupBname2+'.down', excel_title_info_two_group_all_comparison_lnc, annotationformat2)
		data_title2dn = data2dn[0]
		rownum2dn = len(data2dn)
		for i2dn in range(colnum2dn):
			worksheet2.write(0+22+info_addline, i2dn, data_title2dn[i2dn].replace('q_value', 'FDR'), titleformat2)
		for row2dn in range(1, rownum2dn):
			for col2dn in range(colnum2dn):
				worksheet2.write(row2dn+22+info_addline, col2dn, data2dn[row2dn][col2dn],zhengwenformat2)
	workbook2.close()

###################################### 33 All comparison of transcript expression excel (group), over


###############plots

#scatter


def two_samples_scatterplot_protein():
	# for two samples protein, note: scatter (B, A) generate y-axis A, x-axis B
	two_sample_comparisons = get_2_samples_comparisons()
	for two_sample_comparison in two_sample_comparisons:
		sampleA = two_sample_comparison[0]
		sampleB = two_sample_comparison[1]
		fc = float(two_sample_comparison[2])
		# gene_exp.diff.33_vs_S0.txt
		foldname = 'cuffdiff_'+ sampleA.split('_sequence.fastq')[0]+'_vs_'+sampleB.split('_sequence.fastq')[0]
		filename = foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.txt' % (sampleA.replace('_sequence.fastq',''), sampleB.replace('_sequence.fastq',''))
		f = open(filename)
		headlist = f.readline().strip('\n').split('\t')
		sampleAname = headlist[5].replace('_FPKM','')  # 33
		sampleBname = headlist[6].replace('_FPKM','')  # S0
		A = []
		B = []
		colors=[]
		for line in f:
			linelist = line.strip('\n').split('\t')
			if [linelist[6], linelist[5]].count('0')!=2:
				a = float(linelist[5])+1.0
				b = float(linelist[6])+1.0
				A.append(numpy.log2(a))
				B.append(numpy.log2(b))
				if a/b >=fc:
					colors.append('#C6383F')
				elif a/b <= 1.0/fc:
					colors.append('#79AC44')
				else:
					colors.append('#7D74B2')
		f.close()
		fig = figure(figsize=cm2inch(17,17))
		#grid(linestyle='solid',lw=0.5, color='#606060')
		scatter(B, A, marker='o',linewidths=0.75,color= colors, s=4) #, edgecolors='#F0F0F0')
		maxn = max(max(B), max(A))
		edge = round(maxn)+2
		plot([0,edge],[0,edge], c='#606060', lw='1.5')
		x= [0,edge]
		plot(x,[0+numpy.log2(fc),edge+numpy.log2(fc)], c='#606060', lw='1.5')  # up
		plot(x,[0-numpy.log2(fc),edge-numpy.log2(fc)], c='#606060', lw='1.5')  # down
		xlim(0, edge)
		ylim(0, edge)
		ylabel(sampleAname, fontsize=20)  # 6
		xlabel(sampleBname, fontsize=20)  # 7
		savefig(sampleA.replace('_sequence.fastq','')+'_vs_'+sampleB.replace('_sequence.fastq','')+'.mRNA.png', dpi=600.1)
		png2tiff(sampleA.replace('_sequence.fastq','')+'_vs_'+sampleB.replace('_sequence.fastq','')+'.mRNA.png')
		

def two_samples_scatterplot_lnc():
	# for two samples protein, note: scatter (B, A) generate y-axis A, x-axis B
	two_sample_comparisons = get_2_samples_comparisons()
	for two_sample_comparison in two_sample_comparisons:
		sampleA = two_sample_comparison[0]
		sampleB = two_sample_comparison[1]
		fc = float(two_sample_comparison[2])
		# gene_exp.diff.B2_vs_A2.lnctmp.txt
		foldname = 'cuffdiff_'+ sampleA.split('_sequence.fastq')[0]+'_vs_'+sampleB.split('_sequence.fastq')[0]
		filename = foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.txt' % (sampleA.replace('_sequence.fastq',''), sampleB.replace('_sequence.fastq',''))
		f = open(filename)
		headlist = f.readline().strip('\n').split('\t')
		sampleAname = headlist[5].replace('_FPKM','')  # 33
		sampleBname = headlist[6].replace('_FPKM','')  # S0
		A = []
		B = []
		colors=[]
		for line in f:
			linelist = line.strip('\n').split('\t')
			if [linelist[6], linelist[5]].count('0')!=2:
				a = float(linelist[5])+1.0
				b = float(linelist[6])+1.0
				A.append(numpy.log2(a))
				B.append(numpy.log2(b))
				if a/b >=fc:
					colors.append('#C6383F')
				elif a/b <= 1.0/fc:
					colors.append('#79AC44')
				else:
					colors.append('#7D74B2')
		f.close()
		fig = figure(figsize=cm2inch(17,17))
		#grid(linestyle='solid',lw=0.5, color='#606060')
		scatter(B, A, marker='o',linewidths=0.75,color= colors, s=4) #, edgecolors='#F0F0F0')
		maxn = max(max(B), max(A))
		edge = round(maxn)+2
		plot([0,edge],[0,edge], c='#606060', lw='1.5')
		x= [0,edge]
		plot(x,[0+numpy.log2(fc),edge+numpy.log2(fc)], c='#606060', lw='1.5')  # up
		plot(x,[0-numpy.log2(fc),edge-numpy.log2(fc)], c='#606060', lw='1.5')  # down
		xlim(0, edge)
		ylim(0, edge)
		ylabel(sampleAname, fontsize=20)  # 6
		xlabel(sampleBname, fontsize=20)  # 7
		savefig(sampleA.replace('_sequence.fastq','')+'_vs_'+sampleB.replace('_sequence.fastq','')+'.LncRNA.png', dpi=600.1)
		png2tiff(sampleA.replace('_sequence.fastq','')+'_vs_'+sampleB.replace('_sequence.fastq','')+'.LncRNA.png')

def two_groups_scatterplot_protein():
	# for two group protein, note: scatter (B, A) generate y-axis A, x-axis B
	comparisons = get_2group_comparisons()
	#print comparisons
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		fc = float(fc)
		# gene_exp.diff.B_vs_A.protein_coding.txt
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		filename = foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.txt' % (groupAname, groupBname)
		f = open(filename)
		headlist = f.readline().strip('\n').split('\t')
		A = []
		B = []
		colors=[]
		for line in f:
			linelist = line.strip('\n').split('\t')
			if [linelist[6], linelist[5]].count('0')!=2:
				a = float(linelist[5])+1.0
				b = float(linelist[6])+1.0
				A.append(numpy.log2(a))
				B.append(numpy.log2(b))
				if a/b >=fc:
					colors.append('#C6383F')
				elif a/b <= 1.0/fc:
					colors.append('#79AC44')
				else:
					colors.append('#7D74B2')
		f.close()
		fig = figure(figsize=cm2inch(17,17))
		#grid(linestyle='solid',lw=0.5, color='#606060')
		scatter(B, A, marker='o',linewidths=0.75,color= colors, s=4) #, edgecolors='#F0F0F0')
		maxn = max(max(B), max(A))
		edge = round(maxn)+2
		plot([0,edge],[0,edge], c='#606060', lw='1.5')
		x= [0,edge]
		plot(x,[0+numpy.log2(fc),edge+numpy.log2(fc)], c='#606060', lw='1.5')  # up
		plot(x,[0-numpy.log2(fc),edge-numpy.log2(fc)], c='#606060', lw='1.5')  # down
		xlim(0, edge)
		ylim(0, edge)
		ylabel(groupAname, fontsize=20)  # 6
		xlabel(groupBname, fontsize=20)  # 7
		savefig(groupAname.replace('_sequence.fastq','')+'_vs_'+groupBname.replace('_sequence.fastq','')+'.mRNA.png', dpi=600.1)
		png2tiff(groupAname.replace('_sequence.fastq','')+'_vs_'+groupBname.replace('_sequence.fastq','')+'.mRNA.png')

def two_groups_scatterplot_lnc():
	# for two group protein, note: scatter (B, A) generate y-axis A, x-axis B
	comparisons = get_2group_comparisons()
	#print comparisons
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		fc = float(fc)
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		# gene_exp.diff.33_vs_S0.txt
		filename = foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.txt' % (groupAname, groupBname)
		f = open(filename)
		headlist = f.readline().strip('\n').split('\t')
		A = []
		B = []
		colors=[]
		for line in f:
			linelist = line.strip('\n').split('\t')
			if [linelist[6], linelist[5]].count('0')!=2:
				a = float(linelist[5])+1.0
				b = float(linelist[6])+1.0
				A.append(numpy.log2(a))
				B.append(numpy.log2(b))
				if a/b >=fc:
					colors.append('#C6383F')
				elif a/b <= 1.0/fc:
					colors.append('#79AC44')
				else:
					colors.append('#7D74B2')
		f.close()
		fig = figure(figsize=cm2inch(17,17))
		#grid(linestyle='solid',lw=0.5, color='#606060')
		scatter(B, A, marker='o',linewidths=0.75,color= colors, s=4) #, edgecolors='#F0F0F0')
		maxn = max(max(B), max(A))
		edge = round(maxn)+2
		plot([0,edge],[0,edge], c='#606060', lw='1.5')
		x= [0,edge]
		plot(x,[0+numpy.log2(fc),edge+numpy.log2(fc)], c='#606060', lw='1.5')  # up
		plot(x,[0-numpy.log2(fc),edge-numpy.log2(fc)], c='#606060', lw='1.5')  # down
		xlim(0, edge)
		ylim(0, edge)
		ylabel(groupAname, fontsize=20)  # 6
		xlabel(groupBname, fontsize=20)  # 7
		savefig(groupAname.replace('_sequence.fastq','')+'_vs_'+groupBname.replace('_sequence.fastq','')+'.LncRNA.png', dpi=600.1)
		png2tiff(groupAname.replace('_sequence.fastq','')+'_vs_'+groupBname.replace('_sequence.fastq','')+'.LncRNA.png')

def aligned_pairs(file):
	f = open(file)
	line = f.readlines()[-4]
	f.close()
	return line.strip().split()[-1]


def reads_stat(files):
	f=open('read.tongji.xls','a')
	mapping = overall_mapping_rates(files)
	for fqfile in files:
		n = fastq_lines_by_4(fqfile[:-6]+'_1.fastq')
		cutadaptn= fastq_lines_by_4(fqfile[:-6]+'.cutadapt_1.fastq')  # since 1 & 2 are the same
		rate = mapping[fqfile]
		print >>f, '%s\t%d\t%d\t%s' % (fqfile[:-15], n*2, cutadaptn*2, rate)
	f.close()


# both DE add 20160316, remove 20170328
def merge2files(file1, file2):
	f = open(file1)
	all1 = f.read()
	f.close()
	f2=open(file2)
	head=f2.readline()
	all2=f2.read()
	f3=open(file1.replace('.up','.updn'),'w')
	print >>f3, all1+all2,
	f3.close()

def two_sample_comparison_de_transcript_de_genes(): # both DE add 20160316, remove 20170328
	"""
test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
ENSMUST00000000001	ENSMUSG00000000001	Gnai3	chr3:108107279-108146146	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0.175945	inf	0	1	1	no
ENSMUST00000000003	ENSMUSG00000000003	Pbsn	chrX:77837900-77853623	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0	0	0	1	1	no
	
	The name order (A, B) was the same to value order (B, A), but the logfc is opposite
	so if you want to keep logfc, the sample name order should be reversed,
	or just add minus in front of the logfc
	"""
	
	comparisons = get_2_samples_comparisons()  # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		# merge uplnc and dnlnc
		merge2files(foldname+'/isoform_exp.diff.%s_vs_%s.up.lnc.txt' % (samplea, sampleb), foldname+'/isoform_exp.diff.%s_vs_%s.dn.lnc.txt' % (samplea, sampleb))
		lnccomparisonfile = foldname+'/isoform_exp.diff.%s_vs_%s.updn.lnc.txt' % (samplea, sampleb)
		# merge upprotein and dnprotein
		merge2files(foldname+'/gene_exp.diff.%s_vs_%s.up.protein_coding.txt' % (samplea, sampleb), foldname+'/gene_exp.diff.%s_vs_%s.dn.protein_coding.txt' % (samplea, sampleb))
		mRNAcomparisonfile= foldname+'/gene_exp.diff.%s_vs_%s.updn.protein_coding.txt' % (samplea, sampleb)
		fbothde = open(foldname+'/isoform_exp.diff.%s_vs_%s.bothde.txt' % (samplea, sampleb), 'w')
		fbothdego = open(foldname+'/go.LncRNA.%s_vs_%s.txt' % (samplea, sampleb), 'w')
		fbothdepw = open(foldname+'/pathway.LncRNA.%s_vs_%s.txt' % (samplea, sampleb), 'w')
		# DE mRNAs
		DEgenes = {}
		fmRNA = open(mRNAcomparisonfile)
		headmRNA = fmRNA.readline().strip('\n')
		for line in fmRNA:
			linelist = line.strip('\n').split('\t')
			gene = linelist[0].upper()
			#DEgenes[gene]='\t'.join(linelist[1:])
			DEgenes[gene]=line.strip('\n') # for test
		fmRNA.close()
		# lncRNA
		flnc = open(lnccomparisonfile)
		headlnc = flnc.readline().strip('\n')
		print >>fbothde, headlnc+'\t'+headmRNA
		for line2 in flnc:
			line2list = line2.strip('\n').split('\t')
			if line2list[-4].upper() in DEgenes.keys():
				print >>fbothde, line2.strip('\n')+'\t'+DEgenes[line2list[-4].upper()]
				print >>fbothdego, line2list[-4]
				print >>fbothdepw, line2list[-4]+'\torange'
		flnc.close()
		fbothdego.close()
		fbothdepw.close()

def excel_two_sample_comparison_de_transcript_de_genes(): # both DE add 20160316, remove 20170328
	"""
	cuffdiff/isoform_exp.diff.%s_vs_%s.bothde.txt
	"""
	two_sample_comparison_de_transcript_de_genes()
	excel_title_info_two_sample_de_transcript = """

Column A: transcript_id, the transcript id.
Column B: biotype, type of lncRNA, including, pseudogene, antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: log2(fold change), fold change between two samples, differentially expressed lncRNAs were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two samples.
Column J: FDR, adjusted p-value.
Column K ~ S: LncRNA associated gene information.
Column T ~ AM: gene information, including statistics information and annotation information."""

	# write into excel 
	workbook = Workbook('BothDE(2sample).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True})
	titleformat = workbook.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	titleformat2 = workbook.add_format({'bold':True, 'align':'center','font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0, 0.05
	for comparison in comparisons:
		samplea = comparison[0].replace('_sequence.fastq', '')
		sampleb = comparison[1].replace('_sequence.fastq', '')
		foldname = 'cuffdiff_'+ samplea.split('_sequence.fastq')[0]+'_vs_'+sampleb.split('_sequence.fastq')[0]
		#data
		data = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.bothde.txt' % (samplea, sampleb))
		sheetname = '%s_vs_%s' % (samplea, sampleb)
		worksheet = workbook.add_worksheet(sheetname[:30])  # sheetname
		colnum = len(data[0])
		worksheet.merge_range(0, 0, 19+6+info_addline, colnum-1, '', annotationformat)
		worksheet.write_rich_string(0, 0, red, samplea+' vs '+sampleb,'\nFold Change cutoff: ', red, comparison[2], '\n','p-value cutoff: ', red, comparison[3], '\nFPKM >= ', red, str(FACTOR), ' in at least one sample', excel_title_info_two_sample_de_transcript, annotationformat)
		data_title = data[0]
		rownum = len(data)
		worksheet.merge_range(21+6+info_addline, 0, 21+6+info_addline, colnum-1, '%s vs %s both differentially' % (samplea, sampleb), titleformat)
		for i in range(colnum):
			worksheet.write(0+22+6+info_addline,i,data_title[i], titleformat2)
		for row in range(1, rownum):
			for col in range(colnum):
				worksheet.write(row+22+6+info_addline, col, data[row][col],zhengwenformat)
	workbook.close()

# DE both group

def two_group_comparison_de_transcript_de_genes(): # both DE add 20160316, remove 20170328
	"""
test_id	gene_id	gene	locus	sample_1	sample_2	status	value_1	value_2	log2(fold_change)	test_stat	p_value	q_value	significant
ENSMUST00000000001	ENSMUSG00000000001	Gnai3	chr3:108107279-108146146	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0.175945	inf	0	1	1	no
ENSMUST00000000003	ENSMUSG00000000003	Pbsn	chrX:77837900-77853623	6W-WT-FM-1	6W-KO-FM-2	NOTEST	0	0	0	0	1	1	no
	
	The name order (A, B) was the same to value order (B, A), but the logfc is opposite
	so if you want to keep logfc, the sample name order should be reversed,
	or just add minus in front of the logfc
	"""
	comparisons = get_2group_comparisons()
	#print comparisons
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		# gene_exp.diff.B_vs_A.protein_coding.txt
		filename = foldname+'/gene_exp.diff.%s_vs_%s.protein_coding.txt' % (groupAname, groupBname)
		# merge uplnc and dnlnc
		merge2files(foldname+'/isoform_exp.diff.%s_vs_%s.up.lnc.txt' % (groupAname, groupBname), foldname+'/isoform_exp.diff.%s_vs_%s.dn.lnc.txt' % (groupAname, groupBname))
		lnccomparisonfile = foldname+'/isoform_exp.diff.%s_vs_%s.updn.lnc.txt' % (groupAname, groupBname)
		# merge upprotein and dnprotein
		merge2files(foldname+'/gene_exp.diff.%s_vs_%s.up.protein_coding.txt' % (groupAname, groupBname), foldname+'/gene_exp.diff.%s_vs_%s.dn.protein_coding.txt' % (groupAname, groupBname))
		mRNAcomparisonfile= foldname+'/gene_exp.diff.%s_vs_%s.updn.protein_coding.txt' % (groupAname, groupBname)
		fbothde = open(foldname+'/isoform_exp.diff.%s_vs_%s.bothde.txt' % (groupAname, groupBname), 'w')
		fbothdego = open(foldname+'/go.LncRNA.%s_vs_%s.txt' % (groupAname, groupBname), 'w')
		fbothdepw = open(foldname+'/pathway.LncRNA.%s_vs_%s.txt' % (groupAname, groupBname), 'w')
		# DE mRNAs
		DEgenes = {}
		fmRNA = open(mRNAcomparisonfile)
		headmRNA = fmRNA.readline().strip('\n')
		for line in fmRNA:
			linelist = line.strip('\n').split('\t')
			gene = linelist[0].upper()
			#DEgenes[gene]='\t'.join(linelist[1:])
			DEgenes[gene]=line.strip('\n').upper() # for test
		fmRNA.close()
		# lncRNA
		flnc = open(lnccomparisonfile)
		headlnc = flnc.readline().strip('\n')
		print >>fbothde, headlnc+'\t'+headmRNA
		for line2 in flnc:
			line2list = line2.strip('\n').split('\t')
			if line2list[-4].upper() in DEgenes.keys():
				print >>fbothde, line2.strip('\n')+'\t'+DEgenes[line2list[-4].upper()]
				print >>fbothdego, line2list[-4]
				print >>fbothdepw, line2list[-4]+'\torange'
		flnc.close()
		fbothdego.close()
		fbothdepw.close()

def excel_two_group_comparison_de_transcript_de_genes(): # both DE add 20160316, remove 20170328
	"""
	cuffdiff/isoform_exp.diff.%s_vs_%s.bothde.txt
	"""
	two_group_comparison_de_transcript_de_genes()
	excel_title_info_two_group_de_transcript = """

Column A: transcript_id, the transcript id.
Column B: biotype, type of lncRNA, including antisense, ....
Column C: strand, transcription direction.
Column D: gene_id, the Ensembl gene identifier.
Column E: locus, genomic locus. chromosome: start-end.
Column F ~ G, FPKM value of each sample.
Column H: log2(fold change), fold change between two samples, differentially expressed lncRNAs were filtered by this number. inf: infinity, -inf: negative infinity.
Column I: p_value, p-value between two samples.
Column J: FDR, adjusted p-value.
Column K ~ S: LncRNA associated gene information.
Column T ~ AM: gene information, including statistics information and annotation information."""

	# write into excel 
	workbook = Workbook('BothDE(2group).xlsx')  # xlsx's name
	# the following is xlsx properties, press right key --> summary
	annotationformat =  workbook.add_format({'align': 'left', 'fg_color': '#FFFF99','font_name':'Times New Roman'})
	annotationformat.set_text_wrap()  # auto wrapping
	red = workbook.add_format({'color':'red','bold':True})
	titleformat = workbook.add_format({'bold':True,'fg_color':'red', 'align':'center','font_name':'Times New Roman'})
	titleformat2 = workbook.add_format({'bold':True, 'align':'center','font_name':'Times New Roman'})
	zhengwenformat =  workbook.add_format({'align':'left', 'font_name':'Times New Roman'}) # add 20160308
	comparisons = get_2group_comparisons()
	#print comparisons
	for comparison in comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = comparison
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		#data
		data = read_data_content_with_float(foldname+'/isoform_exp.diff.%s_vs_%s.bothde.txt' % (groupAname, groupBname))
		sheetname = '%s_vs_%s' % (groupAname, groupBname)
		worksheet = workbook.add_worksheet(sheetname[:30])  # sheetname
		colnum = len(data[0])
		worksheet.merge_range(0, 0, 19+6+info_addline, colnum-1, '', annotationformat)
		worksheet.write_rich_string(0, 0, red, groupAname+' vs '+groupBname,'\nFold Change cutoff: ', red, comparison[4], '\n','p-value cutoff: ', red, comparison[5], '\nFPKM >= ', red, str(FACTOR), ' in at least one group', excel_title_info_two_group_de_transcript, annotationformat)
		data_title = data[0]
		rownum = len(data)
		worksheet.merge_range(21+6+info_addline, 0, 21+6+info_addline, colnum-1, '%s vs %s both differentially' % (groupAname, groupBname), titleformat)
		for i in range(colnum):
			worksheet.write(0+22+6+info_addline,i,data_title[i], titleformat2)
		for row in range(1, rownum):
			for col in range(colnum):
				worksheet.write(row+22+6+info_addline, col, data[row][col],zhengwenformat)
	workbook.close()

# 20170328, add volcano plot
# volcano plot
 
"""
Volcano Plot of miRNA-seq (mature, p-value -4, fc, -5)
"""
 
def floats(da, listss):
    a= []
    for i in listss:
        a.append(float(da[i]))
    return a

def __volcano_data(inputfile, Pcutoff, FCcutoff):  # old one
    """
    volcano plot drawing
    results_mirnaseq_excel_normalization.txt.trim_groupbranch_group_negative.fc.txt    
    inputfile   ==>     contains p-value (column 4), fc (column 3)
    Pcutoff     ==>     p-value cutoff, such as 0.05
    FCcutoff    ==>     fold change cutoff, such as 2.0
    """    
    f = open(inputfile)
    header = f.readline()
    pindex = header.strip('\n').split('\t').index('p_value')
    fcindex = header.strip('\n').split('\t').index('Fold change')
    pvalues = []
    log2fcs = []
    cols = []
    ecolor=[]
    for line in f:
        linelist = line.strip('\n').split('\t')
        #print float(FCcutoff)
        try:
            fc = float(linelist[fcindex])
            #print fc
            if fc > 0:
                log2fc = numpy.log2(fc)
                log2fcs.append(log2fc)
            else:
                log2fc = -numpy.log2(-fc)
                log2fcs.append(log2fc)
            #print fc
            p=float(linelist[pindex])
            #print p
            if p<=0.000000001:
                p = 0.000000001
            elif p==1.0:
                p = 0.999999999
            pvalues.append(-numpy.log10(p))
          
            if fc >= float(FCcutoff) and p <= float(Pcutoff):
                cols.append("#FF0000")
                ecolor.append("#B30000")
            elif fc <= -float(FCcutoff) and p <= float(Pcutoff):
                cols.append('#FF0000')
                ecolor.append("#B30000")
            else:
                cols.append('#C2C2C2')
                ecolor.append("#8E8E8E")
        except:
            continue
    f.close()
    #print len(fcs)#
    return log2fcs, pvalues, cols, ecolor


def volcano_data(inputfile, Pcutoff, FCcutoff):
    """
    volcano plot drawing
    results_mirnaseq_excel_normalization.txt.trim_groupbranch_group_negative.fc.txt    
    inputfile   ==>     contains p-value (column 4), fc (column 3)
    Pcutoff     ==>     p-value cutoff, such as 0.05
    FCcutoff    ==>     fold change cutoff, such as 2.0
    """    
    f = open(inputfile)
    header = f.readline()
    pindex = header.strip('\n').split('\t').index('p_value')
    fcindex = header.strip('\n').split('\t').index('Fold change')
    newa = fcindex-2
    newb = fcindex-1
    pvalues = []
    log2fcs = []
    cols = []
    ecolor=[]
    for line in f:
        linelist = line.strip('\n').split('\t')
        newa_add1 = float(linelist[newa])+1.0
        newb_add1 = float(linelist[newb])+1.0
        if newa_add1 >= newb_add1:
            fc = newa_add1/newb_add1
            log2fc = numpy.log2(fc)
            log2fcs.append(log2fc)
        else:
            fc = -newb_add1/newa_add1
            log2fc = -numpy.log2(-fc)
            log2fcs.append(log2fc)
            #print fc
        p=float(linelist[pindex])
        #print p
        if p<=0.000000001:
            p = 0.000000001
        elif p==1.0:
            p = 0.999999999
        pvalues.append(-numpy.log10(p))
        if fc >= float(FCcutoff) and p <= float(Pcutoff):
            cols.append("#FF0000")
            ecolor.append("#B30000")
        elif fc <= -float(FCcutoff) and p <= float(Pcutoff):
            cols.append('#FF0000')
            ecolor.append("#B30000")
        else:
            cols.append('#C2C2C2')
            ecolor.append("#8E8E8E")
    f.close()
    #print len(fcs)#
    return log2fcs, pvalues, cols, ecolor 

def merge_updn_files(upfile):
	# gene_exp.diff.sicontrol-DHT_vs_sicontrol-.protein_codingtmp.up.txt
	# gene_exp.diff.sicontrol-DHT_vs_sicontrol-.protein_codingtmp.dn.txt
	f = open(upfile)
	head=f.readline()
	up=f.read()
	f.close()
	
	f2=open(upfile.replace('tmp.up.txt','tmp.dn.txt'))
	f2.readline()
	dn= f2.read()
	f2.close()
	
	fm = open(upfile.replace('tmp.up.txt','tmp.updn.txt'),'w')
	print >>fm, head,
	print >>fm, up+dn,
	fm.close()


def volcano_plot(inputfile, Pcutoff, FCcutoff, A, B):
	merge_updn_files(inputfile)
	fc, pv, co, eco= volcano_data(inputfile.replace('tmp.up.txt','tmp.updn.txt'), Pcutoff, FCcutoff)
	
	fig = figure(figsize=cm2inch(17,17))
	grid(linestyle='solid',lw=1.5, color='#C3C3C3')
	scatter(fc, pv, marker='s',linewidths=0.75,color=co, s=50, edgecolors=eco)
	#scatter(fc, pv, marker='o', color=co, s=10, edgecolors='')
	m = round(max(abs(min(fc)), max(fc)))
	#print m
	ym = max(pv)
	ymin = min(pv)
	plot([-int(m)-10, int(m)+10], [-numpy.log10(float(Pcutoff)), -numpy.log10(float(Pcutoff))], c='#63C363', lw='1.5')
	plot([-numpy.log2(FCcutoff), -numpy.log2(FCcutoff)], [ymin-3, ym+3],  c='#63C363', lw='1.5')
	plot([numpy.log2(FCcutoff),   numpy.log2(FCcutoff)], [ymin-3, ym+3],  c='#63C363', lw='1.5')
	# # gene_exp.diff.sicontrol-DHT_vs_sicontrol-.protein_codingtmp.up.txt
	xlabel('log2(Fold Change)\n%s vs %s' % (A,B))
	ylabel('-log10(pvalue)')
	xlim(-int(m)-0.5, int(m)+0.5)
	ylim(-0.2,ym+1) 
	savefig(os.path.split(inputfile)[-1]+'.png',dpi=600.1)
	print '....'
	png2tiff(os.path.split(inputfile)[-1]+'.png')

def volcano_plots():
	# no plotting for two samples
	#try:
	#	two_sample_comparisons = get_2_samples_comparisons()  # # 6W-KO-FM-2_sequence.fastq, 6W-WT-FM-1_sequence.fastq, 2.0
	#	for two_sample_comparison in two_sample_comparisons:
	#		samplea = two_sample_comparison[0].replace('_sequence.fastq', '')
	#		sampleb = two_sample_comparison[1].replace('_sequence.fastq', '')
	#		foldname2 = 'cuffdiff_'+ samplea2.split('_sequence.fastq')[0]+'_vs_'+sampleb2.split('_sequence.fastq')[0]
	#		print foldname2+'/gene_exp.diff.%s_vs_%s.protein_coding.txt' % (samplea2, sampleb2)
	#		volcano_plot(foldname2+'/gene_exp.diff.%s_vs_%s.protein_coding.txt' % (samplea2, sampleb2), two_sample_comparison[2], two_sample_comparison[3])
	#		volcano_plot(foldname2+'/gene_exp.diff.%s_vs_%s.lnc.txt' % (samplea2, sampleb2), two_sample_comparison[2], two_sample_comparison[3])
	#except:
	#	print "no two sample gene comparisons"
		
	# two groups
	two_group_comparisons = get_2group_comparisons()
	print two_group_comparisons
	for two_group_comparison in two_group_comparisons:
		groupA, groupAname, groupB, groupBname, fc, p = two_group_comparison
		foldname = 'cuffdiff_'+ groupAname+'_vs_'+groupBname
		print foldname, groupA, groupB
		# gene_exp.diff.sicontrol-DHT_vs_sicontrol-.protein_codingtmp.up.txt
		# isoform_exp.diff.sicontrol-DHT_vs_sicontrol-.lnctmp.dn.txt
		volcano_plot(foldname+'/gene_exp.diff.%s_vs_%s.protein_codingtmp.up.txt' % (groupAname, groupBname), float(p), float(fc), groupAname, groupBname)
		volcano_plot(foldname+'/isoform_exp.diff.%s_vs_%s.lnctmp.up.txt' % (groupAname, groupBname), float(p), float(fc), groupAname, groupBname)

#####20161212 add 
# add all genes, all detected genes (>0, for at least one sample), and detected genes (>0, for certain sample)(test)

def total_detected_genes(file):
	"""
	at least one sample >0, just biotype=protein_coding
	"""
	filenum = len(getfastq())
	f=open(file)
	head=f.readline()
	s=0
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[2] == 'protein_coding':
			# gene_short_name	gene_id	biotype	strand	locus	Length	Q-XL-1_FPKM	Q-XL-2_FPKM	Q-XL-3_FPKM	Q-XL-4_FPKM	Q-XL-5_FPKM	Q-XL-6_FPKM	Q-XL-1_count	Q-XL-2_count	Q-XL-3_count	Q-XL-4_count	Q-XL-5_count	Q-XL-6_count
			fpkms=[i for i in linelist[6:6+filenum] if float(i)>0]
			if len(fpkms)>0:
				s+=1
	f.close()
	return s

def total_genes(file):
	"""
	all genes, just biotype=protein_coding
	"""
	f=open(file)
	head=f.readline()
	s=0
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[2] == 'protein_coding':
			s+=1
	f.close()
	return s

def detected_genes(file):
	"""
	one sample >0
	"""
	filenum = len(getfastq())
	#filenum=6
	f=open(file)
	head=f.readline()
	headlist=head.strip().split('\t')
	d=[0]*filenum
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[2] == 'protein_coding':
			for i in range(filenum):
				if float(linelist[i+6])>0:
					d[i]+=1
	f.close()
	return [h.replace('_FPKM','') for h in headlist[6:6+filenum]], d

def detected_gene_stat(file):
	f=open(file[:-4]+'.stat.txt','w')
	samples, numbers = detected_genes(file)
	totaldetectednumber = total_detected_genes(file)
	totalgene = total_genes(file)
	print >>f, 'sample\tdetected\tall_detected(at least 1 sample >0)\tall'
	for i in range(len(samples)):
		print >>f, '%s\t%d\t%d\t%d' % (samples[i], numbers[i], totaldetectednumber, totalgene)
	f.close()
	fig = plt.figure(figsize=cm2inch(17,20))
	ax1 = fig.add_subplot(111)
	ind=range(len(samples))
	colors = ['green']*len(samples) #['#7F7F7F']*9+['#0C71C3']*9+['#AB322F']*9+['#945D8B']*9+['#6CACDE']*9+['#FF2500']*9+['#FFDF79']*9+['#545454']*9
	ax1.set_ylim(0,1)
	ax1.set_ylabel('Detected Gene Ratio',fontsize=16)
	ax2 = ax1.twinx()
	ax2.bar(ind, numbers,color=colors)
	xticks([i+0.4 for i in ind], samples,fontsize=12)
	ax2.set_xlim(0,len(samples))
	ax2.set_ylim(0,totalgene) # total genes here
	ax2.set_ylabel('Total Genes',fontsize=16)
	savefig(file[:-4]+'.png', dpi=600.1)
	png2tiff(file[:-4]+'.png')
 
#### lncRNA stat

def total_detected_lncRNAs(file):
	"""
	at least one sample >0, just biotype in lncRNAclasses
	"""
	filenum = len(getfastq())
	f=open(file)
	head=f.readline()
	s=0
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[3] in lncRNAclasses:
			# gene_short_name	gene_id	biotype	strand	locus	Length	Q-XL-1_FPKM	Q-XL-2_FPKM	Q-XL-3_FPKM	Q-XL-4_FPKM	Q-XL-5_FPKM	Q-XL-6_FPKM	Q-XL-1_count	Q-XL-2_count	Q-XL-3_count	Q-XL-4_count	Q-XL-5_count	Q-XL-6_count
			fpkms=[i for i in linelist[7:7+filenum] if float(i)>0]
			if len(fpkms)>0:
				s+=1
	f.close()
	return s

def total_lncRNAs(file):
	"""
	all genes, just just biotype in lncRNAclasses
	"""
	f=open(file)
	head=f.readline()
	s=0
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[3] in lncRNAclasses:
			s+=1
	f.close()
	return s

def detected_lncRNAs(file):
	"""
	one sample >0
	"""
	filenum = len(getfastq())
	#filenum=6
	f=open(file)
	head=f.readline()
	headlist=head.strip().split('\t')
	d=[0]*filenum
	for line in f:
		linelist=line.strip('\n').split('\t')
		if linelist[3] in lncRNAclasses:
			for i in range(filenum):
				if float(linelist[i+7])>0:
					d[i]+=1
	f.close()
	return [h.replace('_FPKM','') for h in headlist[7:7+filenum]], d

def detected_lncRNAs_stat(file): # bug, need to fix
	f=open(file[:-4]+'.stat.txt','w')
	samples, numbers = detected_lncRNAs(file)
	totaldetectednumber = total_detected_lncRNAs(file)
	totalgene = total_lncRNAs(file)
	print >>f, 'sample\tdetected\tall_detected(at least 1 sample >0)\tall'
	for i in range(len(samples)):
		print >>f, '%s\t%d\t%d\t%d' % (samples[i], numbers[i], totaldetectednumber, totalgene)
	f.close()
	fig = plt.figure(figsize=cm2inch(17,20))
	ax1 = fig.add_subplot(111)
	ind=range(len(samples))
	colors = ['green']*len(samples) #['#7F7F7F']*9+['#0C71C3']*9+['#AB322F']*9+['#945D8B']*9+['#6CACDE']*9+['#FF2500']*9+['#FFDF79']*9+['#545454']*9
	ax1.set_ylim(0,1)
	ax1.set_ylabel('Detected LncRNA Ratio',fontsize=16)
	ax2 = ax1.twinx()
	ax2.bar(ind, numbers,color=colors)
	xticks([i+0.4 for i in ind], samples,fontsize=12)
	ax2.set_xlim(0,len(samples))
	ax2.set_ylim(0,totalgene) # total genes here
	ax2.set_ylabel('Total LncRNA',fontsize=16)
	savefig(file[:-4]+'.png', dpi=600.1)
	png2tiff(file[:-4]+'.png')
 

##### clean and remove tmp

def _md5values(file):
	"""
	output the md5 string for each gzipped fastq file
	"""
	md5file=open(file, 'rb')
	md5=hashlib.md5(md5file.read()).hexdigest()
	return file, md5


def gzipfastq(file):
	os.system('cat %s |gzip >%s.gz' % (file, file))


def pgzipfastq(fastqfiles):
	threads = []
	newfastqfiles = []
	for file in fastqfiles:
		file1 = file[:-6]+'_1.fastq'
		file2 = file[:-6]+'_2.fastq'
		newfastqfiles.append(file1)
		newfastqfiles.append(file2)
	for i in newfastqfiles:
		t = threading.Thread(target=gzipfastq, args=(i,))
		threads.append(t)
	for t in threads:
		t.start()
	for t in threads:
		t.join()


def remove_cutadapt(files):
	for file in files:
		# t33_sequence_1.fastq.cutadapt_1.fastq
		try:
			os.remove(file[:-6]+'.cutadapt_1.fastq')
		except:
			print '%s already removed?!' % (file[:-6]+'.cutadapt_1.fastq')
		try:
			os.remove(file[:-6]+'.cutadapt_2.fastq')
		except:
			print '%s already removed?!' % (file[:-6]+'.cutadapt_2.fastq')
		try:
			os.remove(file[:-6]+'.cutadapt_1.fastq.fa')
		except:
			print '%s already removed?!' % (file[:-6]+'.cutadapt_1.fastq.fa')
		try:
			os.remove(file[:-6]+'.cutadapt_2.fastq.fa')
		except:
			print '%s already removed?!' % (file[:-6]+'.cutadapt_2.fastq.fa')


def md5fastq(files):
	# fastq
	f = open('md5.fastq.txt','a')
	pgzipfastq(files)
	for file in files:
		file1 = file[:-6]+'_1.fastq'
		file2 = file[:-6]+'_2.fastq'
		print >>f, '\t'.join(_md5values(file1+'.gz'))
		print >>f, '\t'.join(_md5values(file2+'.gz'))
	f.close()

def md5bam(files):
	#bam
	f2=open('md5.bam.txt','a')
	for bamfile in files:
		# T-7502_sequence.sorted.bam
		newbamfile = bamfile[:-6]+'.sorted.bam'
		print >>f2, '\t'.join(_md5values(newbamfile))
	f2.close()

def main():
	FASTQFILES = getfastq()  # ok
	getq30(FASTQFILES)       # q30
	pfastqc(FASTQFILES)      # optional # parallel
	pcutadapt(FASTQFILES)  # ok
	run_hisat2_rRNA(FASTQFILES) #OK
	fastq2fa(FASTQFILES)    # ok
	#plot_length_distribution(FASTQFILES,countthreshold, minlength, maxlength) # ok
	rRNA_mapping_rates(FASTQFILES)
	run_hisat2(FASTQFILES)
	reads_stat(FASTQFILES)
	remove_cutadapt(FASTQFILES)
	##QC over
	
	md5fastq(FASTQFILES)  # md5 value of fastq files
	psam2bam(FASTQFILES)
	md5bam(FASTQFILES)
	prun_cuffquant(FASTQFILES)  # cxb files
	 
	##################### expression (by cuffdiff) and 1 vs 1, start
	try:
		run_cuffdiff_1vs1()
	except:
		print '!!!!!!!!!!You only have one sample?'
	##################### expression (by cuffdiff) and 1 vs 1, over

	##################### group vs group, start
	try:
		run_cuffdiff_groupvsgroup()
	except:
		print '!!!!!!!!!!group vs group wrong?'
	##################### group vs group, over
 	#####################<<<<<<<<<<<<<<<<<<<<<< expression excel, start
	# done1
	
	try:
		excel_expression_all_samples_gene()
	except:
		print '!!!!!!!!!!gene expression excel wrong'
	 
	
	try:
		excel_expression_all_samples_transcript()
	except:
		print '!!!!!!!!!!transcript expression excel wrong'
	
	#################### expression excel, over>>>>>>>>>>>>>>>>>>>>>>>>
	####################<<<<<<<<<<<<<<<<<<<<<< all comparison excel, start
	# done2
	# two samples
	 
	try:
		excel_two_sample_comparison_all_gene()
	except:
		print '!!!!!!!!!!two sample all comparison for gene wrong, or without two sample comarison'
	
	
	try:
		excel_two_sample_comparison_all_transcript()
	except:
		print '!!!!!!!!!!two sample all comparison for transcript wrong, or without two sample comarison'
	
	# two groups
	
	try:
		excel_two_group_comparison_all_gene()
	except:
		print '!!!!!!!!!!two group all comparison for gene wrong, or without two group comparison'
	
	
	try:
		excel_two_group_comparison_all_transcript()
	except:
		print '!!!!!!!!!!two group all comparison for transcript wrong, or without two group comparison'
	
	#####################all comparison excel, over >>>>>>>>>>>>>>>>>>>>>
	
	
	#############
	#####################<<<<<<<<<<<<<<<<<<<<<<<de excel, start

	try:
		excel_two_sample_comparison_de_gene()
	except:
		print  '!!!!!!!!!!two sample comparison DE gene wrong, or not exist'
	
	try:
		excel_two_sample_comparison_de_transcript()
	except:
		print '!!!!!!!!!!two sample comparison DE transcript wrong, or not exist'

	try:
		excel_two_group_comparison_de_gene()
	except:
		print '!!!!!!!!!!two group comparison DE gene wrong, or not exist'
	
	try:
		excel_two_group_comparison_de_transcript()
	except:
		print '!!!!!!!!!!two group comparison DE transcript wrong, or not exist' 

	# plot
	try:
		two_samples_scatterplot_protein()
	except:
		print '!!!!!!!!!!two sample scatter protein wrong'
	
	try:
		two_samples_scatterplot_lnc()
	except:
		print '!!!!!!!!!!two sample scatter lnc wrong'
	
	try:
		two_groups_scatterplot_protein()
	except:
		print '!!!!!!!!!!two group scatter protein wrong'
	
	try:
		two_groups_scatterplot_lnc()
	except:
		print '!!!!!!!!!!two group scatter lnc wrong'
	
	try:
		volcano_plots()
	except:
		'volcano plot wrong'
 
	
	# remove the following bothde
	#try:
	#	excel_two_sample_comparison_de_transcript_de_genes()
	#except:
	#	print '!!!!!!!!!two sample both DE wrong'

	#try:
	#	excel_two_group_comparison_de_transcript_de_genes()
	#except:
	#	print '!!!!!!!!!two group both DE wrong'

	# detected genes
	try:
		detected_gene_stat('cuffdiff/Gene_Level_Profile.txt')
	except:
		print 'gene stat wrong'

	# detected LncRNAs # bug here
	#try:
	#	detected_lncRNAs_stat('cuffdiff/Transcript_Level_Profile_lnc.txt')
	#except:
	#	print 'LncRNA stat wrong'
 
	print 'Congratulations, ALL DONE'

main()

"""
# TODO:
# 1, how to calculate fpkm? use an example
# 2, images (cluster)

# BUG:
1, once a group comparison. for group comparison, 20160425 fixed
2, gene symbol not the same, e.g. ATF7 & Atf7 in mouse, so just upper, 20161212 fixed
3, two group miss at least one group >0.1 (only for p=1, since we always using 0.05, so no effect on report)
4, remove gene symbol '-' in go, pathway results!
5, gene_id not in the first part, for non model species, fixed and to test 20180404
6, 2 sample de lost one line, line: 2181, fixed in 20180723


# UPDATE
# 20181024
transform into cloud3, excel update: rich string format.

# 20180320
add go, pathway checking (remove - symbol)

# 20180211
add md5 for bam

# 20170913
fold change -- log2 fc, no add 1, just 2** log2fc

# 20170410
1, change q_value into FDR, lw
2, change png into tiff

# 20170405
change all the figures into 17cm x 17cm (or other height), by liwei

# 20170328, to test
1, add fold change into Differentially expressed lncRNA & mRNAs, and regulations
2, seperate up and down of all comparison
3, remove both DE
4, add volcano plot (add 1, then calculate FC, keep p)


# 20161212
add all genes, all detected genes (>0, for at least one sample), and detected genes (>0, for certain sample)(done)

# 20161125
add Ensembl for human if no lncRNAsource, done

# 20161116
generate read_statistic.xls directly

#20161105
if gtf with # lines, skip that line

# 20161025
1, using first-stranded in cuffquant to keep the same
#20160817
1, changed rRNA mapping into hista2, and redirect to file, and calculate rRNA mapping rate
2, hista2 mapping redirect to file, and calculate mapping rate


#20160608
1, q30 to one file

# 20160425
1, add Q30 value of each sample

# 20160421 (totally changed?)
1, change tophat2 into hisat2
2, using cuffquant
3, using cxb as the input of cuffdiff

# 20160420
1, run fastqc parallel (run every time)
2, run cutadapt parallel

# 20160316
1, combined lncRNA & mRNA (both DE), and both DE go, pw.

# 20160310
1, gene only gene level + pathwayid,name, goid, name
2, lncRNA only transcript level + associated mRNA (300kb)

# 20160308
1, excel using times new roman
# 20160307
1, add go, pathway id and name into the gene file (to test)

# 20160229
1, add gene info to gene expression, data from ncbi ftp # mouse Mus_musculus.gene_info (to test)
2, add align summary into statistics (to test)

# 20160224
1, CHANGE cutoff from q to p (default 0.05)

# 20160120
1, ADD scatter plots, both protein, lncRNA for two samples, two groups comparison

# 20160118
1, ADD lncRNA-associated mRNAs (data from pz)

# 20160113
1, parse gtf bugfixed (ok), lnc: in ['3prime_overlapping_ncrna','antisense','lincRNA','retained_intron','sense_intronic','processed_transcript','misc_RNA', 'long_noncoding']

# 20151230
1, add adaptor statistics
2, todo: add tophat statistics

# 20151228
1, all use classic-fpkm

#20151130
1, extract lncRNA

#20151028
1, add QC of fastq, adaptor, rRNA (ok)
2, kc & y: two functions

# 20151026
1, write FACTOR into config file
2, 1vs1: geometric (default); group vs group: classic-fpkm to avoid explanation problem
3, extract protein coding results for 2 group comparison
4, go, pathway for protein_coding de for 2 groups
5, transcript filter by FACTOR (ok)


20151016
1, go, pathway for protein_coding de for 2 samples


# 20151015
1, extract protein coding results for 2 sample comparison

# 20151014
1, filter excel results by FACTOR, ok

ref: 25124925, page 14:
Eight hundred and forty-five upregulated genes between sham and 7 dpa hearts were 
identified using Cuffdiff  (1.3.0), on the basis of three criteria:
(1) P (un-corrected) < 0.05; (2) fold change > 2; and (3) RPKM > 0.2 in at least one sample.

# 20151009
1, add update info
"""
