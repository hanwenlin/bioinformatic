#!/usr/bin/python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from PIL import Image
import sys
import os

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

def readdata(file):
	f = open(file)
	s=[]
	for line in f:
		if line.strip():
			s.append(line.strip())
	f.close()
	return set(s)

def png2tiff(file):
	pngfile = Image.open(file)
	# Save as TIFF
	pngfile.save(file.replace('.png','.tiff'), dpi=(600,600))
	#os.remove(file)

def _venn3(filea,fileb,filec):
	set1 = readdata(filea+'.txt')
	set2 = readdata(fileb+'.txt')
	set3 = readdata(filec+'.txt')
	fig = plt.figure(figsize=(cm2inch(17,17)), dpi=600)
	v=venn3_unweighted([set1, set2,set3 ], set_labels =(filea, fileb,filec ))
	#v.get_patch_by_id('100').set_color('red')
	#v.get_patch_by_id('100').set_color('blue')
	fig.savefig(filea+' n '+fileb+' n '+filec+'.png', dpi=600)
	png2tiff(filea+' n '+fileb+' n '+filec+'.png')
	os.remove(filea+' n '+fileb+' n '+filec+'.png')



filea='AIR1'
fileb='AIR2'
filec='AIR3'
_venn3(filea,fileb,filec)

filea='ETH1'
fileb='ETH2'
filec='ETH3'
_venn3(filea,fileb,filec)





