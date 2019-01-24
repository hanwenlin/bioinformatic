#!/usr/bin/python
# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2, venn2_circles
from PIL import Image
import sys

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


def _venn2(filea,fileb):
	set1 = readdata(filea+'.txt')
	set2 = readdata(fileb+'.txt')
	fig = plt.figure(figsize=(cm2inch(17,17)), dpi=600)
	venn2([set1, set2 ], (filea, fileb ))
	fig.savefig(filea+' n '+fileb+'.png', dpi=600)
	png2tiff(filea+' n '+fileb+'.png')

#filea=sys.argv[1]
#fileb=sys.argv[2]

filea='N'
fileb='T'
_venn2(filea,fileb)
