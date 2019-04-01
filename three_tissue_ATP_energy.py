#!/usr/bin/env python
#-*- coding:utf8 -*-
import math


while True:
	tissue_cell = input("please choose one tissue('muscle':'m','liver':'l','brain':'b'): ")
	if tissue_cell in ['m','l','b']:
		if tissue_cell=='m':
			ATP,ADP,Pi = 8.0,0.9,8.0
		elif tissue_cell=='l':
			ATP,ADP,Pi = 3.5,1.8,5.0
		else:
			ATP,ADP,Pi = 2.6,0.7,2.7
		R,T = 0.00831,298
		deltaGO = -30.5
		energy = deltaGO+R*T*math.log(ADP*Pi/ATP)
		print(energy)
	else:
		break
