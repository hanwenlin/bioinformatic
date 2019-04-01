#!/usr/bin/env python
#-*- coding:utf-8 -*-
from math import sqrt

def distenceOfTwoPoints(p1,p2):
	'''
	p1,p2三维空间点，用数组表示
	'''
	return sqrt(sum([(k-v)**2 for k,v in dict(zip(p1,p2)).items()]))

