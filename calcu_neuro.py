#!/usr/bin/env python
#-*- coding:utf-8 -*-

data =[]
with open('neuro_data.txt','r') as f:
	for i in f:
		data.append(float(i.strip()))
	
mins,maxs = min(data),max(data)
average = sum(data)/len(data)
std = sum([(i-average)**2 for i in data])/len(data)

with open('calculation_neuro.txt','w') as f:
	f.write('the minimum of the data is %.2f\n' %mins)
	f.write('the maximum of the data is %.2f\n' %maxs)
	f.write('the average is %.2f\n' %average)
	f.write('the std is %.2f\n' %std)
	
	
