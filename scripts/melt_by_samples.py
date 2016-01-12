#!/usr/bin/python

import sys
import csv

def melt_by():
	try:
		return { 'samples': -3, 'HET': -2, 'HOM_ALT': -1 }[sys.argv[1]]
	except:
		return -3

lines = sys.stdin.readlines()
print '\t'.join(['sample']+lines.pop(0).split('\t')[:-3])
for line in lines:
	cols = line.split('\t')
	for sample in map(lambda s: s.strip(), cols[melt_by()].split(',')):
		if sample == '': continue
		print '\t'.join([sample]+cols[:-3])
