#!/usr/bin/env python
import sys
#import os
from math import sqrt
from optparse import OptionParser


def median(x):
	l = len(x)
	if l % 2 == 0:
		return (x[l/2] + x[l/2 +1])/2
	else:
		return x[l/2]

def avg(x):
	return sum(x)/len(x)
	
def variance(x):
	a = avg(x)
	return avg(map(lambda y: (y-a)**2, x))
	
def stddev(x):
	return sqrt(variance(x))

def print_statistics(x):
	# print("Length: %d"%(len(x)))
	print("Min: %f"%(min(x)))
	print("Median: %f"%(median(x)))
	print("Average: %f"%(avg(x)))
	print("Standard Deviation: %f"%(stddev(x)))
	print("Max: %f"%(max(x)))


def compareGDT(a,b):
	return compareEnergy(a,b,1)

def compareEnergy(a,b,i=0):
	x = a[i]
	y = b[i]
	if x < y:
		return -1
	if x == y:
		return 0
	else:
		return 1

#Add options for the script here!
def add_options( parser ):
	#Add options one at a time
	parser.add_option("--score_file", type="string", dest="score_file", help="Name of the score file")
	# parser.add_option("--native", type="string", dest="native", help="Name of the native .pdb file")

	#Parse args and get options
	options, args = parser.parse_args()
	return options, args


def main():
	# create a OptionParser object
	parser = OptionParser()
	options, args = add_options(parser)

	# load score file input
	score_file = open(options.score_file)
	lines_list = score_file.readlines()
	score_file.close()
	field_list = [str_.split() for str_ in lines_list]
	# build a list of dictionaries (with keys from first line in file)
	field_list = [l for l in field_list if l[0] == "SCORE:"]
	
	field_dict_list = [dict(zip(field_list[0], l)) for l in field_list[1:]]
	
	# prints statistics
	gdt_list = [float(d["gdtmm_full"]) for d in field_dict_list]
	energy_list = [float(d["score"]) for d in field_dict_list]
	
	energy_list.sort()
	gdt_list.sort()
	print("Energy")
	print_statistics(energy_list)
	print("######\nGDT")
	print_statistics(gdt_list)


if __name__ == '__main__':
	sys.exit(main())
else:
	print("Loaded as a module!")   
