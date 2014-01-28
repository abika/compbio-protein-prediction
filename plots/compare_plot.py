#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys
import matplotlib.pyplot as plt
import os

from pylab import *
from AnnoteFinder import *
from PymolLauncher import *
from pymol import *

def _read_scorefile(score_file):
    # load score file input
    cont_str = open(score_file, 'r').read()
    lines_list = [url_str.rstrip() for url_str in cont_str.splitlines()]
    # create list of lists with all fields
    field_list = [str_.split() for str_ in lines_list]
    # filter irrelevant lines
    field_list = [l for l in field_list if l[0] == "SCORE:"] 
    # create list of dictionaries (with keys from first line in file)
    field_dict_list = [dict(zip(field_list[0], l)) for l in field_list[1:]]
    
    # sort list
    field_dict_list.sort(key=lambda d: float(d["score"])) # smaller score better
    '''
    print('decoy scores:')
    print('score gdtmm_full description')
    for d in field_dict_list:
        print(d["score"], d["gdtmm_full"], d["description"])
    '''    
    return field_dict_list

def plot_scores(scoreFile,label):
    field_dict_list = _read_scorefile(scoreFile)
    
    # plot gdt vs. energy
    energy_gdt_list = [(float(d["score"]), float(d["gdtmm_full"])) for d in field_dict_list]
    energy_list, gdt_list = zip(*energy_gdt_list)
    plt.plot(gdt_list, energy_list,'o',label=label)

def helpMe():
    print("usage:")
    print("%s <.fsc file> <label> ..."%(sys.argv[0]))

def main(argv=sys.argv):
    # logging.getLogger().setLevel(logging.INFO)
    
    if(argv[1] == "-h" or argv[1] == "--help"):
        helpMe()
        exit(0)
        
    plt.xlabel("GDT_mm")
    plt.ylabel("energy")
    plt.axis(xmin=0, xmax=1)
    for i in range(1,len(argv),2):
        a = argv[i]
        b = argv[i+1]
        plot_scores(a,b)
    
    #field_dict_list = _read_scorefile(score_file_path)
    #plt.plot(gdt_list[5:], energy_list[5:], 'ko') # else in black
    #plt.title("title")
    
    plt.legend()
    plt.show()
    
main()