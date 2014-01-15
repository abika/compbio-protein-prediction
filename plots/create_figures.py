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

def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("score_file", type=str, help="Rosetta score file (in model directory)")
    parser.add_argument("native_pdb", type=str, help="Native pdb file")
    return parser.parse_args()

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
    
    print('decoy scores:')
    print('score gdtmm_full description')
    for d in field_dict_list:
        print(d["score"], d["gdtmm_full"], d["description"])
        
    return field_dict_list

def main(argv=sys.argv):
    logging.getLogger().setLevel(logging.INFO)
    
    args = _arguments()
    score_file_path = os.path.abspath(args.score_file)
    model_dir = os.path.dirname(score_file_path)
    target_dir, base_dir_name = os.path.split(os.path.abspath(os.path.join(model_dir, '..')))
    target_dir_name = os.path.basename(target_dir)

    field_dict_list = _read_scorefile(score_file_path)
    
    # plot gdt vs. energy
    field_dict_list.sort(key=lambda d: float(d["score"])) # smaller score better
    energy_gdt_list = [(float(d["score"]), float(d["gdtmm_full"])) for d in field_dict_list]
    energy_list, gdt_list = zip(*energy_gdt_list)
    plt.plot(gdt_list[:5], energy_list[:5], 'ro') # 5 best in red
    plt.plot(gdt_list[5:], energy_list[5:], 'ko') # else in black
    #plt.title("title")
    plt.xlabel("GDT_mm")
    plt.ylabel("energy")
    plt.axis(xmin=0, xmax=1)         
    
    fname = os.path.join(os.getcwd(), target_dir_name+'_'+base_dir_name+'_energy_gdt.png')
    plt.savefig(fname)
    #plt.show()
    
    # args for PymolLauncher
    gdts, scores, pdbs = gdt_list, energy_list, [d["description"] for d in field_dict_list]

    """ Code that links the data to the PymolLauncher Object
        Uncomment this if you have edited this file and are ready 
        to link the PymolLauncher to the scatter plot"""
    pl1 = PymolLauncher(gdts, scores, pdbs)
    pl1.set_native(args.native_pdb)
    pdb_dir = model_dir
    print('pdb_dir: '+pdb_dir)    
    pl1.set_pdb_dir(pdb_dir)
    connect('button_press_event', pl1)
    gca().set_autoscale_on(False)
    
    import pymol
    pymol.finish_launching()
    #from pymol import cmd
    pl1.select_best()
    show()

    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())

