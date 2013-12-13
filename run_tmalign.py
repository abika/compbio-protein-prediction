#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Dec 13 16:21:09 2013

@author: Alexander Bikadorov

Run TMalign for all templates in a directory against one target and...
TODO: do nothing

TMalign needs to be in PATH

"""

import sys
import argparse
import os
import glob
import subprocess 

def _files_in_dir(dir_, regex='*.*'):
    abs_path = os.path.abspath(dir_)
    return glob.glob(os.path.join(abs_path, regex))

def _arguments():
    parser = argparse.ArgumentParser()
    #Add options one at a time
    parser.add_argument("target_pdb_file", type=str, help="Target pdb file")
    parser.add_argument("pdb_dir", type=str, help="Path to pdb database (directory with pdb files)")

    #Parse args and get options
    return parser.parse_args()

def _run_tmalign(target_pdb, template_pdb):
    #args = ['TMalign', target_pdb, template_pdb, '-a', '|', 'grep ^TM-score']
    args = ['TMalign '+target_pdb+" "+template_pdb+' -a | grep ^TM-score | tail -n 1']
    out = subprocess.check_output(args, shell=True)
    print(out)
    return out
    
def main(argv=sys.argv):
    
    args = _arguments()
    
    scores = [_run_tmalign(args.target_pdb_file, template) for template in _files_in_dir(args.pdb_dir, '*.pdb')]
    
    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())

