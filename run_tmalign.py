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
import logging

def _files_in_dir(dir_, regex='*.*'):
    """Return a list of all files in 'dir_' matching 'regex' with absolute path"""
    abs_path = os.path.abspath(dir_)
    if not os.path.isdir(abs_path):
        logging.warning('does not exist/is not a directory: '+abs_path)
    return glob.glob(os.path.join(abs_path, regex))

def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_pdb_file", type=str, help="Target pdb file")
    parser.add_argument("pdb_dir", type=str, help="Path to pdb database (directory with pdb files)")
    return parser.parse_args()

def _run_tmalign(target_pdb, template_pdb):
    args = ["TMalign "+target_pdb+" "+template_pdb+" -a | grep ^TM-score "
            """| tail -n 1 | sed -n 's|^.*TM-score= \([0-9.]\+\).*|\\1|p'"""]
    out = float(subprocess.check_output(args, shell=True))
    #print(os.path.basename(template_pdb), out)
    return out

def _group_it(l, n):
    """Yield successive 'n'-sized chunks from sequence 'l'.
       The last chunk contains len('l') modulo 'n' elements.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

def main(argv=sys.argv):
    
    args = _arguments()
    
    scores = []
    CHUNK_SIZE = 100
    files = _files_in_dir(args.pdb_dir, '*.pdb')
    for chunk in _group_it(files, CHUNK_SIZE):
        scores += [_run_tmalign(args.target_pdb_file, template) for template in chunk]
        print("%d/%d" % (len(scores), len(files)))
    
    scores.sort(key=lambda t: t[1])
    
    print(scores)
    
    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())

