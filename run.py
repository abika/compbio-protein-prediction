#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Dec 13 16:21:09 2013

@author: Alexander Bikadorov

TMalign needs to be in PATH!

"""

import sys
import logging
import argparse
import os
import subprocess

import _myutils

def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_pdb_file", type=str, help="Target pdb file")
    parser.add_argument("pdb_dir", type=str, help="Path to pdb database (directory with pdb files)")
    return parser.parse_args()

def _read_scorefile(target_base_dir):
    # load score file input
    score_file_str = os.path.join(target_base_dir, 'models', 'score.fsc')
    lines_list = _myutils.read_file_lines(score_file_str)
    
    # create list of lists with all fields
    field_list = [str_.split() for str_ in lines_list]
    # filter irrelevant lines
    field_list = [l for l in field_list if l[0] == "SCORE:"] 
    # create list of dictionaries (with keys from first line in file)
    field_dict_list = [dict(zip(field_list[0], l)) for l in field_list[1:]]
    
    # sort list
    field_dict_list.sort(key=lambda d: d["score"]) 
    
    print('scores:')
    for d in field_dict_list:
        print(d["score"], d["gdtmm_full"], d["description"])
        
    return field_dict_list

def _run_tmalign(target_pdb, template_pdb):
    args = ["TMalign "+target_pdb+" "+template_pdb+" -a | grep ^TM-score "
            """| tail -n 1 | sed -n 's|^.*TM-score= \([0-9.]\+\).*|\\1|p'"""]
    out = subprocess.check_output(args, shell=True).rstrip().decode("utf-8")
    return os.path.basename(template_pdb), out

def _tmalign(model_pdb_file, pdb_dir):
    print('running TMAlign for '+model_pdb_file+'...')    
    scores = []
    CHUNK_SIZE = 100
    files = _myutils.files_in_dir(pdb_dir, '*.pdb')
    # run tmalign for each template in database
    for chunk in _myutils.group_it(files, CHUNK_SIZE):
        scores += [_run_tmalign(model_pdb_file, template) for template in chunk]
        print("%d/%d" % (len(scores), len(files)))
    
    scores.sort(key=lambda t: t[1], reverse=True)
    
    return scores
    
    
def main(argv=sys.argv):
    logging.getLogger().setLevel(logging.INFO)
    
    args = _arguments()
    target_base_dir = os.path.join(*_myutils.split_path(os.path.abspath(args.target_pdb_file))[:-2])
    
    # step 1: run compbio_app
    # TODO: currently done manually
    
    # step 2: read scorefile 
    field_dict_list = _read_scorefile(target_base_dir)
    
    # step 3: run tmalign
    model_pdb_file = os.path.join(target_base_dir, 'models', field_dict_list[0]["description"] + '.pdb')
    scores = _tmalign(model_pdb_file, args.pdb_dir)
    
    print('\n'.join(str(t) for t in scores[:10]))

    output_file_str = os.path.splitext(model_pdb_file)[0]+'_tmalign_scores.txt'
    _myutils.write_file(output_file_str , '\n'.join(' '.join(t) for t in scores))    
    
    # step 4: get sequence scores
    # TODO
    
    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())

