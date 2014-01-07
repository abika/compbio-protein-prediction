#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Dec 13 16:21:09 2013

@author: Alexander Bikadorov

TMalign AND blastp need to be in PATH!

"""

import sys
import logging
import argparse
import os
import subprocess

import _myutils

NUM_MODELS = 20 # number of used models per decoy

def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("target_pdb_file", type=str, help="Target pdb file")
    parser.add_argument("pdb_dir", type=str, help="Path to structure database (directory with pdb files)")
    parser.add_argument("fasta_dir", type=str, help="Path to sequence database (directory with fasta files)")
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
    field_dict_list.sort(key=lambda d: float(d["score"])) 
    
    print('scores:')
    print('score gdtmm_full description')
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

#BLAST_OUTP = '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
BLAST_OUTP = ["10 bitscore score"]
def _run_blast(target_fasta_file, model_fasta_file):
    args = ['blastp', '-query', target_fasta_file, '-subject', model_fasta_file, '-outfmt'] + BLAST_OUTP
    out = subprocess.check_output(args).strip().decode("utf-8")
    return float(out.split(',')[0]) if len(out) else 0

def main(argv=sys.argv):
    logging.getLogger().setLevel(logging.INFO)
    
    args = _arguments()
    target_pdb_path = os.path.abspath(args.target_pdb_file)
    target_base_dir = os.path.join(*_myutils.split_path(target_pdb_path)[:-2])
    fasta_dir = os.path.abspath(args.fasta_dir)
    
    # step 1: run compbio_app
    # TODO: currently done manually
    #compbio_app.linuxgccdebug @flags -database "$ROSETTA_DIR"/rosetta_database --nstruct 10 
    
    # step 2: read scorefile 
    field_dict_list = _read_scorefile(target_base_dir)
    
    # step 3: get structure scores with tmalign
    decoy_dict = field_dict_list[0] # TODO: use more than 1 decoy
    print('decoy selected: '+decoy_dict['description'])
    
    decoy_pdb_file = os.path.join(target_base_dir, 'models', decoy_dict['description'] + '.pdb')
    struc_scores = _tmalign(decoy_pdb_file, args.pdb_dir)

    output_file_str = os.path.splitext(decoy_pdb_file)[0]+'_tmalign_scores.txt'
    _myutils.write_file(output_file_str , '\n'.join(' '.join(t) for t in struc_scores))    

    struc_scores = struc_scores[:NUM_MODELS]
    
    # step 4: get sequence scores
    target_fasta_file = os.path.splitext(target_pdb_path)[0]+'.fasta'
    scores = []
    for model_pdb_file, struc_score in struc_scores:
        model_fasta_file = os.path.join(fasta_dir, model_pdb_file.rstrip('.pdb') + '.fasta')
        seq_score = _run_blast(target_fasta_file, model_fasta_file)
        scores.append((model_pdb_file, struc_score, seq_score))
    
    print('scores :')
    print('\n'.join(str(t) for t in scores))    
    
    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())

