#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Dec 13 16:21:09 2013

@author: Alexander Bikadorov, Marcelo Millani

Dependencies:
- TMalign needs to be in PATH
- blastp needs to be in PATH
- grep needs to be in PATH
- head needs to be in PATH

"""

import sys
import logging
import argparse
import os
import subprocess
import itertools

import _myutils

NUM_STRUC_MODELS = 20 # number of used models per decoy
NUM_CONSTR_MODELS = 3 # number of models used for constraint extraction per iteration
WEIGHT_STRUC = 1
WEIGHT_SEQ = 0.01
D = 1 # lower bound distance in angstrom for constraint extraction
SD = '0.2' # scaling value for constraint distance
TMALIGN_SEPARATOR = "# SEPARATOR #"
PROCESSES = 8

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
    field_dict_list.sort(key=lambda d: float(d["score"])) # smaller score better
    
    print('scores:')
    print('score gdtmm_full description')
    for d in field_dict_list:
        print(d["score"], d["gdtmm_full"], d["description"])
        
    return field_dict_list

def _run_para_tmalign(target_pdb,templates_pdbs):
    # parser = r"grep ^TM-score | grep [0-9]*\\.[0-9]* -o | head -n 1"
    args = ["paratmalign %d %s %s %s"%(PROCESSES,target_pdb, TMALIGN_SEPARATOR ," ".join(templates_pdbs))]
    all = subprocess.check_output(args, shell=True).rstrip().decode("utf-8")
    executions = all.split(TMALIGN_SEPARATOR)
    name = os.path.basename(template_pdb)
    scores = []
    # gets the score of each execution
    for e in executions:
      score = subprocess.check_output("echo $s | $s"$(e,parser), shell=True).rstrip().decode("utf-8")
      scores.append((name,float(score)))
      
    return scores


def _run_tmalign(target_pdb, template_pdb):
    # TODO: not sure if the correct score is extracted, '-a' maybe wrong too
    # args = ["TMalign "+target_pdb+" "+template_pdb+" -a | grep ^TM-score "
    #        """| tail -n 1 | sed -n 's|^.*TM-score= \([0-9.]\+\).*|\\1|p'"""]
    # TMalign says:
    # There are two TM-scores reported. You should use the one normalized by
    # the length of the protein you are interested in.
    # I think we should stick to that
    parser = r"grep ^TM-score | grep [0-9]*\\.[0-9]* -o | head -n 1"
    args = ["TMalign %s %s | "%(target_pdb, template_pdb) + parser]
    score = subprocess.check_output(args, shell=True).rstrip().decode("utf-8")
    return os.path.basename(template_pdb), float(score)

def _para_tmalign(model_pdb_file, pdb_dir):
    print('running paraTMAlign for '+model_pdb_file+'...')    
    scores = []
    CHUNK_SIZE = 100
    files = _myutils.files_in_dir(pdb_dir, '*.pdb')
    # run tmalign for each template in database
    for chunk in _myutils.group_it(files, CHUNK_SIZE):
        scores += _run_para_tmalign(model_pdb_file, chunk)
        print("%d/%d" % (len(scores), len(files)))
    
    scores.sort(key=lambda t: t[1], reverse=True) # higher score better
    return scores

def _tmalign(model_pdb_file, pdb_dir):
    print('running TMAlign for '+model_pdb_file+'...')    
    scores = []
    CHUNK_SIZE = 100
    files = _myutils.files_in_dir(pdb_dir, '*.pdb')
    # run tmalign for each template in database
    for chunk in _myutils.group_it(files, CHUNK_SIZE):
        scores += [_run_tmalign(model_pdb_file, template) for template in chunk]
        print("%d/%d" % (len(scores), len(files)))
    
    scores.sort(key=lambda t: t[1], reverse=True) # higher score better
    return scores

#BLAST_OUTP = '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'
BLAST_OUTP = ["10 bitscore score"]
def _run_blast(target_fasta_file, model_fasta_file):
    args = ['blastp', '-query', target_fasta_file, '-subject', model_fasta_file, '-outfmt'] + BLAST_OUTP
    out = subprocess.check_output(args).strip().decode("utf-8")
    return float(out.split(',')[0]) if len(out) else 0

def _residue_indices(res_list, bound_list):
    it = itertools.count(1)
    return [a for a, x in zip((next(it) if a != '-' else None for a in res_list), bound_list) if x == ':']
    
def _run_tmalign_constr(target_pdb_file, template_pdb_file, d):
    """Run TMalign to get distance constraints between atom pairs.
       A list of tuples is returned containing the indices for amino acid pairs.
       d: minimum distance in angstrom.
    """
    args = ["TMalign "+target_pdb_file+" "+template_pdb_file+" -d "+str(d)+" | tail -n 4"]
    outp_lines = subprocess.check_output(args, shell=True).rstrip().decode("utf-8").splitlines()
    return list(zip(_residue_indices(outp_lines[0], outp_lines[1]), _residue_indices(outp_lines[2], outp_lines[1])))

def main(argv=sys.argv):
    logging.getLogger().setLevel(logging.INFO)
    
    args = _arguments()
    target_pdb_path = os.path.abspath(args.target_pdb_file)
    target_base_dir = os.path.join(*_myutils.split_path(target_pdb_path)[:-2])
    pdb_dir = os.path.abspath(args.pdb_dir)
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
    struc_scores = _tmalign(decoy_pdb_file, pdb_dir)

    output_file_str = os.path.splitext(decoy_pdb_file)[0]+'_tmalign_scores.txt'
    _myutils.write_file(output_file_str , '\n'.join(' '.join(str(t)) for t in struc_scores))    

    struc_scores = struc_scores[:NUM_STRUC_MODELS]
    
    # step 4: get sequence scores
    # TODO: between target and model, or between decoy and model?
    target_fasta_file = os.path.splitext(target_pdb_path)[0]+'.fasta'
    scores = []
    for model_pdb_file, struc_score in struc_scores:
        model_fasta_file = os.path.join(fasta_dir, model_pdb_file.rstrip('.pdb') + '.fasta')
        seq_score = _run_blast(target_fasta_file, model_fasta_file)
        scores.append((model_pdb_file, struc_score, seq_score, WEIGHT_STRUC * struc_score + WEIGHT_SEQ * seq_score))
    
    scores.sort(key=lambda t: t[3], reverse=True)
    print('scores :')
    print('\n'.join(str(t) for t in scores))    

    # step 5: get constraints from selected models    
    constr_models = scores[:1] # TODO: how many models?
    res_pairs = []
    for model_file, strs, seqs, fs in constr_models:
        model_pdb_path = os.path.join(pdb_dir, model_file)
        res_pairs += _run_tmalign_constr(target_pdb_path, model_pdb_path, D)
    
    # remove inconsistent pairs
    res_pairs = _myutils.remove_dups(res_pairs, comp_item_index=0)
    res_pairs = _myutils.remove_dups(res_pairs, comp_item_index=1)
    
    # save constraint file
    
    # biopython not required, the residue numbering should be the same for tmalign and rosetta
    #parser = Bio.PDB.PDBParser(PERMISSIVE=1)
    #structure = parser.get_structure('id', target_pdb_file)

    # AtomPair: Atom1_Name Atom1_ResNum Atom2_Name Atom2_ResNum Func_Type Func_Def
    # AtomPair SG 5 V1 32 HARMONIC 0.0 0.2
    ros_constr = ['AtomPair CA '+str(n1)+' CA '+str(n2)+' HARMONIC 0.0 '+SD for n1, n2 in res_pairs]
    constr_file_path = os.path.join(target_base_dir,'inputs', 'ros_constraints.txt')
    _myutils.write_file(constr_file_path, '\n'.join(ros_constr) + '\n')
    
    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())
