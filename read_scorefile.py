#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Fri Dec 13 16:01:42 2013

@author: Alexander Bikadorov

Read scorefile and output sorted list

"""

import sys
import argparse


def _arguments():
    parser = argparse.ArgumentParser()
    #Add options one at a time
    parser.add_argument("score_file", type=str, help="Name of the score file")
    #parser.add_option("--native", type="string", dest="native", help="Name of the native .pdb file")

    #Parse args and get options
    return parser.parse_args()


def main(argv=sys.argv):
    
    args = _arguments()

    # load score file input
    score_file = open(args.score_file)
    lines_list = score_file.readlines()
    score_file.close()
    
    # create list of lists with all fields
    field_list = [str_.split() for str_ in lines_list]
    # filter irrelevant lines
    field_list = [l for l in field_list if l[0] == "SCORE:"] 
    # create list of dictionaries (with keys from first line in file)
    field_dict_list = [dict(zip(field_list[0], l)) for l in field_list[1:]]
    
    # sort list
    field_dict_list.sort(key=lambda d: d["score"])    
    
    for d in field_dict_list:
        print(d["score"], d["gdtmm_full"], d["description"])
    
    print("DONE!")

if __name__ == "__main__":
    sys.exit(main())

