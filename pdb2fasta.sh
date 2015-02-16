#! /bin/bash

# 2012 - P. Poulain

usage="Usage: $0 file.pdb [[file2.pdb] ...]"

#=============================================================================
# input data
#=============================================================================
# check number of arguments
if [ ! $# -ge 1 ]
then
    echo "Argument error" 1>&2
    echo $usage 1>&2
    exit 1
fi

#=============================================================================
# functions
#=============================================================================
# list chains in PDB
list_chain() {
    awk '/^ATOM/ && $3 == "CA" {print substr($5,1,1)}' | uniq # If the number of residue is greater than 999, 
    # the chain number increase with the new residue number , for example: A1000,A1002,A1003, it would be a new chains.
    }

# extract residue sequence from ATOM lines
# take residue from CA atom since
# there is one CA atom per residue
extract_seq_chain() {
    awk -v ch=$1 '/^ATOM/ && $3 == "CA" && $5 == ch {print $4}'
}

extract_seq() {
    awk '/^ATOM/ && $3 == "CA" {print $4}'
}

# convert newline by space
# for sed lowers this works too: sed ':a;N;$!ba;s/\n/ /g'
remove_newline() {
    tr '\n' ' ' 
}

# convert 3-letter residue code to 1-letter
convert_aa() {
    sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g'
}

# remove space between residues
remove_space() {
    sed 's/ //g' 
}

# split fasta sequence at 60 characters (easier to read than 80)
split_60() {
    fold -w 60
}

split_80() {
    fold -w 80
}

#=============================================================================
# list chain in PDB file
#=============================================================================
#chains=$(cat $name | list_chain)

#=============================================================================
# try to extract a sequence
#=============================================================================
#for chain in $chains
#do

#name=$1

for name in "$@"
do  
    # check first argument is an existing regular file
    if [ ! -f $name ]
    then
        echo "$name is not a regular file" 1>&2
        echo $usage 1>&2
        exit 1
    fi
		
		destination=$(echo "$name" | sed s/\\..*/\\.fasta/)
		echo $destination

    #sequence=$(cat $name | extract_seq_chain  $chain | remove_newline | convert_aa | remove_space)
    sequence=$(cat $name | extract_seq | remove_newline | convert_aa | remove_space)
    if [[ -n $sequence ]]
    then
        size=$(echo $sequence | wc -c)
        size=$((size-1))
        #echo ">${name%.pdb} | chain $chain | $size aa"
        echo ">${name%.pdb}|PDBID|CHAIN|SEQUENCE" > "$destination"
        #echo $sequence | split_80
        echo $sequence >> "$destination"
    fi
done
