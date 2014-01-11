#!/bin/bash

# $1 - number of processes
# $2 - pdb file of the protein of interest
# $3 - string that should be placed after every execution of TMalign
# $n, n>3 - pdb file of proteins to which we should compare $2

separator=$3

# runs TMalign for between a given protein and the rest of the arguments, writing $separator to the end of each output
f()
{
	p=$1
	shift
	for i in $@
	do
		TMalign $p $i
		echo "$separator"
		echo "$i"
		echo "$separator"
	done
}

processes=$1
protein=$2
shift 3
echo "Processes: \"$processes\""
echo "Protein: \"$protein\""
echo "Separator: \"$separator\""
part=$(($# / processes))
echo $#
echo "part: \"$part\""
# Distributes the input into equally-sized parts, each being executed sequentially
for (( i=0; i < $processes - 1; i++))
do
	l=${@:1:part}
	f "$protein" "$l" &
	shift $part
done

f "$protein" "$@"

wait