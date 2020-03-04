#!/usr/bin/env bash

#Clean jobs output files
echo "Are you sure you want to delete the contents of the out/ and log/ folders? [y/n]"
while true; do
    read in
    if [[ "${in}" == "y" ]]; then 
	rm -rf out/*;
	rm -rf log/*;
	break;
    elif [[ "${in}" == "n" ]]; then 
	break;
    else
	echo "Please provide a valid answer. [y/n]"
    fi
done;

#Generate file with the IDs of the Ntuples to be processed
#awk #1: removes everything up to 'ntuples_'
#awk #2: removes everything up to after the number
#awk #3: removes all blank lines
echo "Generating new 'ntuple_ids.txt' file..."
ls -l ../ntuples/ | awk '{print substr($9,8)}' | awk '{print substr($1,1,length($1)-5)}' | awk '/./' > ntuple_ids.txt
