#!/usr/bin/env bash

#Clean jobs output files
echo "Are you sure you want to delete the contents of the out/ and log/ folders? [y/n]"
while true; do
    read in
    if [[ "${in}" == "y" ]]; then 
	rm -rf out/*;
	rm -rf log/*;
	exit 0;
    elif [[ "${in}" == "n" ]]; then 
	exit 0;
    else
	echo "Please provide a valid answer. [y/n]"
    fi
done;

#Generate file with the IDs of the Ntuples to be processed
ls -l ../ntuples/ | awk '{print substr($9,8)}' | awk '{print substr($1,1,length($1)-5)}' > ntuple_ids.txt
