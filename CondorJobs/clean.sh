#!/usr/bin/env bash
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
