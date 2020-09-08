#!/usr/bin/env bash
echo "Are you sure you want to delete the contents of the out/ and log/ folders? [y/n]"
while true; do
    read in
    if [[ "${in}" == "y" ]]; then 
	rm -rf out/*;
	rm -rf log/*;
	rm clue_*.out;
	rm clue_*.err;
	rm clue_*.log;
	rm clue_*.sub;
	rm clue*rescue*;
	rm clue_*dag;
	rm clue_*metrics;
	exit 0;
    elif [[ "${in}" == "n" ]]; then 
	exit 1;
    fi
done;

