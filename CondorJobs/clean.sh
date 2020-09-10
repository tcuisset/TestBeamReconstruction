#!/usr/bin/env bash
echo "Are you sure you want to delete the contents of the out/ and log/ folders? [y/n]"
CLEANPATH="$HOME/${CMSSW_VERSION}/src/UserCode/CondorJobs/"
while true; do
    read in
    if [[ "${in}" == "y" ]]; then 
	rm -rf $CLEANPATH/out/*;
	rm -rf $CLEANPATH/log/*;
	rm $CLEANPATH/clue_*.out;
	rm $CLEANPATH/clue_*.err;
	rm $CLEANPATH/clue_*.log;
	rm $CLEANPATH/clue_*.sub;
	rm $CLEANPATH/clue*rescue*;
	rm $CLEANPATH/clue_*dag;
	rm $CLEANPATH/clue_*metrics;
	exit 0;
    elif [[ "${in}" == "n" ]]; then 
	exit 1;
    fi
done;

