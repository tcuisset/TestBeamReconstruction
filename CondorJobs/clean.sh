#!/usr/bin/env bash
echo "Are you sure you want to delete the contents of the out/, log/ and submission/ folders? [y/n]"
eval `scramv1 runtime -sh`;

CLEANPATH="$CMSSW_BASE/src/UserCode/CondorJobs"
declare -a PREFIXES=("clue" "sim_cmssw")

while true; do
    read in
    if [[ "${in}" == "y" ]]; then 
	rm -rf "$CLEANPATH/out/"*;
	rm -rf "$CLEANPATH/log/"*;
	for pref in "${PREFIXES[@]}"; do
	    rm "$CLEANPATH/${pref}"*".out";
	    rm "$CLEANPATH/${pref}"*".err";
	    rm "$CLEANPATH/${pref}"*".log";
	    if [ $pref == "clue" ]; then
		rm "$CLEANPATH/${pref}"*".sub";
	    fi
	    rm "$CLEANPATH/${pref}"*"sub.condor.sub";
	    rm "$CLEANPATH/${pref}"*"rescue"*;
	    rm "$CLEANPATH/${pref}"*"dag";
	    rm "$CLEANPATH/${pref}"*"metrics";
	done;
	rm "$CLEANPATH/submission/selection/"*"sub";
	rm "$CLEANPATH/submission/analysis/"*"sub";
	exit 0;
    elif [[ "${in}" == "n" ]]; then 
	exit 1;
    fi
done;

