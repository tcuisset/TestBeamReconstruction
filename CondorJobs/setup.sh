#!/usr/bin/env bash

#Generate file with the IDs of the Ntuples to be processed
#awk #1: removes everything up to 'ntuples_'
#awk #2: removes everything up to after the number
#awk #3: removes all blank lines
MAINPATH=$HOME/$CMSSW_VERSION"/src/TestBeamReconstruction/"
echo "Generating new 'ntuple_ids.txt' file..."
ls -l $MAINPATH/ntuples_in/ | awk '{print substr($9,8)}' | awk '{print substr($1,1,length($1)-5)}' | awk '/./' > $MAINPATH"/CondorJobs/"ntuple_ids.txt
