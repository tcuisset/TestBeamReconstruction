#!/usr/bin/env bash
declare -a TAGS=( "W2p3_dpos1p3"
		  "W2p9_dpos1p3"
		  "W4p0_dpos1p3"
		  "W5p0_dpos1p3"
		  "W2p9_dpos3p4"
		  "W2p9_dpos999" )
for t in "${TAGS[@]}"; do
    bash DataProcessing/join_ntuples.sh --datatype data --showertype em --analysistype clusterdep --tag ${t} &
done
