#!/usr/bin/env bash

varExists() { 
    # Checks whether a certain environment variable already exists
    # Arguments:
    # 1. Variable being checked
    local flag=false;
    if [ -z "${1}" ]; then
	flag=true;
    fi
    echo $flag;
}

export XRD_NETWORKSTACK=IPv4
export CMSSWVER="CMSSW_11_1_0_pre2"
export SCRAM_ARCH="slc7_amd64_gcc820"

if [ $(varExists "${INIT_FOLDER}") = true ] && [ $(varExists "${CMSSW_PATH}") = true ] &&
    [ $(varExists "${HOME_DIR}") = true ] && [ $(varExists "${FULL_PATH}") = true ]; then
    INIT_FOLDER=$(pwd);
    CMSSW_PATH="/${CMSSWVER}/src/";
    HOME_DIR="/afs/cern.ch/user/b/bfontana";
    FULL_PATH="${HOME_DIR}""${CMSSW_PATH}";
else
    echo "Use different variable names.";
    exit 0;
fi

cd "${FULL_PATH}";

source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh` #cmsenv substitute

#back to the job folder
cd "${INIT_FOLDER}";

clue "${1}" "${2}" "${3}" "${4}" "${5}" 0 1 1

outfolder="/afs/cern.ch/work/b/bfontana/TestBeamReconstruction/"
outfolder_move="/eos/user/b/bfontana/TestBeamReconstruction/"
outfile="${1}_${2}_${3}_${4}_${5}.csv"
if [ -r "${outfolder}${outfile}" ]; then
    mv "${outfolder}${outfile}" "${outfolder_move}""${outfile}";
    echo "Output file stored in ${outfolder_move}${outfile}.";
    exit 0;
else
    echo "File ${outfolder}${outfile} was not produced by the configuration file.";
    exit 1;
fi
