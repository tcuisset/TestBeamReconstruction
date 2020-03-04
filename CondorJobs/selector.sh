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

INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_${1}.root";
echo "Input file: ${INFILE}"
OUTFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${1}.txt";
echo "Output file: ${OUTFILE}"
process_data_exe "${INFILE}" "${OUTFILE}";
