#!/usr/bin/env bash
declare -a ENERGIES=("20" "30" "50" "80" "100" "120" "150" "200" "250" "300")
declare -a DATATYPES=("data" "sim_proton" "sim_noproton")

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

##########################
########PARSING###########
##########################
ARGS=`getopt -o "" -l ",data_type:,ntupleid:,energy:" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi
eval set -- "$ARGS"
echo "##### Input options: #####"
while true; do
    case "$1" in
	--ntupleid)
	    if [ -n "$2" ]; then
		NTUPLEID="${2}";
		echo "Ntuple id: ${NTUPLEID}";
	    fi
	    shift 2;;
	
	--data_type)
	    if [ -n "$2" ]; then
		if [[ " ${DATATYPES[@]} " =~ " ${2} " ]]; then
		    DATATYPE="${2}";
		    echo "Data type: ${DATATYPE}";
		else
		    echo "'--data_type' can be one of the following:"
		    echo "sim_proton / sim_noproton / data"
		    exit 1;
		fi
	    fi
	    shift 2;;

	--energy)
	    if [ -n "$2" ]; then
		if [[ ! " ${ENERGIES[@]} " =~ " ${2} " ]]; then
		    echo "Energy with value ${2} is not part of the analysis."
		    exit 1;
		else
		    ENERGY="${2}";
		    echo "Beam energy: ${ENERGY} GeV";
		fi
	    fi
	    shift 2;;
	
	--)
	    shift
	    break;;
    esac
done
echo "##########################"
##########################
##########################
##########################

##########################
########ARG CHECKS########
##########################
if [[ -z "${DATATYPE}" ]]; then
    echo "Please specify the data type."
    printf "Accepted values are: "
    printf "%s " "${DATATYPES[@]}"
    printf "\n"
    exit 1;
fi  
if [[ ( "${DATATYPE}" == *"sim"* ) && ( "${NTUPLEID}" -gt 4 ) ]]; then
    echo "Simulation data has Ntuples numbered from 0 to 4."
    exit 1;
fi
if [[ ( "${DATATYPE}" == *"sim"* ) && ( -z "${ENERGY}" ) ]]; then
    echo "Simulation data requires specifying the beam energy."
    printf "Accepted values are: "
    printf "%s " "${ENERGIES[@]}"
    printf "[GeV].\n"
    exit 1;
fi
if [[ ( "${DATATYPE}" == "data" ) && ( ! -z "${ENERGY}" ) ]]; then
    echo "Real data does not require specifying the beam energy. Please remove it."
    exit 1;
fi
##########################
##########################
##########################

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

EOS_PATH="/eos/user/b/bfontana/TestBeamReconstruction/job_output/"
OUTNAME="outEcut" #extract ntuple number
if [[ "${DATATYPE}" == "data" ]]; then
    INFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${NTUPLEID}.root";
    OUTFILE1="${EOS_PATH}hit_dependent/${OUTNAME}.csv"; 
    OUTFILE2="${EOS_PATH}layer_dependent/${OUTNAME}.root";
    OUTFILE3="${EOS_PATH}cluster_dependent/${OUTNAME}.root";
elif [[ "${DATATYPE}" == "sim_noproton" ]]; then
    INFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.root"
    OUTFILE1="${EOS_PATH}hit_dependent/${OUTNAME}_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.csv"; 
    OUTFILE2="${EOS_PATH}layer_dependent/${OUTNAME}_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.root";
    OUTFILE3="${EOS_PATH}cluster_dependent/${OUTNAME}_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.root";
elif [[ "${DATATYPE}" == "sim_proton" ]]; then
    INFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.root"
    OUTFILE1="${EOS_PATH}hit_dependent/${OUTNAME}_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.csv"; 
    OUTFILE2="${EOS_PATH}layer_dependent/${OUTNAME}_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.root";
    OUTFILE3="${EOS_PATH}cluster_dependent/${OUTNAME}_${DATATYPE}_beamen${ENERGY}_${NTUPLEID}.root";
fi

echo "Input file: ${INFILE}"
echo -e "Output files:\n${OUTFILE1}\n${OUTFILE2}\n${OUTFILE3}"
analyze_data_exe "${INFILE}" "${OUTFILE1}" "${OUTFILE2}" "${OUTFILE3}";
