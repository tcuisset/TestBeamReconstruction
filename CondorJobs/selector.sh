#!/usr/bin/env bash
declare -a ENERGIES=("20" "30" "50" "80" "100" "120" "150" "200" "250" "300")
declare -a DATATYPES=("data" "sim_proton" "sim_noproton")
declare -a SHOWERTYPES=("em" "had")

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
ARGS=`getopt -o "" -l ",showertype:,datatype:,ntupleid:,energy:" -n "getopts_${0}" -- "$@"`

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
		
	--datatype)
	    if [ -n "$2" ]; then
		if [[ " ${DATATYPES[@]} " =~ " ${2} " ]]; then
		    DATATYPE="${2}";
		    echo "Data type: ${DATATYPE}";
		else
		    echo "'--datatype' can be one of the following:"
		    echo "sim_proton / sim_noproton / data"
		    exit 1;
		fi
	    fi
	    shift 2;;

	--showertype)
	    if [ -n "$2" ]; then
		if [[ " ${SHOWERTYPES[@]} " =~ " ${2} " ]]; then
		    DATATYPE="${2}";
		    echo "Data type: ${SHOWERTYPE}";
		else
		    echo "'--showertype' can be one of the following:"
		    echo "em / had"
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
if [[ -z "${SHOWERTYPE}" ]]; then
    echo "Please specify the shower type."
    printf "Accepted values are: "
    printf "%s " "${SHOWERTYPES[@]}"
    printf "\n"
    exit 1;
fi
if [[ ( "${DATATYPE}" == "sim_noproton" ) && ( "${SHOWERTYPE}" == "had" ) ]]; then
    echo "There is no proton-free sample for hadronic showers."
    exit 1;
fi
if [[ -z "${ENERGY}" ]]; then
    echo "Please specify the beam energy."
    printf "Accepted values are: "
    printf "%s " "${ENERGIES[@]}"
    printf "[GeV].\n"
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

if [[ "${DATATYPE}" == "data" ]]; then
    INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_${NTUPLEID}.root";
    OUTFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_${NTUPLEID}.root";
elif [[ "${DATATYPE}" == "sim_noproton" ]]; then
    INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v5/electrons/ntuple_sim_config22_pdgID11_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root";
    OUTFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
elif [[ "${DATATYPE}" == "sim_proton" ]]; then
    if [[ "${SHOWERTYPE}" == "em" ]]; then
	INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v3/electrons/ntuple_sim_config22_pdgID11_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root";
    elif [[ "${SHOWERTYPE}" == "had" ]]; then
       	INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v44_VtxBeam_v3/CorrectFHLay10/pions/ntuple_sim_config22_pdgID211_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root"
    fi
    OUTFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
fi
echo "Input file: ${INFILE}"
echo "Output file: ${OUTFILE}"
process_data_exe "${INFILE}" "${OUTFILE}" "${DATATYPE}" "${SHOWERTYPE}" "${ENERGY}";
