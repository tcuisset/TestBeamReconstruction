#!/usr/bin/env bash
declare -a ENERGIES=("20" "30" "50" "80" "100" "120" "150" "200" "250" "300")
declare -a DATATYPES=("data" "sim_proton_v3" "sim_proton_v7" "sim_noproton_v5" "sim_noproton_v6")
declare -a SHOWERTYPES=("em" "had")
declare -a STEPS=("selection" "analysis")
#Note : putting /grid_mnt/.... is necessary, just using /home/llr leads to failures
export X509_USER_PROXY=/grid_mnt/vol_home/llr/cms/cuisset/.t3/proxy.cert
OUTPUT_FOLDER="/grid_mnt/data_cms_upgrade/cuisset/testbeam18/ntuple-selection/v2" #Default output folder

varExists() { 
    # Checks whether a certain environment variable already exists
    # Arguments:
    # 1. Variable being checked
	#Returns true if variable does *not* exist
    local flag=false;
    if [ -z "${1}" ]; then
	flag=true;
    fi
    echo $flag;
}

##########################
########PARSING###########
##########################
ARGS=`getopt -o "" -l ",ntupleid:,step:,datatype:,showertype:,energy:,tag:,w0:,dpos:,outputfolder::" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi
eval set -- "$ARGS"
echo "##### Input options: #####"
while true; do
    case "$1" in
	--step)
	    if [ -n "$2" ]; then
		if [[ " ${STEPS[@]} " =~ " ${2} " ]]; then
		    STEP="${2}";
		    echo "Step: ${STEP}";
		else
		    echo "'--step' can be one of the following:"
		    printf "%s " "${STEPS[@]}"
		    exit 1;
		fi
	    fi
	    shift 2;;

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
		    printf "%s " "${DATATYPES[@]}"
		    exit 1;
		fi
	    fi
	    shift 2;;

	--showertype)
	    if [ -n "$2" ]; then
		if [[ " ${SHOWERTYPES[@]} " =~ " ${2} " ]]; then
		    SHOWERTYPE="${2}";
		    echo "Data type: ${SHOWERTYPE}";
		else
		    echo "'--showertype' can be one of the following:"
		    printf "%s " "${SHOWERTYPES[@]}"
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

	--tag)
	    if [ -n "$2" ]; then
		TAG="${2}";
		echo "Tag: ${TAG}";
	    fi
	    shift 2;;

	--w0)
	    if [ -n "$2" ]; then
		W0="${2}";
		echo "w0 (cluster position measurement): ${W0}";
	    fi
	    shift 2;;

	--dpos)
	    if [ -n "$2" ]; then
		DPOS="${2}";
		echo "dpos (cluster position measurement): ${DPOS}";
	    fi
	    shift 2;;
	
	--outputfolder)
		if [ -n "$2" ]; then
		OUTPUT_FOLDER="$2"
		echo "Output folder: ${OUTPUT_FOLDER}";
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
if [[ -z "${STEP}" ]]; then
    echo "Please specify the step to run."
    printf "Accepted values are: "
    printf "%s " "${STEPS[@]}"
    printf "\n"
    exit 1;
fi
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
if [[ (( "${DATATYPE}" == "sim_noproton_v5" ) || ( "${DATATYPE}" == "sim_noproton_v6" ) ) && ( "${SHOWERTYPE}" == "had" ) ]]; then
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
if [[ ( -z "${TAG}" ) && ( "${STEP}" == "analysis" ) ]]; then
    echo "Please specify the tag when running the analysis step."
    exit 1;
fi
if [[ ( -z "${W0}" ) && ( "${STEP}" == "analysis" ) ]]; then
    echo "Please specify the 'w0' parameter when running the analysis step, to fully specify the cluster position measurement algorithm."
    exit 1;
fi
if [[ ( -z "${DPOS}" ) && ( "${STEP}" == "analysis" ) ]]; then
    echo "Please specify the 'dpos' parameter when running the analysis step, to fully specify the cluster position measurement algorithm."
    exit 1;
fi

##########################
##########################
##########################
export XRD_NETWORKSTACK=IPv4
export q="slc7_amd64_gcc820"

# What is this supposed to do ?
# if [ $(varExists "${INIT_FOLDER}") = true ] && [ $(varExists "${CMSSW_PATH}") = true ] &&
#     [ $(varExists "${HOME_DIR}") = true ] && [ $(varExists "${ANALYSIS_PATH}") = true ]; then
    
#     ANALYSIS_PATH="/afs/cern.ch/user/${USER:0:1}/${USER}/TestBeamAnalysis/src/";
# else
#     echo "Use different variable names.";
#     exit 0;
# fi

# cd "${ANALYSIS_PATH}";
# source /afs/cern.ch/cms/cmsset_default.sh
# eval `scramv1 runtime -sh` #cmsenv substitute

# #back to the job folder
# cd "${INIT_FOLDER}";

#Either use eos mounted on /eos (need kerberos ticket on LLR T3)
#INPUT_FILE_FOLDER="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/"
#Or use xrootd (need grid certificate on the node)
INPUT_FILE_FOLDER="root://eoscms.cern.ch///eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/"

if [[ "${STEP}" == "selection" ]]; then

    if [[ "${DATATYPE}" == "data" ]]; then
		if [[ "${SHOWERTYPE}" == "em" ]]; then
			INFILE="$INPUT_FILE_FOLDER/ntuples/v16/ntuple_${NTUPLEID}.root"; #HGCAL only
		elif [[ "${SHOWERTYPE}" == "had" ]]; then
			INFILE="$INPUT_FILE_FOLDER/ahcal-hgcal-merged-ntuples/ahcal_v8-hgcal_v16/merged_ntuple_${NTUPLEID}.root"; #HGCAL+AHCAL
		fi
		OUTFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_${NTUPLEID}.root";
	elif [[ "${DATATYPE}" == "sim_noproton_v5" ]]; then
		INFILE="$INPUT_FILE_FOLDER/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v5/electrons/ntuple_sim_config22_pdgID11_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root";
		OUTFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
	elif [[ "${DATATYPE}" == "sim_noproton_v6" ]]; then
		INFILE="$INPUT_FILE_FOLDER/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v6_33m/electrons/ntuple_sim_config22_pdgID11_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root";
		OUTFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
	elif [[ "${DATATYPE}" == "sim_proton_v3" ]]; then
		if [[ "${SHOWERTYPE}" == "em" ]]; then
			INFILE="$INPUT_FILE_FOLDER/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v3/electrons/ntuple_sim_config22_pdgID11_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root";
		elif [[ "${SHOWERTYPE}" == "had" ]]; then
			INFILE="$INPUT_FILE_FOLDER/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v44_VtxBeam_v3/CorrectFHLay10/pions/ntuple_sim_config22_pdgID211_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root"
		fi
		OUTFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
    elif [[ "${DATATYPE}" == "sim_proton_v7" ]]; then
		if [[ "${SHOWERTYPE}" == "em" ]]; then
			INFILE="$INPUT_FILE_FOLDER/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v7_33m/electrons/ntuple_sim_config22_pdgID11_beamMomentum${ENERGY}_listFTFP_BERT_EMN_0000_${NTUPLEID}.root";
		elif [[ "${SHOWERTYPE}" == "had" ]]; then
			echo "sim_proton_v7 and hadronic showers : I do not know which simulation to use"
			exit 1
		fi
		OUTFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
	fi
    echo "Input file: ${INFILE}"
    echo "Output file: ${OUTFILE}"
    process_data_exe "${INFILE}" "${OUTFILE}" "${DATATYPE}" "${SHOWERTYPE}" "${ENERGY}";

elif [[ "${STEP}" == "analysis" ]]; then

    OUTNAME="outEcut"
    if [[ "${DATATYPE}" == "data" ]]; then
		INFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_${NTUPLEID}.root";
    else
		INFILE="${OUTPUT_FOLDER}/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root"
    fi

    EOS_PATH="${OUTPUT_FOLDER}/${TAG}/"
    mkdir -p "${EOS_PATH}"

    HITFOLDER="hit_dependent/"
    LAYERFOLDER="layer_dependent/"
    CLUSTERFOLDER="cluster_dependent/"
    mkdir -p "${EOS_PATH}${HITFOLDER}"
    mkdir -p "${EOS_PATH}${LAYERFOLDER}"
    mkdir -p "${EOS_PATH}${CLUSTERFOLDER}"
    
    OUTFILE1="${EOS_PATH}${HITFOLDER}${OUTNAME}_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.csv"; 
    OUTFILE2="${EOS_PATH}${LAYERFOLDER}${OUTNAME}_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";
    OUTFILE3="${EOS_PATH}${CLUSTERFOLDER}${OUTNAME}_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root";

    echo "Input file: ${INFILE}"
    echo -e "Output files:\n${OUTFILE1}\n${OUTFILE2}\n${OUTFILE3}"
    analyze_data_exe "${INFILE}" "${OUTFILE1}" "${OUTFILE2}" "${OUTFILE3}" "${SHOWERTYPE}" "${W0}" "${DPOS}";

fi
