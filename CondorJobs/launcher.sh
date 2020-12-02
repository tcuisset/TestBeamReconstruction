#!/usr/bin/env bash
declare -a ENERGIES=("20" "30" "50" "80" "100" "120" "150" "200" "250" "300")
declare -a DATATYPES=("data" "sim_proton" "sim_noproton" "sim_cmssw")
declare -a SHOWERTYPES=("em" "had")
declare -a STEPS=("selection" "analysis")
declare -a CELLTYPES=("LD" "HD")

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
ARGS=`getopt -o "" -l ",ntupleid:,step:,datatype:,showertype:,energy:,tag:,w0:,dpos:,celltype:" -n "getopts_${0}" -- "$@"`

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

	--celltype)
	    if [ -n "$2" ]; then
		if [[ " ${CELLTYPES[@]} " =~ " ${2} " ]]; then
		    CELLTYPE="${2}";
		    echo "Cell type: ${CELLTYPE}";
		else
		    echo "'--celltype' can be one of the following:"
		    printf "%s " "${CELLTYPES[@]}"
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
if [[ ( ( "${DATATYPE}" == "sim_proton" ) || ( "${DATATYPE}" == "sim_noproton" ) )
	&& ( "${NTUPLEID}" -gt 4 ) ]]; then
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
if [[ ( "${DATATYPE}" == "sim_cmssw" ) && ( "${STEP}" == "selection" ) ]]; then
    echo "The locally produced CMSSW simulation samples do not require any selection step."
    exit 1;
fi
if [[ ( "${DATATYPE}" == "sim_noproton" ) && ( "${SHOWERTYPE}" == "had" ) ]]; then
    echo "There is no proton-free sample for hadronic showers."
    exit 1;
fi
if [[ ( "${DATATYPE}" == "sim_cmssw" ) && ( "${SHOWERTYPE}" != "em" ) ]]; then
    echo "The locally produced CMSSW simulation samples only took into account the EE section."
    exit 1;
fi
if [[ ( ( "${DATATYPE}" == "sim_cmssw" ) && ( -z "${CELLTYPE}" ) )
       || ( ( "${DATATYPE}" != "sim_cmssw" ) && ( ! -z "${CELLTYPE}" ) ) ]]; then
    echo "The cell type should only and always be specified with '--datatype=sim_cmssw'."
    exit 1;
fi
if [[ -z "${ENERGY}" ]]; then
    echo "Please specify the beam energy."
    printf "Accepted values are: "
    printf "%s " "${ENERGIES[@]}"
    printf "[GeV].\n"
    exit 1;
fi
if [[ ( "${DATATYPE}" == "sim_cmssw" ) && ( "${ENERGY}" -ne 50 ) ]]; then
    echo "The locally produced CMSSW simulation samples were produced for 50GeV em showers."
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
export SCRAM_ARCH="slc7_amd64_gcc820"

if [ $(varExists "${INIT_FOLDER}") = true ] && [ $(varExists "${CMSSW_PATH}") = true ] &&
    [ $(varExists "${HOME_DIR}") = true ] && [ $(varExists "${ANALYSIS_PATH}") = true ]; then
    INIT_FOLDER=$(pwd);
    ANALYSIS_PATH="${CMSSW_BASE}/src/";
else
    echo "Use different variable names.";
    exit 0;
fi

cd "${ANALYSIS_PATH}";
source /afs/cern.ch/cms/cmsset_default.sh
eval `scramv1 runtime -sh` #cmsenv substitute

#back to the job folder
cd "${INIT_FOLDER}";

if [[ "${STEP}" == "selection" ]]; then

    if [[ "${DATATYPE}" == "data" ]]; then
	if [[ "${SHOWERTYPE}" == "em" ]]; then
	    INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_${NTUPLEID}.root"; #HGCAL only
	elif [[ "${SHOWERTYPE}" == "had" ]]; then
	    INFILE="/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ahcal-hgcal-merged-ntuples/ahcal_v8-hgcal_v16/merged_ntuple_${NTUPLEID}.root"; #HGCAL+AHCAL
	fi
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

    printf "Command:\n"
    COMMAND="process_data_exe ${INFILE} ${OUTFILE} ${DATATYPE} ${SHOWERTYPE} ${ENERGY}";
    echo ${COMMAND};
    ${COMMAND};

elif [[ "${STEP}" == "analysis" ]]; then

    OUTNAME="outEcut"
    if [[ "${DATATYPE}" == "data" ]]; then
	INFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_${NTUPLEID}.root";
	INTREE="relevant_branches"
	CLEAN=1
    elif [[ "${DATATYPE}" == "sim_noproton" ]]; then
	INFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root"
	INTREE="relevant_branches"
	CLEAN=1
    elif [[ "${DATATYPE}" == "sim_proton" ]]; then
	INFILE="/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGY}_${NTUPLEID}.root"
	INTREE="relevant_branches"
	CLEAN=1
    elif [[ "${DATATYPE}" == "sim_cmssw" ]]; then
	INFILE="/eos/user/b/bfontana/SinglePhoton/${CELLTYPE}/sim_cmssw_${NTUPLEID}.root"
	INTREE="ntuplizer/relevant_branches"
	CLEAN=0
    fi

    EOS_PATH="/eos/user/b/bfontana/TestBeamReconstruction/${TAG}/"
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
    echo "Input tree: ${INTREE}"
    echo -e "Output files:\n${OUTFILE1}\n${OUTFILE2}\n${OUTFILE3}"

    printf "Command:\n"
    COMMAND="analyze_data_exe ${INFILE} ${OUTFILE1} ${OUTFILE2} ${OUTFILE3} ${INTREE} ${SHOWERTYPE} ${W0} ${DPOS} ${CLEAN}";
    echo ${COMMAND};
    ${COMMAND};

fi
echo "##########################"
