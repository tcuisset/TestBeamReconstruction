#!/usr/bin/env bash
declare -a ENERGIES=("20" "30" "50" "80" "100" "120" "150" "200" "250" "300")
declare -a DATATYPES=("data" "sim_proton" "sim_noproton")
declare -a SHOWERTYPES=("em" "had")
declare -a ANALYSISTYPES=("layerdep" "clusterdep")

##########################
########PARSING###########
##########################
ARGS=`getopt -o "" -l ",datatype:,showertype:,analysistype:" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi
eval set -- "$ARGS"
echo "##### Input options: #####"
while true; do
    case "$1" in
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

	--analysistype)
	    if [ -n "$2" ]; then
		if [[ " ${ANALYSISTYPES[@]} " =~ " ${2} " ]]; then
		    ANALYSISTYPE="${2}";
		    echo "Analysis type: ${ANALYSISTYPE}";
		else
		    echo "'--analysistype' can be one of the following:"
		    printf "%s " "${ANALYSISTYPES[@]}"
		    exit 1;
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
if [[ -z "${ANALYSISTYPE}" ]]; then
    echo "Please specify the analysis type."
    printf "Accepted values are: "
    printf "%s " "${ANALYSISTYPES[@]}"
    printf "\n"
    exit 1;
fi  
##########################
##########################
##########################

if [[ "${ANALYSISTYPE}" == "layerdep" ]]; then
    JOBSFOLDER="layer_dependent"
elif [[ "${ANALYSISTYPE}" == "clusterdep" ]]; then
    JOBSFOLDER="cluster_dependent"
fi
len="${#ENERGIES[@]}"
for(( j=0; j<${len}; j++ )); do
    IN="/eos/user/${USER:0:1}/${USER}/TestBeamReconstruction/job_output/${JOBSFOLDER}/outEcut_${DATATYPE}_${SHOWERTYPE}_beamen${ENERGIES[j]}_";
    if [[ $(ls "${IN}"*root -A) ]]; then #in case the input files do exist
	hadd -f /eos/user/${USER:0:1}/${USER}/TestBeamReconstruction/job_output/"${JOBSFOLDER}"/hadd_"${ANALYSISTYPE}"_"${DATATYPE}"_"${SHOWERTYPE}"_beamen${ENERGIES[j]}.root "${IN}"*root;
    fi
done
