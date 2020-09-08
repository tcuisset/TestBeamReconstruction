#!/usr/bin/env bash
declare -a ENERGIES=("20" "30" "50" "80" "100" "120" "150" "200" "250" "300")
declare -a DATATYPES=("data" "sim_proton" "sim_noproton")
declare -a ANALYSISTYPES=("layerdep" "clusterdep")

##########################
########PARSING###########
##########################
ARGS=`getopt -o "" -l ",datatype:,analysistype:" -n "getopts_${0}" -- "$@"`

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
		    echo "sim_proton / sim_noproton / data"
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
		    echo "layerdep / clusterdep"
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
    hadd -f /eos/user/b/bfontana/TestBeamReconstruction/job_output/"${JOBSFOLDER}"/hadd_"${ANALYSISTYPE}"_"${DATATYPE}"_beamen${ENERGIES[j]}.root /eos/user/b/bfontana/TestBeamReconstruction/job_output/"${JOBSFOLDER}"/outEcut_"${DATATYPE}"*_beamen${ENERGIES[j]}_*.root;
done	
