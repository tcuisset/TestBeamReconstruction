#!/usr/bin/env bash
declare -a TAGS=( "W2p3_dpos1p3"
		  "W2p9_dpos1p3"
		  "W4p0_dpos1p3"
		  "W5p0_dpos1p3"
		  "W2p9_dpos3p4"
		  "W2p9_dpos999" )
declare -a VARS=( "dx"
		  "dy"
		  "dxdy"
		  "posx_posy" )
declare -a BEAM_ENERGY=( 20 30 50 80 100 120 150 200 250 300 )
USE_SAVED_DATA=0
BEAM_ENERGY=50

######################################
##Agument parsing#####################
######################################
ARGS=`getopt -o "" -l ",chosen_energy:,use_saved_data:,var:" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi
eval set -- "$ARGS"
echo "##### Input options: #####"
while true; do
    case "$1" in
	--var)
	    if [ -n "$2" ]; then
		if [[ ! " ${VARS[@]} " =~ " ${2} " ]]; then
		    echo "Var ${2} is not accepted."
		    exit 1;
		else
		    VAR="${2}";
		    echo "variable to consider: ${VAR}";
		fi
	    fi
	    shift 2;;

	--use_saved_data)
	    USE_SAVED_DATA="${2}"
	    if [ ${USE_SAVED_DATA} -eq 1 ]; then
		echo "Use saved data.";
	    fi
	    shift 2;;

	--chosen_energy)
	    if [ -n "$2" ]; then
		if [[ ! " ${BEAM_ENERGY[@]} " =~ " ${2} " ]]; then
		    echo "Beam energy ${2} is not part of the analysis."
		    exit 1;
		else
		    BEAM_ENERGY="${2}";
		    echo "variable to consider: ${BEAM_ENERGY}";
		fi
	    fi
	    shift 2;;
	
	--)
	    shift
	    break;;
    esac
done

######################################
##Run the code########################
######################################
for i in $(seq 0 $(expr "${#TAGS[@]}" - 1)); do
    COMMAND="python DataProcessing/python/cluster_dep.py --datatype data --showertype em --${VAR} --chosen_energy ${BEAM_ENERGY} --tag ${TAGS[${i}]}"
    if [ ${USE_SAVED_DATA} -eq 1 ]; then
	${COMMAND} --use_saved_data &
    else
	${COMMAND} &
    fi
    pids[${i}]=$!
done

######################################
##Wait for jobs to finish#############
######################################
for pid in ${pids[*]}; do
    wait $pid
done
echo "All jobs finished."
