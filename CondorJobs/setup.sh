#!/usr/bin/env bash
declare -a CELLTYPES=("LD" "HD")

ARGS=`getopt -o "" -l ",celltype:" -n "getopts_${0}" -- "$@"`

#Bad arguments
if [ $? -ne 0 ];
then
  exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
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

	--)
	    shift
	    break;;
    esac
done

MAINPATH="${CMSSW_BASE}/src/UserCode/CondorJobs/"
OUTFILE="${MAINPATH}ntuple_sim_cmssw_${CELLTYPE}_ids.txt"
echo "Generating '${OUTFILE}' file..."
ls -l /eos/user/b/bfontana/SinglePhoton/${CELLTYPE}/ | awk 'NR>1 {sub("\\.root", "", $9); sub("sim_cmssw_", "", $9); print $9}' > ${OUTFILE}
