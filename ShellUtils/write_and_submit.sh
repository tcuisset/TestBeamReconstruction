declare -a TAGS=( "W5p0_dpos999" )
#                  "W2p9_dpos1p3"
#		  "W4p0_dpos1p3"
#		  "W5p0_dpos1p3"
#                 "W2p9_dpos3p4"
#                  "W2p9_dpos999" )

for i in `seq 0 $(expr "${#TAGS[@]}" - 1)`; do
    W0=${TAGS[i]%_dpos*} #remove everything starting from '_dpos'
    W0=${W0:1:${#W0}} #remove the initial 'W'
    W0=`echo ${W0} | sed 's/p/./g'` #replace 'p' by '.'v
    DPOS=${TAGS[i]##*_dpos} #remove everything before, and including, '_dpos'
    DPOS=`echo ${DPOS} | sed 's/p/./g'` #replace 'p' by '.'

    write_dag --datatype data --showertype em --w0 ${W0} --dpos ${DPOS} --tag ${TAGS[i]} --last_step_only;
    condor_submit_dag CondorJobs/clue_data_em_"${TAGS[i]}"_analysis_only.dag;

done
