write_dag --datatype data --showertype em --w0 2.3 --dpos 1.3 --tag W2p3_dpos1p3 --last_step_only;
condor_submit_dag CondorJobs/clue_data_em_W2p3_dpos1p3_analysis_only.dag;

write_dag --datatype data --showertype em --w0 2.9 --dpos 1.3 --tag W2p9_dpos1p3 --last_step_only;
condor_submit_dag CondorJobs/clue_data_em_W2p9_dpos1p3_analysis_only.dag;

write_dag --datatype data --showertype em --w0 4.0 --dpos 1.3 --tag W4p0_dpos1p3 --last_step_only;
condor_submit_dag CondorJobs/clue_data_em_W4p0_dpos1p3_analysis_only.dag;

write_dag --datatype data --showertype em --w0 5.0 --dpos 1.3 --tag W5p0_dpos1p3 --last_step_only;
condor_submit_dag CondorJobs/clue_data_em_W5p0_dpos1p3_analysis_only.dag;

write_dag --datatype data --showertype em --w0 2.9 --dpos 3.4 --tag W2p9_dpos3p4 --last_step_only;
condor_submit_dag CondorJobs/clue_data_em_W2p9_dpos3p4_analysis_only.dag;

write_dag --datatype data --showertype em --w0 2.9 --dpos 999.0 --tag W2p9_dpos999 --last_step_only;
condor_submit_dag CondorJobs/clue_data_em_W2p9_dpos999_analysis_only.dag;
