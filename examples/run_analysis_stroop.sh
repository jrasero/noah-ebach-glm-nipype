analysis_type=$1

python /home/javier/Documentos/software/noah-ebach-analysis/noah_ebach_analysis.py \
-b /media/javier/data/NOAH/BIDS \
-f /media/javier/data/NOAH/preprocessed/fmriprep_syn_sdc_w_reconall \
-o /media/javier/data/NOAH/analyses/motion_acompcor_cosines_fwd_syn_sdc \
--ncpus 10 --config_file /media/javier/data/NOAH/analyses/config_stroop.json stroop $analysis_type

