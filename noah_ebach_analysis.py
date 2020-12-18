#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from pathlib import Path
from glob import glob

def get_parser():
    """Define the command line interface"""
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

    parser = ArgumentParser(description='NOAH/eBACH Analysis',
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-b', '--bids_dir', action='store', type=Path,
                        required=True, dest = "bids_dir",
                        help='BIDS folder')
    parser.add_argument('-f','--fmriprep_dir', action='store', type=Path,
                        required=True, dest = "fmriprep_dir",
        help= 'fMRIPrep folder')
    parser.add_argument('-o', '--output_dir', action='store', type=Path,
                        required=True, dest = "output_dir",
                        help='the output directory')
    parser.add_argument('task_id',  type=str,
                        choices=['stroop', 'msit', 'emoreap'],
                        help='select a specific task to be processed')
    parser.add_argument('analysis_level', 
                        choices=['participant', 'group'], 
                        help='type of analysis')    
    parser.add_argument('--participant_label', action='store', type=str,
                        nargs='*', help='process only particular subjects')
    parser.add_argument('--ncpus', action='store', type=int,
                         help='number of cpus')    
    parser.add_argument('-w', '--work_dir', action='store', type=Path,
                        dest = "work_dir",
                         help='path where intermediate results should be stored')
    
    parser.add_argument('--config_file', action='store', type=Path,
                        dest = "config_file",
                         help='path to config file')


    return parser    
    

def main():
    
    from bids.layout import BIDSLayout
    from templateflow.api import get as tpl_get
    
    from nipype.pipeline.engine import Workflow
    import json
    from os.path import join as opj
    
    #from nilearn import image
    from first_level import create_first_level_wf 
    from group_level import create_group_level_wf 
    from utils import (create_workflow_name, create_output_dir, 
                        get_contrasts, get_data_info, default_task_config)
    opts = get_parser().parse_args()
    
    
    bids_layout = BIDSLayout(opts.bids_dir.absolute().as_posix(), 
                             derivatives = opts.fmriprep_dir.absolute().as_posix()
                             )
    
    data_info = get_data_info(bids_layout)
    subject_list = data_info.subject_list
    session_list = data_info.session_list
    repetition_time = data_info.TR
    
    if opts.ncpus:
        run_config = dict(plugin = 'MultiProc',
                          plugin_args =  {'n_procs': opts.ncpus,
                                          'raise_insufficient': False,
                                          'maxtasksperchild': 1}
                          )
    else:
        run_config = dict(plugin = 'Linear')
            
    task_id = opts.task_id
    print("RUNNIN TASK = %s" % task_id)
    
    # If subject ids are supplied, do only for those
    if opts.participant_label:
        subject_list = list(set(subject_list) & set(opts.participant_label))
    
    output_dir = opts.output_dir
    
    #if output_dir.exists() is False:
     #   raise #
        
    if opts.work_dir:
        work_dir = opts.work_dir
    else:
        work_dir = Path(output_dir).joinpath("work")
        work_dir.mkdir(parents=True, exist_ok=True)
        work_dir = work_dir.absolute().as_posix()
        
    log_dir = Path(output_dir).joinpath("log/task-%s" % task_id)
    log_dir.mkdir(parents=True, exist_ok=True)

    contrasts = get_contrasts(task_id)
    
    with open(log_dir.joinpath("contrast.log"), "w") as f:
        for ii, contrast in enumerate(contrasts):
            f.write(("condition %d = ") % (ii + 1))
            f.write(str(contrasts[ii]) + "\n")
        
    if opts.config_file:
        #TODO: write function that checks the extension and field of this file
        with open(opts.config_file, "r") as f:
            config_task = json.load(f)
    else:
        config_task = default_task_config(task_id)
        
    with open(log_dir.joinpath("config.log"), "w") as f:
        for key, value in config_task.items():
            f.write(key + " = " )
            f.write(str(value) + "\n")
        
    ################### FIRST LEVEL PART##################
        
    query_task = {'preproc_bold': dict(suffix='bold', 
                             extension= ['nii', 'nii.gz'],
                             desc='preproc'), 
                  'confounds_file':dict(suffix='regressors', 
                                     extension= ['tsv', 'tsv.gz'],
                                     desc='confounds'),
                  'brain_mask': dict(suffix='mask', 
                                  datatype='func',
                                  extension= ['nii', 'nii.gz'],
                                  desc='brain'),
                  'events_file':dict(suffix='events', extension=['tsv'])
                      }
    
    config_first = config_task["config_first"]
    
    first_level_dir = output_dir.joinpath("first_level")
    first_level_dir.mkdir(parents=True, exist_ok=True)
    
    first_level_wf = Workflow(name="First-level")
    
    # this loops add 
    for subject_id in subject_list:
        for session_id in session_list:
        
            inputs_files = {}
            
            try:
                for key, query in query_task.items():
                    args = query.copy()
                    filelist = bids_layout.get(task=task_id, 
                                               subject=subject_id, 
                                               session=session_id, 
                                               **args)
                    inputs_files[key] = filelist[0].path
            except:
                print("some file no present for subject %s, task %s" % (subject_id, task_id))
                continue
    
            preproc_bold = inputs_files['preproc_bold']
            brain_mask = inputs_files['brain_mask']
            confounds_file = inputs_files['confounds_file']
            events_file = inputs_files['events_file']
        
            run_name = create_workflow_name(task_id, 
                                            subject_id, 
                                            session_id, 
                                            None)
            
            output_first_dir = create_output_dir(first_level_dir, 
                                                 task_id, 
                                                 subject_id, 
                                                 session_id, 
                                                 None)
            
            output_first_dir = output_first_dir.absolute().as_posix()
            
            individual_wf = create_first_level_wf(name=run_name,
                                                  output_dir=output_first_dir,
                                                  preproc_bold=preproc_bold,
                                                  brain_mask=brain_mask,
                                                  events_file=events_file,
                                                  confounds_file=confounds_file,
                                                  contrasts=contrasts,
                                                  repetition_time=repetition_time,
                                                  **config_first)
        
            if os.path.exists(output_first_dir) is False:
                print("adding first-level %s " % run_name)
                first_level_wf.add_nodes([individual_wf])
    
    # Run this
    first_level_wf.base_dir = work_dir
    first_level_wf.run(**run_config)
    
    print("first level analysis done!")
    
    ################### GROUP LEVEL PART##################
    if opts.analysis_level == "group":
        config_group = config_task["config_group"]

        group_level_wf = Workflow(name="Group-level")
    	
	#TODO:See how to pass this as input
        generate_group_mask = False
        if generate_group_mask is True:
	  #config_group["generate_group_mask"]:
            print("TODO")
           
            #TODO: Use what, the segmentation probs, the masks directly?
            
            # list_gm_masks = [glob(opj(fmriprep_dir, 
            #                           subj, 
            #                           "anat", 
            #                           "*_space-MNI152NLin2009cAsym_dseg.nii.gz"))[0] \
            #       for subj in subject_list]
    
            # gm_masks = [image.math_img("(img==2).astype(int)", img=gm_mask) for \
            #                 gm_mask in list_gm_masks]
            # mean_gm_mask = image.mean_img(gm_masks)
            
        else:
            # Use mask in standard space from template
            group_mask_file = tpl_get("MNI152NLin2009cAsym", 
                                      resolution=2, 
                                      desc='brain',suffix='mask').as_posix()
        
        
        #select copes and varcopes
        for ii, contrast in enumerate(contrasts):
            copes  = glob(opj(output_dir, 
                         "first_level", "task-" + task_id, 
                         "*", "ses-01", "copes",
                         "cope" + str(ii+1) + ".nii.gz"))
            
            varcopes = glob(opj(output_dir, 
                            "first_level", "task-" + task_id,
                            "*", "ses-01", "varcopes",
                            "varcope" + str(ii+1) + ".nii.gz"))
            
            group_level_dir = output_dir.joinpath("group_level/task-%s/cond_%d" % (task_id, 
                                                                                   (ii+1)))
            group_level_dir.mkdir(parents=True, exist_ok=True)
            group_level_dir = group_level_dir.absolute().as_posix()
        
            name = "group_level_task_%s_cond_%d_wf" % (task_id, (ii+1))
            
            cond_group_wf = create_group_level_wf(name,
                                                  group_level_dir,
                                                  copes,
                                                  varcopes, 
                                                  group_mask_file, 
                                                  **config_group)
            
            group_level_wf.add_nodes([cond_group_wf])
            
        
    
        
    # Run group level workflow
    group_level_wf.base_dir = work_dir
    group_level_wf.run(**run_config)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
    
