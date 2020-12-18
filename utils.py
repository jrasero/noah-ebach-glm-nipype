def read_contrast(filepath):
    
    import json, os
    filepath = os.path.abspath(filepath)
    #print(filepath)
    try:    
        with open(filepath, 'r') as f:
            datafile = json.load(f)
            
            contrasts = []
            for key in datafile.keys():
                contrasts.append(datafile[key])
            
            return contrasts
    except NameError:
        print("File does not exist")    
        
        

def create_subject_info(bids_evs_file, 
                        confounds_file, 
                        confounds, 
                        start_ix,
                        twenty_four=False):
    
    import pandas as pd
    import numpy as np
    from nipype.algorithms.modelgen import bids_gen_info
    
    
    info = bids_gen_info([bids_evs_file], condition_column = "trial_type")
    
    confounds_df = pd.read_csv(confounds_file, sep="\t", na_values="n/a")
    
    try:
        confounds_df = confounds_df.loc[:, confounds]
        if confounds_df.isna().any(None):
            print("There are {:d} rows with NaNs that we fill with a 0 value".format(sum(confounds_df.isna().any(axis=1))))
            confounds_df = confounds_df.fillna(0)
    except:
        RuntimeError("Some confounding variables do not exist")
        
    if twenty_four:
        
        motion_cols = ['trans_x','trans_y','trans_z', 
                     'rot_x', 'rot_y', 'rot_z']
        # Take only motion parameters' data
        data = confounds_df.loc[:, motion_cols].values
        # Compute square of this
        motion_cols_sq = [col + "_sq" for col in motion_cols]
        data_sq_df = pd.DataFrame(data**2, columns = motion_cols_sq)
        # Roll one position
        data_roll = np.roll(data, 1, axis=0)
        data_roll[0] = 0
        motion_roll_cols = [col + "_dt" for col in motion_cols]
        data_roll_df = pd.DataFrame(data_roll, columns = motion_roll_cols)
        # Compute square of rolled data
        motion_roll_sq_cols = [col + "_sq_dt" for col in motion_cols]
        data_roll_sq_df = pd.DataFrame(data_roll**2, columns = motion_roll_sq_cols)
        # Concatenate new data to original data frame 
        confounds_df = pd.concat([confounds_df,data_sq_df,data_roll_df, data_roll_sq_df],
                                 axis=1)
        # Redefine confound names
        confounds = confounds + motion_cols_sq + motion_roll_cols + motion_roll_sq_cols
        
    info[0].update(regressors = [confounds_df.loc[start_ix: , name].to_list() \
                                 for name in confounds],
                   regressor_names = confounds)
    
    return info


def create_workflow_name(task_id, subject_id, session_id, run_id):
    
    name = 'task_' + task_id + '_sub_'  + subject_id 
    
    if session_id is not None:
        name += '_ses_' + session_id
        
        if run_id is not None:
             name += '_run_' + run_id
    else:
        if run_id is not None:
             name += '_run_' + run_id
    return name
   
def create_output_dir(base_dir, task_id, subject_id, session_id, run_id):
    
    output_dir = base_dir.joinpath('task-' + task_id, 'sub-' + subject_id)
    
    if session_id is not None:
        
        output_dir = output_dir.joinpath('ses-' + session_id)
        
        if run_id is not None:
            output_dir = output_dir.joinpath('run-' + run_id)
    else:
        if run_id is not None:
            output_dir = output_dir.joinpath('run-' + run_id)

    return output_dir
   
def get_data_info(bids_layout):
    
    from collections import namedtuple

    raw_subject_list = bids_layout.get_subjects(scope='raw')
    fmriprep_subject_list = bids_layout.get_subjects(scope='derivatives')
    subject_list = list(set(raw_subject_list) & set(fmriprep_subject_list)) 

    # take all sessions both raw and preprocessed
    raw_session_list = bids_layout.get_sessions(scope='raw')
    fmriprep_session_list = bids_layout.get_sessions(scope='derivatives')
    session_list = list(set(raw_session_list) & set(fmriprep_session_list)) 
            
    #raw_run_list = bids_layout.get_runs(scope='raw')
    #fmriprep_run_list = bids_layout.get_runs(scope='derivatives')
    #run_list = list(set(raw_run_list) & set(fmriprep_run_list)) 

    #if not session_list:
     #   session_list = [None]
    #if not run_list:
      #  run_list = [None]

    repetition_time = bids_layout.get_tr()
    
    data_info = namedtuple("DataInfo", ["subject_list", 
                                        "session_list", 
                                        "TR"])
    
    return data_info(subject_list, session_list, repetition_time)

def get_contrasts(task_id):
    if task_id == "msit":
        
        # Task conditions
        cong_cond = ['Congruent', 'T', ['Congruent'], [1]]
        incong_cond = ['Incongruent', 'T', ['Incongruent'], [1]]
        
        cong_vs_inc = ['Congruent > Incongruent', 'T',  
                       ['Congruent', 'Incongruent'], [1, -1]]
        
        contrasts= [cong_cond, incong_cond, cong_vs_inc]
    
        
    elif task_id == "stroop":
        
        cong_cond = ['Congruent', 'T', ['Congruent'], [1]]
        incong_cond = ['Incongruent', 'T', ['Incongruent'], [1]]
        
        cong_vs_inc = ['Congruent > Incongruent', 'T',  
                       ['Congruent', 'Incongruent'], [1, -1]]
        
        contrasts= [cong_cond, incong_cond, cong_vs_inc]
    
    elif task_id == "emoreap":
        
        lookneg_cond = ['LookNeg', 'T', ['LookNeg'],[1]]
        lookneut_cond = ['LookNeut', 'T', ['LookNeut'],[1]]
        regneg_cond = ['RegNeg', 'T', ['RegNeg'],[1]]
        
        lookneg_vs_lookneut = ['LookNeg > LookNeut', 'T', ['LookNeg', 'LookNeut'], [1, -1]]
        regneg_vs_lookneg = ['RegNeg > LookNeg', 'T', ['RegNeg', 'LookNeg'], [1, -1]]
        regneg_vs_lookneut = ['RegNeg > LookNeut', 'T', ['RegNeg', 'LookNeut'], [1, -1]]
        
        contrasts= [lookneg_cond, lookneut_cond, regneg_cond,
                    lookneg_vs_lookneut, regneg_vs_lookneg, regneg_vs_lookneut
                    ]
    
    else:
        ValueError("")
        
    return contrasts


def default_task_config(task_id):
    """
    
    Function that will give a default configuration pipeline in case 
    no json file is supplied. Based on previously run analyses
    
    """
    # Default confounders = motion parameters + acompcor components (+cosines?) + fwd
    confounds =  ["trans_x","trans_y","trans_z", 
                  "rot_x", "rot_y", "rot_z",
                  "a_comp_cor_00", "a_comp_cor_01", "a_comp_cor_02",
                  "a_comp_cor_03", "a_comp_cor_04", "a_comp_cor_05",
		  "cosine00", "cosine01", "cosine02", "cosine03",
		  "cosine04", "cosine05","cosine06","cosine07","cosine08",
                  "framewise_displacement"]
    
    RANDOM_STATE = 69
    config = dict()
    
    if task_id=="emoreap":
        config["config_first"] = {"thigh_pass" : 128.0,
                                  "fwhm" : 6.0, 
                                  "start_ix" : 0,
                                  "twenty_four" : False,
                                  "confounds" : confounds}
    else:
        config["config_first"] = {"thigh_pass" : 187.0,
                                  "fwhm" : 6.0, 
                                  "start_ix" : 0,
                                  "twenty_four" : False,
                                  "confounds" : confounds}
        
        
    config["config_group"] = {"flame_mode" : "flame1",
                               "randomise" : True,
                                "n_perms" : 10000, 
                                "seed" : RANDOM_STATE}
        
    return config



