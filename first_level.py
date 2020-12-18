def create_first_level_wf(name,
                          output_dir,
                          preproc_bold,
                          brain_mask,
                          events_file,
                          confounds_file,                          
                          contrasts,
                          repetition_time,
                          start_ix,# = None,
                          fwhm,# = None,
                          thigh_pass,# = None,
                          confounds,
                          twenty_four=False):
    
    import nipype.pipeline.engine as pe 
    from nipype.interfaces import utility, io

    from nipype.algorithms import modelgen 
    from nipype.interfaces import fsl
    from niflow.nipype1.workflows.fmri.fsl.preprocess import create_susan_smooth
    
    from utils import create_subject_info
    
    first_level_wf = pe.Workflow(name = name)
    
    # Node to collect the inputs as explained above
    inputNode = pe.Node(interface=utility.IdentityInterface(fields=["bids_evs_file",
                                                                    "confounds_file",
                                                                    "in_func",
                                                                    "brain_mask"]),
                     name = "inputSource")
    
    # set inputs
    inputNode.inputs.in_func = preproc_bold
    inputNode.inputs.brain_mask = brain_mask
    inputNode.inputs.confounds_file = confounds_file
    inputNode.inputs.bids_evs_file = events_file
    
    
    # Node to drop first volumes (if none, set start_ix=0?)
    #if start_ix is None:
     #   start_ix = 0
        
    extractroi = pe.Node(fsl.ExtractROI(t_size=-1, t_min = start_ix),
                             name = "dropFirstVols")
    
    first_level_wf.connect(inputNode, "in_func", extractroi, "in_file")
    
    # Node to produce subject info needed to specify the model
    subject_info  = pe.Node(name="subject_info",
                            interface= utility.Function(input_names=["bids_evs_file", 
                                                                     "confounds_file",
                                                                     "confounds",
                                                                     "start_ix",
                                                                     "twenty_four"],
                                                        output_names = ["out_subj_info"],
                                                        function = create_subject_info)
                            )
    
    subject_info.inputs.start_ix = start_ix
    subject_info.inputs.confounds = confounds
    subject_info.inputs.twenty_four = twenty_four

   
    first_level_wf.connect(inputNode, "bids_evs_file", subject_info, "bids_evs_file")
    first_level_wf.connect(inputNode, "confounds_file", subject_info, "confounds_file")
    
    
    mask_bold = pe.Node(interface= fsl.maths.ApplyMask(),
                        name="mask_bold")

    first_level_wf.connect(extractroi, "roi_file", mask_bold, "in_file")
    first_level_wf.connect(inputNode, "brain_mask", mask_bold, "mask_file")
    
    # SUSAN smoothing using nipye predefined workflow
    susan = create_susan_smooth()
    susan.inputs.inputnode.fwhm = fwhm
    
    first_level_wf.connect(mask_bold, "out_file", susan, "inputnode.in_files")
    first_level_wf.connect(inputNode, "brain_mask", susan, "inputnode.mask_file")
     
    # Node to specify the FSL Model
    modelspec = pe.Node(modelgen.SpecifyModel(parameter_source='FSL',
					      input_units='secs',
                                              high_pass_filter_cutoff=thigh_pass,
                                              time_repetition = repetition_time),
                     name="modelspec")
    
    first_level_wf.connect(susan, "outputnode.smoothed_files", modelspec, "functional_runs")
    first_level_wf.connect(subject_info, "out_subj_info", modelspec, "subject_info")
        
    # Node to generate the FEAT specific files
    #TODO: derivatives (try on again 1 sep)?
    level1design = pe.Node(fsl.Level1Design(bases={'dgamma':{'derivs': True}},
                                            model_serial_correlations=True,
                                            contrasts = contrasts,
                                            interscan_interval = repetition_time),
                        name="level1design")
        
    first_level_wf.connect(modelspec, "session_info", level1design, "session_info")
    
    # Node to generate design.mat files
    feat = pe.Node(fsl.FEATModel(), 
                   name="FEAT")
    
    first_level_wf.connect([(level1design, feat, [("ev_files",  "ev_files"),
                                                  ("fsf_files",  "fsf_file")
                                                  ])
    ])
    
    # This node converts a list into a single file
    select_smooth_file= pe.Node(utility.Select(index=0), 
                                name='select_smooth_file')
    first_level_wf.connect(susan,"outputnode.smoothed_files", select_smooth_file, "inlist")

    # Node that uses FSL film_gls command to fit a design matrix to voxel timeseries and compute the contrast maps
    modelestimate = pe.Node(interface= fsl.FILMGLS(threshold=100.),
        name='modelestimate')
    
    first_level_wf.connect(select_smooth_file, "out", modelestimate, "in_file")
    first_level_wf.connect(feat, "design_file", modelestimate, "design_file")
    first_level_wf.connect(feat, "con_file", modelestimate, "tcon_file")
    first_level_wf.connect(feat, "fcon_file", modelestimate, "fcon_file")
    
    # Node to output the contrast maps
    outputNode = pe.Node(interface=utility.IdentityInterface(fields=["betasmat",
                                                                     "fstats",
                                                                     "zstats",
                                                                     "tstats",
                                                                     "zfstats",
                                                                     "copes",
                                                                     "varcopes"]),
                     name = "outputSource")
    
    first_level_wf.connect(modelestimate, "param_estimates", outputNode, "betasmat")
    first_level_wf.connect(modelestimate, "zstats", outputNode, "zstats")
    first_level_wf.connect(modelestimate, "fstats", outputNode, "fstats")
    first_level_wf.connect(modelestimate, "tstats", outputNode, "tstats")
    first_level_wf.connect(modelestimate, "zfstats", outputNode, "zfstats")
    first_level_wf.connect(modelestimate, "copes", outputNode, "copes")
    first_level_wf.connect(modelestimate, "varcopes", outputNode, "varcopes")
    
    # Output node for results
    datasink = pe.Node(io.DataSink(base_directory = output_dir),
                    name="datasink")
        
     # Dump outputs to output_dir
    first_level_wf.connect([(outputNode, datasink, [("betasmat", "betas"),
                                                    ("tstats", "stats"),
                                                    ("zstats", "stats.@zstats"),
                                                    ("fstats", "stats.@fstats"),
                                                    ("zfstats", "stats.@zfstats"),
                                                    ("copes", "copes"),
                                                    ("varcopes", "varcopes")
                                              ])
        ])

    # This is just to have the design matrix used
    first_level_wf.connect(feat, 'design_image', datasink, 'design_image')

    return first_level_wf
