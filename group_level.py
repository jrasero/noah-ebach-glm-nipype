def create_group_level_wf(name,
                          output_dir,
                          copes, 
                          varcopes, 
                          mask, 
                          flame_mode = "flame1",
                          randomise=True,
                          n_perms=1000,
                          seed = None):
    
    """ 
    
    This function creates the workflow for doing the group level
    analysis from previous computations. In particular, we have to pass it
    the copes, varcopes, a mask at least
    
    """
    
    import nipype.pipeline.engine as pe 
    from nipype.interfaces import utility, io
    from nipype.interfaces import fsl

    group_level_wf = pe.Workflow(name = name)
    
    # Node to collect the inputs as explained above
    inputNode = pe.Node(utility.IdentityInterface(fields=["copes",
                                                          "varcopes",
                                                          "mask"]),
                     name = "inputSource")
    
    inputNode.inputs.copes = copes
    inputNode.inputs.varcopes = varcopes
    inputNode.inputs.mask = mask

    merge_copes = pe.Node(fsl.Merge(dimension="t"), 
                          name="merge_copes")
    merge_varcopes = pe.Node(fsl.Merge(dimension="t"),
                             name="merge_varcopes")
    
    design_matrix = pe.Node(fsl.L2Model(num_copes = len(copes)), 
                            name = "design_matrix")

    flame = pe.Node(fsl.FLAMEO(run_mode=flame_mode),
                    name= "flame")
    
    datasink = pe.Node(io.DataSink(base_directory=output_dir), 
                       name="datasink")
    
    group_level_wf.connect(inputNode, 'copes', merge_copes, 'in_files')
    group_level_wf.connect(inputNode, 'varcopes', merge_varcopes, 'in_files')
    group_level_wf.connect(merge_copes, 'merged_file', flame, 'cope_file')
    group_level_wf.connect(merge_varcopes, 'merged_file', flame, 'var_cope_file')
    group_level_wf.connect(design_matrix, 'design_mat', flame, 'design_file')
    group_level_wf.connect(design_matrix, 'design_con', flame, 't_con_file')
    group_level_wf.connect(design_matrix, 'design_grp', flame, 'cov_split_file')
    group_level_wf.connect(inputNode, 'mask', flame, 'mask_file')
    
    group_level_wf.connect([(flame, datasink, [('pes', 'flame'),
                                               ('tstats', 'flame.@tstats'),
                                               ('zstats', 'flame.@zstats'),
                                               ('zfstats', 'flame.@zfstats')])
                            ])
    
    # Correct these images?
    if randomise:
        
        #if seed is None:
        # Randomise using a threshold-free procedure    
        randomise = pe.Node(fsl.Randomise(num_perm = n_perms, 
                                          tfce=True,
                                          one_sample_group_mean=True), 
                            name="tfce")
        
        # Since this is a pure one-sided test, we to use the abs, 
        # and then retrieve the sign 
        abs_img = pe.Node(fsl.ImageMaths(op_string="-abs"), name="abs")
        
        group_level_wf.connect(merge_copes, 'merged_file', abs_img, 'in_file')
        group_level_wf.connect(abs_img, 'out_file', randomise, 'in_file')
        group_level_wf.connect(inputNode, 'mask', randomise, 'mask')
            
        group_level_wf.connect([(randomise, datasink, 
                                 [("t_corrected_p_files",'randomise'),
                                  ("tstat_files",'randomise.@tstat_files')])
                                ])
    
        
    return group_level_wf
