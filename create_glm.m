function glm = create_glm(glm_var)
glm.subrun_length = 2*glm_var.restdur+2*glm_var.blocks*(glm_var.blockdur+glm_var.restdur); 
glm.hrf_distri = spm_hrf(glm_var.TR)';
contexts = {'taskdplus', 'taskdminus'};
for j=1:2 %2 design matrices for TaskD+ and TaskD-
    % Construct GLM matrix 
    subrun = cell(glm_var.n_subruns,1);
    design_mat = cell(glm_var.n_subruns,1);
    SVM_mat = cell(glm_var.n_subruns,1);
    for k = 1:glm_var.n_subruns   
        %randomize order of blocks within subruns (1 is faces, 2 is houses)
        design = zeros(2,glm.subrun_length);
        SVM_design = zeros(2*glm_var.blocks,glm.subrun_length);
        subrun{k} = [repmat(1,1,glm_var.blocks),repmat(2,1,glm_var.blocks)];
        subrun{k} = subrun{k}(randperm(length(subrun{k}))); 
        %Structure is long_rest-subrun
        %Define start times for each subrun and input into design matrix
        sr1_start = glm_var.restdur+1;

        cur_timepoint = sr1_start;
        for i=1:length(subrun{k})
            design(subrun{k}(i),cur_timepoint:cur_timepoint+glm_var.blockdur-1)=1;
            cur_timepoint = cur_timepoint+glm_var.blockdur+glm_var.restdur;
        end
        %convolve with HRF to generate the design matrix
        design_mat{k}=conv2(design,glm.hrf_distri);
        % trim back to run length
        design_mat{k} = design_mat{k}(:, 1:glm.subrun_length);
        
        %Generate individual columns for SVM classification
        cur_timepoint = sr1_start;
        c1_count = 1;
        c2_count = 1;
        for i=1:length(subrun{k})
            if subrun{k}(i)==1
                SVM_design(c1_count,cur_timepoint:cur_timepoint+glm_var.blockdur-1)=1;
                c1_count = c1_count+1;
            else
                SVM_design(glm_var.blocks+c2_count,cur_timepoint:cur_timepoint+glm_var.blockdur-1)=1;
                c2_count = c2_count+1;
                
            end
            cur_timepoint = cur_timepoint+glm_var.blockdur+glm_var.restdur;
        end
        %convolve with HRF to generate the design matrix
        SVM_mat{k} = conv2(SVM_design,glm.hrf_distri);
        SVM_mat{k} = SVM_mat{k}(:,1:glm.subrun_length);
    end
    
    glm.design(j).design_mat = design_mat;
    glm.design(j).SVM_mat = SVM_mat;
    % should not be implicit
    glm.design(j).context = contexts{j};
end
