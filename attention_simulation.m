% Simulate attention and layer effects for Deming regression paper. 
% Objectives:
% 1. Simulate effects of attention in voxels for two conditions TaskD+ (dplus) and TaskD- (dminus)
% 2. Generate simulated fMRI data with superficial bias, thermal and physiological noise (measured_response) 
% 3. Apply sinosodial, linear and mean detrending to the data
% 4. Generate scatterplots for comparison with real data to obtain estimates of noise for simulation
% 5. Demonstrate that our multiplicative gain factor mimics real data   
% 6. Utilize various metrics to recover the underlying attention modulation in response to different parameter levels
% 7. Generate a tabulated list of data to be queried by a seperate function and generate laminar profiles of different attention modulations. 


function attention_simulation

% Add SPM to path (needed for the hrf function)
% addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/


rng('default')

% number of voxels in the ROI
n_voxels = 500;
%number of iterations to try for each conditions
iter = 10000;

%Reduce the n_voxels or iter to reduce computational time

%Declare GLM for this study
GLM_var.TR = 2.5; %TR
GLM_var.blocks = 10; %no of blocks for each condition
GLM_var.blockdur = 6; %block dur in terms of TR
GLM_var.restdur = 1; %rest dur between blocks in terms of TR

% two noise sources at voxel level 
% Physiological noise that scales with laminar bias
physio_sigma_list = [8];
% Thermal noise that is consistent across layers
thermal_sigma_list = [12];

% strength of modulation (attentional)
attention_list = [2,3];
% strength of superficial bias
superficial_bias_list = [1, 1.5, 2];

% pflag = 1 generates scatterplot, otherwise no scatterplot.
pflag = 1; 



% Make directory for output of scatterplots
if ~exist('sim_scatter', 'dir') && pflag == 1
    mkdir('sim_scatter')
end


%loop over all variables that we are interested in and call the model function
all_res = cell(iter,1);

% if no parallel computing available, use: 
% for iterind=1:iter  
parfor iterind=1:iter     
    glm = create_glm(GLM_var);
    for noiseind = 1:numel(physio_sigma_list)
        for thermalind = 1:numel(thermal_sigma_list)
            for biasind = 1:numel(superficial_bias_list)
                for attind = 1:numel(attention_list)

                    thermal_sigma = thermal_sigma_list(thermalind);
                    physio_sigma = physio_sigma_list(noiseind);
                    superficial_bias = superficial_bias_list(biasind);
                    attention_modulation = attention_list(attind);
                    
                    estimate = attention_simulation_iteration(iterind,n_voxels,GLM_var,thermal_sigma,physio_sigma,superficial_bias,attention_modulation,glm,pflag);
                    all_res{iterind}(end+1,:) = [physio_sigma,thermal_sigma,superficial_bias,attention_modulation,estimate]; 
                end
            end
        end
    end
end
save('att_sim_results_same_pref.mat','all_res');

end

function voxel_response = calc_voxel_response(firing_rate, n_neurons, attention)
    % The response of a voxel depends on two populations of cells: face cells and house cells.
    % Each cell population can be active or not (firing_rate - 1 or 0)
    % Each voxel varies in how many neuron cells there are of a given type (n_neurons, ie how much they contribute to the voxel)
    % Each voxel varies in how much it is modulated by the task (attention)
    % noise is added in a later step
    voxel_response = firing_rate.face .* n_neurons.face .* attention.face + firing_rate.house .* n_neurons.house .* attention.house;
end

function glm = create_glm(GLM_var)
    glm.total_length = 2*GLM_var.restdur+2*GLM_var.blocks*(GLM_var.blockdur+GLM_var.restdur);
    
    for j=1:2 %2 design matrices for TaskD+ and TaskD-
        % Construct GLM matrix 
        subrun = cell(4,1);
        SVM_subrun = cell(4,1);
        run_mat = cell(4,1);
        for k = 1:4     % 4 runs per condition
            %randomize order of blocks within subruns (1 is faces, 2 is houses)
            design = zeros(2,glm.total_length);
            SVM_design = zeros(2*GLM_var.blocks,glm.total_length);
            subrun{k} = [repmat(1,1,GLM_var.blocks),repmat(2,1,GLM_var.blocks)];
            subrun{k} = subrun{k}(randperm(length(subrun{k}))); 
            %Structure is long_rest-subrun
            %Define start times for each subrun and input into design matrix
            sr1_start = GLM_var.restdur;

            cur_timepoint = sr1_start;
            for i=1:length(subrun{k})
                design(subrun{k}(i),cur_timepoint:cur_timepoint+GLM_var.blockdur-1)=1;
                cur_timepoint = cur_timepoint+GLM_var.blockdur+GLM_var.restdur;
            end
            run_mat{k}=design;
            
            
            %Generate individual columns for SVM classification
            cur_timepoint = sr1_start;
            c1_count = 1;
            c2_count = 1;
            for i=1:length(subrun{k})
                if subrun{k}(i)==1
                    SVM_design(c1_count,cur_timepoint:cur_timepoint+GLM_var.blockdur-1)=1;
                    c1_count = c1_count+1;
                else
                    SVM_design(10+c2_count,cur_timepoint:cur_timepoint+GLM_var.blockdur-1)=1;
                    c2_count = c2_count+1;
                    
                end
                cur_timepoint = cur_timepoint+GLM_var.blockdur+GLM_var.restdur;
            end
            SVM_subrun{k} = SVM_design;
        end
        design_mat = horzcat(run_mat{1:4});
        %convolve with HRF to generate the design matrix
        glm.hrf_distri = spm_hrf(GLM_var.TR)';
        conv_design = conv2(design_mat,glm.hrf_distri,'same');
        
        glm.final_design{j} = conv_design;
        
        SVM_design = blkdiag(SVM_subrun{1:4});
        conv_SVM = conv2(SVM_design,glm.hrf_distri,'same');
        
        glm.SVM_design{j} = conv_SVM;
    end
end


function estimate = attention_simulation_iteration(count,n_voxels,GLM_var,thermal_sigma,physio_sigma,superficial_bias,attention_modulation,glm,pflag)
    
    % let's suppose this example ROI has more face cells than house cells (let's call it FFA)
    % Calculating this here so that each iteration can have different frequencies
    n_neurons.face = abs(normrnd(0,2,[1, n_voxels]));
    n_neurons.house = abs(normrnd(0,1,[1, n_voxels]));
    
    %======================================================================
    % Task with distractor - TaskD+ (dplus)
    %======================================================================
    
    % face stimulus, attention task, with distractor
    attention.face = attention_modulation;
    attention.house = 1;
    firing_rate.face = 1;
    firing_rate.house = 1;  
    dplus.face_response = calc_voxel_response(firing_rate, n_neurons, attention);

    % house stimulus, attention task, with distractor
    attention.face = 1;
    attention.house = attention_modulation;
    firing_rate.face = 1;
    firing_rate.house = 1;
    dplus.house_response = calc_voxel_response(firing_rate, n_neurons, attention);
    
    dplus.raw_response = [dplus.face_response', dplus.house_response']*glm.final_design{1};
    dplus.thermal_noise = normrnd(0,1,size(dplus.raw_response))*thermal_sigma;
    dplus.physio_noise = normrnd(0,1,size(dplus.raw_response))*physio_sigma;
    dplus.mesured_response = dplus.thermal_noise+superficial_bias*(dplus.raw_response+dplus.physio_noise);
    
    % Linear, first order sinosodial and mean detrending
    nvols = size(dplus.mesured_response,2);
    linear_trend=1:nvols;
    sin_trend=sin((linear_trend)*2*pi/(nvols-1));
    cos_trend=cos((linear_trend)*2*pi/(nvols-1));
    mean_trend = ones(1,nvols);
    
    dt_design = [linear_trend;sin_trend;cos_trend;mean_trend];
    
    dplus.measured_trend = dplus.mesured_response/dt_design;
    dplus.measured_est = dplus.measured_trend*dt_design;
    dplus.detrended_response = dplus.mesured_response - dplus.measured_est;
    
    dplus.design_trend = glm.final_design{1}/dt_design;
    dplus.design_est = dplus.design_trend*dt_design;
    dplus.detrended_design = glm.final_design{1} - dplus.design_est;
    
    dplus.SVM_design_trend = glm.SVM_design{1}/dt_design;
    dplus.SVM_design_est = dplus.SVM_design_trend*dt_design;
    dplus.detrended_SVM = glm.SVM_design{1} - dplus.SVM_design_est;
    
    %Method 1: Z-scoring
    dplus.zbeta_estimates = zscore(dplus.detrended_response')'/dplus.detrended_design;
    dplus.zscore = mean(dplus.zbeta_estimates(:,1)-dplus.zbeta_estimates(:,2));
    
    %Fitting the measured response to the GLM for the other methods
    dplus.beta_estimates = dplus.detrended_response/dplus.detrended_design;
    dplus.contrast_estimates = dplus.beta_estimates(:,1)-dplus.beta_estimates(:,2);
    
    % Split data into sub runs for SVM and LDC cross-validation
    subrun_data = cell(4,1);
    SVM_mat = cell(4,1);
    LDC_mat = cell(4,1);
    for i=1:4
        subrun_data{i}=dplus.detrended_response(:,(i-1)*glm.total_length+1:i*glm.total_length);
        SVM_mat{i}=dplus.detrended_SVM((i-1)*20+1:i*20,(i-1)*glm.total_length+1:i*glm.total_length);
        LDC_mat{i}=dplus.detrended_design(:,(i-1)*glm.total_length+1:i*glm.total_length);
    end    

    %Method 2: SVM classifier for TaskD+
    correct_class = 0;
    for i=1:4
        train_ind = setdiff(1:4,i);
        test_ind = i;
        train_design = blkdiag(SVM_mat{train_ind});
        test_design = SVM_mat{test_ind};
        
        train_data = horzcat(subrun_data{train_ind});
        test_data = subrun_data{test_ind};
        
        train_betas = train_data/train_design;
        test_betas = test_data/test_design;
        
        test_labels = [repmat(0,1,GLM_var.blocks),repmat(1,1,GLM_var.blocks)]; 
        train_labels = [test_labels,test_labels,test_labels];
        svm_model = fitcsvm(train_betas',train_labels);
        predict_labels = predict(svm_model,test_betas');
        for s = 1:size(predict_labels)
            if predict_labels(s) == test_labels(s)
                correct_class= correct_class+1;
            end
        end
    end
    dplus.SVM = correct_class/8*GLM_var.blocks;
    
    %Method 3: LDC for TaskD+
    LDC_dist = nan(4,1);
    for i=1:4
        train_ind = setdiff(1:4,i);
        test_ind = i;
        train_design = horzcat(LDC_mat{train_ind});
        test_design = LDC_mat{test_ind};
        conv_train_design = conv2(train_design,glm.hrf_distri,'same');
        conv_test_design = conv2(test_design,glm.hrf_distri,'same');
        train_data = horzcat(subrun_data{train_ind});
        test_data = subrun_data{test_ind};
        
        train_betas = train_data/conv_train_design;
        train_con = train_betas *[1;-1];
        train_predict = train_betas*conv_train_design;
        train_res = train_data - train_predict;
        train_cov = covdiag(train_res');
        train_cov = train_cov./size(train_cov,1);
        weights = train_con' * pinv(train_cov);
        weights = weights/norm(weights);
        
        test_betas = test_data/conv_test_design;
        test_con = test_betas * [1;-1];
        
        LDC_dist(i) = test_con'*weights';
    end
    
    dplus.LDC = nanmean(LDC_dist);
    
    %Store real contrast for reference
    dplus.real_bold_response = dplus.face_response'-dplus.house_response';
    

    %======================================================================
    % next, attention task without distractor - TaskD- (dminus)
    %======================================================================
    
    % face stimulus, attention task, no distractor
    attention.face = attention_modulation;
    attention.house = 1;
    firing_rate.face = 1;
    firing_rate.house = 0;
    dminus.face_response = calc_voxel_response(firing_rate, n_neurons, attention);

    % house stimulus, attention task, no distractor
    attention.face = 1;
    attention.house = attention_modulation;
    firing_rate.face = 0;
    firing_rate.house = 1;
    dminus.house_response = calc_voxel_response(firing_rate, n_neurons, attention);
    
    dminus.raw_response = [dminus.face_response', dminus.house_response']*glm.final_design{2};
    dminus.thermal_noise = normrnd(0,2,size(dminus.raw_response))*thermal_sigma;
    dminus.physio_noise = normrnd(0,2,size(dminus.raw_response))*physio_sigma;
    dminus.measured_response = dminus.thermal_noise+superficial_bias*(dminus.raw_response+dminus.physio_noise);   
    
    dminus.measured_trend = dminus.measured_response/dt_design;
    dminus.measured_est = dminus.measured_trend*dt_design;
    dminus.detrended_response = dminus.measured_response - dminus.measured_est;
    
    dminus.design_trend = glm.final_design{2}/dt_design;
    dminus.design_est = dminus.design_trend*dt_design;
    dminus.detrended_design = glm.final_design{2} - dminus.design_est;
    
    %Fitting the measured response to the GLM 
    dminus.beta_estimates = dminus.detrended_response/dminus.detrended_design;
    dminus.contrast_estimates = dminus.beta_estimates(:,1)-dminus.beta_estimates(:,2);
    
    %======================================================================
    % Evaluating Methods 4,5 and 6
    %======================================================================
    
    % Method 4: Estimating Attention Modulation with mean(dplus/dminus)
    dplus_dminus_ratio = mean(dplus.contrast_estimates./dminus.contrast_estimates);
    dplus_dminus_ratio_est =  1/(1-dplus_dminus_ratio);
    
    % Method 5: Estimating Attention Modulation with mean(dplus)/mean(dminus)
    dplus_dminus_ROI_ratio = mean(dplus.contrast_estimates)/mean(dminus.contrast_estimates);
    dplus_dminus_ROI_ratio_est =  1/(1-dplus_dminus_ROI_ratio);
   
    % Method 6: Estimating Attention Modulation with dplus/dminus Deming Regression
    deming_dplus_dminus = deming(dminus.contrast_estimates,dplus.contrast_estimates);
    deming_dplus_dminus_est = 1/(1-deming_dplus_dminus(2));
    
    estimate = [dplus.zscore, dplus.SVM, dplus.LDC, dplus_dminus_ratio_est, dplus_dminus_ROI_ratio_est, deming_dplus_dminus_est, mean(dplus.real_bold_response), mean(dplus.contrast_estimates)];
    
    %Generate deming plots for comparison with real data
    fname = sprintf('sim_scatter/dplus_dminus_physio_sigma_%g_thermal_sigma_%g_salient_%g_bias_%g_att_%g.png',physio_sigma,thermal_sigma,1,superficial_bias,attention_modulation);
    if count==1 && pflag
        figure
        scatter(dminus.contrast_estimates,dplus.contrast_estimates)
        axis_limits = floor(max([dminus.contrast_estimates;dplus.contrast_estimates]))+1;
        ylim([-axis_limits axis_limits]);
        xlim([-axis_limits axis_limits]);
        ylabel('TaskD+')
        xlabel('TaskD-')
        xFit = linspace(-axis_limits ,axis_limits, 1000);
        yFit = polyval([deming_dplus_dminus(2),deming_dplus_dminus(1)] , xFit);
        hold on;
        plot(xFit, yFit, 'r-')
        saveas(gcf,fname,'png');
        close all
    end  
    
end
