% Simulate attention and layer effects for Deming regression paper. 
% Objectives:
% 1. Simulate effects of attention in voxels for two conditions TaskD+ (dplus) and TaskD- (dminus)
% 2. Generate simulated fMRI data with superficial bias, thermal and physiological noise (measured_response) 
% 3. Apply sinusoidal, linear and mean detrending to the data
% 4. Generate scatterplots for comparison with real data to obtain estimates of noise for simulation
% 5. Demonstrate that our multiplicative gain factor mimics real data   
% 6. Utilize various metrics to recover the underlying attention modulation in response to different parameter levels
% 7. Generate a tabulated list of data to be queried by a seperate function and generate laminar profiles of different attention modulations. 


function attention_simulation(iter, pflag)

% Add SPM to path (needed for the hrf function)
if isempty(which('spm_hrf'))
    error('ensure SPM12 is on the Matlab path');
end

%number of iterations to try for each conditions
if ~exist('iter','var') || isempty(iter)
    iter = 10000;
end

% pflag = true generates scatterplot, otherwise no scatterplot.
if ~exist('pflag','var') || isempty(pflag)
    pflag = true;
end

%Reduce the n_voxels or iter to reduce computational time
% number of voxels in the ROI
n_voxels = 2500;


%Declare GLM for this study
glm_var.TR = 2.5; %TR
glm_var.blocks = 10; %no of blocks of house/face for each subrun
glm_var.blockdur = 6; %block dur in terms of TR
glm_var.restdur = 1; %rest dur between blocks in terms of TR
glm_var.n_subruns = 4; %no of subruns per attention condition (TaskD+/TaskD-)

% two noise sources at voxel level 
% Physiological noise that scales with laminar bias
physio_sigma_list = [8];
% Dimensionality of physiological noise
physio_vects = 20;
% Thermal noise that is consistent across layers
thermal_sigma_list = [15];

% strength of modulation (attentional)
attention_list = [2,3];
% strength of superficial bias
superficial_bias_list = [1, 1.5, 2];


% Make directory for output of scatterplots
if ~exist('sim_scatter', 'dir') && pflag
    mkdir('sim_scatter')
end

results(iter) = struct();

% custom seed but preserve deterministic behaviour over iterations
seed = 123456;

% if no parallel computing available, use: 
% for iterind=1:iter
parfor iterind=1:iter
    % Declaring random seed so that results are reproducible
    % Seed needs to be declared inside parfor loop otherwise matlab will seed individual workers randomly
    count = 1;
    rng(iterind * seed)
    glm = create_glm(glm_var);
    for physio_sigma = physio_sigma_list
        for thermal_sigma = thermal_sigma_list
            for superficial_bias = superficial_bias_list
                for attentional_modulation = attention_list
                    params = struct('physio_sigma',physio_sigma,'physio_vects',physio_vects,'thermal_sigma',thermal_sigma,'superficial_bias',superficial_bias,'attentional_modulation',attentional_modulation);
                    [estimates, plot_vars] = attention_simulation_iteration(n_voxels,glm_var,params,glm);
                    if iterind==1 && pflag            % only plot first iteration
                        scatter_plot(params,plot_vars.dplus,plot_vars.dminus,plot_vars.deming_regression);
                    end
                    results(iterind).params(count) = params;
                    results(iterind).estimates(count) = estimates; 
                    count = count+1;
                end
            end
        end
    end
end
save('att_sim_results.mat','results');

end

function [estimates, plot_vars] = attention_simulation_iteration(n_voxels,glm_var,params,glm)
    
    % let's suppose this example ROI has more face cells than house cells (let's call it FFA)
    % Calculating this here so that each iteration can have different frequencies
    % For Section 4.4 results, replace "n_neurons.face = abs(normrnd(0,1.1,[1, n_voxels]));" with "n_neurons.face = abs(normrnd(0,0.5,[1, n_voxels]));"
    n_neurons.face = abs(normrnd(0,1.45,[1, n_voxels]));
    n_neurons.house = abs(normrnd(0,0.7,[1, n_voxels]));
    
    
    %======================================================================
    % Task with distractor - TaskD+ (dplus)
    %======================================================================
    
    % face stimulus, attention task, with distractor
    attention.face = params.attentional_modulation;
    attention.house = 1;
    firing_rate.face = 1;
    firing_rate.house = 1;  
    dplus.face_response = calc_voxel_response(firing_rate, n_neurons, attention);

    % house stimulus, attention task, with distractor
    attention.face = 1;
    attention.house = params.attentional_modulation;
    firing_rate.face = 1;
    firing_rate.house = 1;
    dplus.house_response = calc_voxel_response(firing_rate, n_neurons, attention);
    
    % construct timecourses by scaling design matrix according to the face and house responses
    dplus.measured_response_array = calc_measured_response(dplus,glm.design(1),params);
    
    
    %detrending data
    dplus.detrended_response_array = detrend(dplus.measured_response_array);
    dplus.detrended_design_array = detrend(glm.design(1).design_mat);
    dplus.detrended_SVM_array = detrend(glm.design(1).SVM_mat);
    
    dplus.detrended_response = horzcat(dplus.detrended_response_array{:});
    dplus.detrended_design = horzcat(dplus.detrended_design_array{:});
    
    %Method 1: Z-scoring
    dplus.zbeta_estimates = zscore(dplus.detrended_response')'/dplus.detrended_design;
    estimates.zscore = mean(dplus.zbeta_estimates(:,1)-dplus.zbeta_estimates(:,2));
    
    %Method 2: SVM classifier for TaskD+
    correct_class = 0;
    for i=1:4
        train_ind = setdiff(1:4,i);
        test_ind = i;
        train_design = blkdiag(dplus.detrended_SVM_array{train_ind});
        test_design = dplus.detrended_SVM_array{test_ind};
        
        train_data = horzcat(dplus.detrended_response_array{train_ind});
        test_data = dplus.detrended_response_array{test_ind};
        
        train_betas = train_data/train_design;
        test_betas = test_data/test_design;
        
        test_labels = [repmat(0,1,glm_var.blocks),repmat(1,1,glm_var.blocks)]; 
        train_labels = [test_labels,test_labels,test_labels];
        svm_model = fitcsvm(train_betas',train_labels);
        predict_labels = predict(svm_model,test_betas');
        for s = 1:size(predict_labels)
            if predict_labels(s) == test_labels(s)
                correct_class= correct_class+1;
            end
        end
    end
    estimates.SVM = 100*correct_class/(2*glm_var.n_subruns*glm_var.blocks);
    
    %Method 3: LDC for TaskD+
    LDC_dist = nan(4,1);
    for i=1:4
        train_ind = setdiff(1:4,i);
        test_ind = i;
        train_design = horzcat(dplus.detrended_design_array{train_ind});
        test_design = dplus.detrended_design_array{test_ind};
        train_data = horzcat(dplus.detrended_response_array{train_ind});
        test_data = dplus.detrended_response_array{test_ind};
        
        train_betas = train_data/train_design;
        train_con = train_betas *[1;-1];
        train_predict = train_betas*train_design;
        train_res = train_data - train_predict;
        train_cov = covdiag(train_res');
        train_cov = train_cov./size(train_cov,1);
        weights = train_con' * pinv(train_cov);
        
        test_betas = test_data/test_design;
        test_con = test_betas * [1;-1];
        
        LDC_dist(i) = test_con'*weights';
    end
    
    estimates.LDC = mean(LDC_dist);
    
    
    %Fitting the measured response to the GLM for the other methods
    dplus.beta_estimates = dplus.detrended_response/dplus.detrended_design;
    dplus.contrast_estimates = dplus.beta_estimates(:,1)-dplus.beta_estimates(:,2);
    
    %Store tstat for comparison with real fMRI data
    model_fit = dplus.beta_estimates *dplus.detrended_design;
    dplus.residual = dplus.detrended_response - model_fit;
    dplus.std_error = std(dplus.residual,0,2)/sqrt(size(dplus.residual,2));
    dplus.tstat = dplus.contrast_estimates./dplus.std_error;
    estimates.tstat_dplus_mean = mean(dplus.tstat);
    estimates.tstat_dplus_std = std(dplus.tstat);
    
    %Store BOLD responses for reference
    dplus.real_bold_response = dplus.face_response'-dplus.house_response';
    estimates.dplus_measured_response = mean(dplus.contrast_estimates);
    

    %======================================================================
    % next, attention task without distractor - TaskD- (dminus)
    %======================================================================
    
    % face stimulus, attention task, no distractor
    attention.face = params.attentional_modulation;
    attention.house = 1;
    firing_rate.face = 1;
    firing_rate.house = 0;
    dminus.face_response = calc_voxel_response(firing_rate, n_neurons, attention);

    % house stimulus, attention task, no distractor
    attention.face = 1;
    attention.house = params.attentional_modulation;
    firing_rate.face = 0;
    firing_rate.house = 1;
    dminus.house_response = calc_voxel_response(firing_rate, n_neurons, attention);
    
    dminus.measured_response_array = calc_measured_response(dminus,glm.design(2),params);
    
    dminus.detrended_response_array = detrend(dminus.measured_response_array);
    dminus.detrended_design_array = detrend(glm.design(2).design_mat);
    
    dminus.detrended_response = horzcat(dminus.detrended_response_array{:});
    dminus.detrended_design = horzcat(dminus.detrended_design_array{:});
    
    
    %Fitting the measured response to the GLM 
    dminus.beta_estimates = dminus.detrended_response/dminus.detrended_design;
    dminus.contrast_estimates = dminus.beta_estimates(:,1)-dminus.beta_estimates(:,2);
    
    %Store tstat for comparison with real fMRI data
    model_fit = dminus.beta_estimates *dminus.detrended_design;
    dminus.residual = dminus.detrended_response - model_fit;
    dminus.std_error = std(dminus.residual,0,2)/sqrt(size(dminus.residual,2));
    dminus.tstat = dminus.contrast_estimates./dminus.std_error;
    
    estimates.tstat_dminus_mean = mean(dminus.tstat);
    estimates.tstat_dminus_std = std(dminus.tstat);


    
    
    %======================================================================
    % Evaluating Methods 4,5 and 6
    %======================================================================
    
    % Method 4: Estimating Attention Modulation with mean(dplus/dminus)
    raw_ratio = mean(dplus.contrast_estimates./dminus.contrast_estimates);
    estimates.raw_ratio =  raw_ratio;
    estimates.raw_ratio_est =  1/(1-raw_ratio);
    
    % Method 5: Estimating Attention Modulation with mean(dplus)/mean(dminus)
    ROI_ratio = mean(dplus.contrast_estimates)/mean(dminus.contrast_estimates);
    estimates.ROI_ratio =  ROI_ratio;
    estimates.ROI_ratio_est =  1/(1-ROI_ratio);
   
    % Method 6: Estimating Attention Modulation with dplus/dminus Deming Regression
    deming_regression = deming(dminus.contrast_estimates,dplus.contrast_estimates);
    estimates.deming_regression =  deming_regression(2);
    estimates.deming_est = 1/(1-deming_regression(2));
    
    
    % Assigning variables needed for plotting 
    plot_vars=struct('dplus',dplus,'dminus',dminus,'deming_regression',deming_regression);

    
end

function voxel_response = calc_voxel_response(firing_rate, n_neurons, attention)
    % The response of a voxel depends on two populations of cells: face cells and house cells.
    % Each cell population can be active or not (firing_rate - 1 or 0)
    % Each voxel varies in how many neuron cells there are of a given type (n_neurons, ie how much they contribute to the voxel)
    % Each voxel varies in how much it is modulated by the task (attention)
    % noise is added in a later step
    voxel_response = firing_rate.face .* n_neurons.face .* attention.face + firing_rate.house .* n_neurons.house .* attention.house;
end

function glm = create_glm(glm_var)
    glm.subrun_length = 2*glm_var.restdur+2*glm_var.blocks*(glm_var.blockdur+glm_var.restdur); 
    glm.hrf_distri = spm_hrf(glm_var.TR)';
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
            design_mat{k}=conv2(design,glm.hrf_distri,'same');
            
            
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
            SVM_mat{k} = conv2(SVM_design,glm.hrf_distri,'same');
        end
        
        glm.design(j).design_mat = design_mat;
        glm.design(j).SVM_mat = SVM_mat;
    end
end


function measured_response = calc_measured_response(cond,design,params)
    measured_response = cell(size(design.design_mat,1),1);
    noise_proj=normrnd(0,1,size(cond.face_response,2),params.physio_vects); % projection vector of physio noise vector on data
    for i=1:size(design.design_mat,1)
        raw_response = [cond.face_response', cond.house_response']*design.design_mat{i};
        thermal_noise = normrnd(0,1,size(raw_response))*params.thermal_sigma;
        noise_vect=normrnd(0,1,params.physio_vects,size(raw_response,2));   % physio noise vector that is projected across different voxels
        physio_rnd = noise_proj*noise_vect;
        physio_rnd = physio_rnd/sqrt(params.physio_vects/2);
        physio_noise = physio_rnd * params.physio_sigma;
        measured_response{i} = thermal_noise+params.superficial_bias*(raw_response+physio_noise);
    end
end

function detrended_data = detrend(data)
     detrended_data = cell(size(data,1),1);
    for i=1:size(data,1)
        % Linear, first order sinusoidal and mean detrending
        nvols = size(data{i},2);
        linear_trend=1:nvols;
        sin_trend=sin((linear_trend)*2*pi/(nvols-1));
        cos_trend=cos((linear_trend)*2*pi/(nvols-1));
        mean_trend = ones(1,nvols);
        dt_design = [linear_trend;sin_trend;cos_trend;mean_trend];

        trend = data{i}/dt_design;
        est = trend*dt_design;
        detrended_data{i} = data{i} - est;
    end
end

function scatter_plot(params,dplus,dminus,deming_regression)
    %Generate deming plots for comparison with real data
    fname = sprintf('sim_scatter/dplus_dminus_physio_sigma_%g_thermal_sigma_%g_bias_%g_att_%g.png',params.physio_sigma,params.thermal_sigma,params.superficial_bias,params.attentional_modulation);
    figure
    scatter(dminus.contrast_estimates,dplus.contrast_estimates)
    axis_limits = floor(max([dminus.contrast_estimates;dplus.contrast_estimates])/5)*5;
    ylim([-axis_limits axis_limits]);
    xlim([-axis_limits axis_limits]);
    ylabel('TaskD+')
    xlabel('TaskD-')
    xFit = linspace(-axis_limits ,axis_limits, 1000);
    yFit = polyval([deming_regression(2),deming_regression(1)] , xFit);
    hold on;
    plot(xFit, yFit, 'r-')
    ax = gca
    ax.FontSize = 10;
    saveas(gcf,fname,'png');
    close all
end
