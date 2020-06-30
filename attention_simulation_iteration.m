function [estimates, plot_vars] = attention_simulation_iteration(n_voxels,sim_par,iter_par,glm)
    
    % let's suppose this example ROI has more face cells than house cells (let's call it FFA)
    % Calculating this here so that each iteration can have different frequencies
    n_neurons.face = abs(normrnd(0,sim_par.n_neurons_face_sigma,[1, n_voxels]));
    n_neurons.house = abs(normrnd(0,sim_par.n_neurons_house_sigma,[1, n_voxels]));
    
    
    %======================================================================
    % Task with distractor - TaskD+ (dplus)
    %======================================================================
    
    % face stimulus, attention task, with distractor
    dplus.face_response = calc_voxel_response(struct('face', 1, 'house', 1), ...
        n_neurons, struct('face', iter_par.attentional_modulation, 'house', 1));

    % house stimulus, attention task, with distractor
    dplus.house_response = calc_voxel_response(struct('face', 1, 'house', 1), ...
        n_neurons, struct('face', 1, 'house', iter_par.attentional_modulation));
    
    % construct timecourses by scaling design matrix from dplus context (1) according to
    % the face and house responses
    dplus.measured_response_array = calc_measured_response(dplus,glm.design(1),iter_par);

    dplus = calc_fmri_model(dplus, glm.design(1));
    
    %Method 1: Z-scoring
    dplus.zbeta_estimates = zscore(dplus.detrended_response')'/dplus.detrended_design;
    estimates.zscore = mean(dplus.zbeta_estimates(:,1)-dplus.zbeta_estimates(:,2));
    
    %Method 2: SVM classifier for TaskD+
    svm_correct = {};
    for test_ind=1:numel(dplus.detrended_SVM_array)
        train_ind = setdiff(1:numel(dplus.detrended_SVM_array),test_ind);
        train_design = blkdiag(dplus.detrended_SVM_array{train_ind});
        test_design = dplus.detrended_SVM_array{test_ind};
        
        train_data = horzcat(dplus.detrended_response_array{train_ind});
        test_data = dplus.detrended_response_array{test_ind};
        
        train_betas = train_data/train_design;
        test_betas = test_data/test_design;
        
        test_labels = [zeros(1,sim_par.blocks),ones(1,sim_par.blocks)]; 
        train_labels = repmat(test_labels, [1, numel(train_ind)]);
        svm_model = fitcsvm(train_betas',train_labels);
        predict_labels = predict(svm_model,test_betas');
        svm_correct{test_ind} = predict_labels' == test_labels;
    end
    estimates.SVM = 100 * mean(horzcat(svm_correct{:}));
    
    %Method 3: LDC for TaskD+
    LDC_dist = nan(sim_par.n_subruns,1);
    for test_ind=1:sim_par.n_subruns
        train_ind = setdiff(1:sim_par.n_subruns,test_ind);
        train_design = horzcat(dplus.detrended_design_array{train_ind});
        test_design = dplus.detrended_design_array{test_ind};
        train_data = horzcat(dplus.detrended_response_array{train_ind});
        test_data = dplus.detrended_response_array{test_ind};
        
        train_betas = train_data/train_design;
        train_con = train_betas *[1;-1];
        train_predict = train_betas*train_design;
        train_res = train_data - train_predict;
        train_cov = covdiag(train_res');
        % what does this do?
        train_cov = train_cov./size(train_cov,1);
        weights = train_con' / train_cov;
        
        test_betas = test_data/test_design;
        test_con = test_betas * [1;-1];
        
        LDC_dist(test_ind) = weights * test_con;
    end
    estimates.LDC = mean(LDC_dist);
    
    estimates.tstat_dplus_mean = mean(dplus.tstat);
    estimates.tstat_dplus_std = std(dplus.tstat);
    estimates.dplus_measured_response = mean(dplus.contrast_estimates);

    %======================================================================
    % next, attention task without distractor - TaskD- (dminus)
    %======================================================================
    
    % face stimulus, attention task, no distractor
    dminus.face_response = calc_voxel_response(struct('face', 1, 'house', 0), ...
        n_neurons, struct('face', iter_par.attentional_modulation, 'house', 1));

    % house stimulus, attention task, no distractor
    dminus.house_response = calc_voxel_response(struct('face', 0, 'house', 1), ...
        n_neurons, struct('face', 1, 'house', iter_par.attentional_modulation));
    
    % nb glm.design(2) is the design matrix for the dminus context (independently
    % randomised block sequence)
    dminus.measured_response_array = calc_measured_response(dminus,glm.design(2),iter_par);

    dminus = calc_fmri_model(dminus, glm.design(2));
    
    estimates.tstat_dminus_mean = mean(dminus.tstat);
    estimates.tstat_dminus_std = std(dminus.tstat);
    
    %======================================================================
    % Evaluating Methods 4,5 and 6
    %======================================================================
    
    % Method 4: Estimating Attention Modulation with mean(dplus/dminus)
    estimates.raw_ratio = mean(dplus.contrast_estimates./dminus.contrast_estimates);
    estimates.raw_ratio_est =  1/(1-estimates.raw_ratio);
    
    % Method 5: Estimating Attention Modulation with mean(dplus)/mean(dminus)
    estimates.ROI_ratio = mean(dplus.contrast_estimates)/mean(dminus.contrast_estimates);
    estimates.ROI_ratio_est =  1/(1-estimates.ROI_ratio);
   
    % Method 6: Estimating Attention Modulation with dplus/dminus Deming Regression
    deming_regression = deming(dminus.contrast_estimates,dplus.contrast_estimates);
    estimates.deming_regression =  deming_regression(2);
    estimates.deming_est = 1/(1-estimates.deming_regression);
    
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

function measured_response = calc_measured_response(condition,design,params)
    measured_response = cell(size(design.design_mat,1),1);
    noise_proj=normrnd(0,1,size(condition.face_response,2),params.physio_vects); % projection vector of physio noise vector on data
    for i=1:size(design.design_mat,1)
        raw_response = [condition.face_response', condition.house_response']*design.design_mat{i};
        thermal_noise = normrnd(0,1,size(raw_response))*params.thermal_sigma;
        noise_vect=normrnd(0,1,params.physio_vects,size(raw_response,2));   % physio noise vector that is projected across different voxels
        physio_noise = noise_proj*noise_vect;
        physio_noise = physio_noise/sqrt(params.physio_vects);
        physio_noise = physio_noise * params.physio_sigma;
        measured_response{i} = thermal_noise+params.superficial_bias*(raw_response+physio_noise);
    end
end

function d = calc_fmri_model(d, design)

%detrending data
d.detrended_response_array = detrend_timecourses(d.measured_response_array);
d.detrended_design_array = detrend_timecourses(design.design_mat);
d.detrended_SVM_array = detrend_timecourses(design.SVM_mat);

d.detrended_response = horzcat(d.detrended_response_array{:});
d.detrended_design = horzcat(d.detrended_design_array{:});

%Fitting the measured response to the GLM for the other methods
d.beta_estimates = d.detrended_response/d.detrended_design;
d.contrast_estimates = d.beta_estimates(:,1)-d.beta_estimates(:,2);

%Store tstat for comparison with real fMRI data
model_fit = d.beta_estimates * d.detrended_design;
residual = d.detrended_response - model_fit;
df = size(d.residual, 2) - rank(d.detrended_design');
mrss = sum(d.residual.^2, 2) / df;
cmat = inv(d.detrended_design * d.detrended_design');
convec = [1 -1];
se = sqrt(diag(convec * cmat * convec') * mrss);
d.tstat = d.contrast_estimates ./ se;

%Store BOLD responses for reference
d.real_bold_response = d.face_response'-d.house_response';

end
