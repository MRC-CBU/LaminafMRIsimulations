% Simulate attention and layer effects for 7T project. 
% Objectives:
% 1. Simulate effects of attention in voxels for two conditions TaskD+ and TaskD-
% 2. Generate simulated fMRI data with superficial bias, thermal and
% physiological noise (measured_response) 
% 3. Apply sinosodial, linear and mean detrending to the data
% 4. Generate scatterplots for comparison with real data to obtain
% estimates of noise for simulation
% 5. Demonstrate that our multiplicative gain factor mimics real data   
% 6. Utilize our potential metrics to recover the underlying attention
% modulation in response to different parameter levels with 500 iterations to generate a stable
% estimate of mean and variance
% 7. Generate a tabulated list of data to be queried by a seperate function
% and generate laminar profiles of different attention modulations. 

function attention_simulation_new

% The response of a voxel depends on two populations of cells: face cells and house
% cells.
% Each cell population can be active or not (firingrate - 1 or 0)
% Each voxel varies in how many cells there are of a given type (frequency, ie how much they contribute to the voxel -
% equivalently you can think of this as response gain in equally-sized populations)
% Each voxel varies in how much it is modulated by the task (attention)
% noise is added in a later step
voxelsignal = @ (firingrate, frequency, attention) ...
    firingrate.face .* frequency.face_cells .* attention.face + ...
    firingrate.house .* frequency.house_cells .* attention.house;


% based on this we can simulate the response of an ROI to our task
nvoxel = 500;
%number of iterations to try for each conditions
iter = 10000;

%Reduce the nvoxel or iter to reduce computational time

%Declare GLM for this study
GLM_var.TR = 2.5; %TR
GLM_var.blocks = 10; %no of blocks for each condition
GLM_var.blockdur = 6; %block dur in terms of TR
GLM_var.restdur = 1; %rest dur between blocks in terms of TR



% two noise sources at voxel level
% Physiological noise that scales with laminar bias
noiselevel.physio = [4];
% Thermal noise that is consistent across layers
noiselevel.thermal = [6];

% strength of modulation (attentional)
attentional_modulation = [2,3];
% strength of superficial bias
superficial_bias = [1, 1.5, 2];

%Declaring Results and Variance Table 
reps=numel(noiselevel.physio)*numel(attentional_modulation)*numel(superficial_bias)*numel(noiselevel.thermal);
cols=4+8;
sz = [reps,cols];
var_types = repmat({'double'},1,cols);
var_names = {'Physio_Noise','Thermal_Noise','Superficial_Bias','Attention_Modulation','Zscore_TaskD','SVM_TaskD','LDC_TaskD','Mean_TaskD_TaskND','Mean_ROI_TaskD_TaskND','Deming_TaskD_TaskND','Real_BOLD_TaskD','Measured_BOLD_TaskD'};

res_table = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names);
var_table = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names);


% Make directory for output of scatterplots
if ~exist('sim_scatter', 'dir')
    mkdir('sim_scatter')
end


%loop over all variables that we are interested in and call the model function
cur_row=1;
for noiseind = 1:numel(noiselevel.physio)
    for thermalind = 1:numel(noiselevel.thermal)
        for biasind = 1:numel(superficial_bias)
            for attind = 1:numel(attentional_modulation)
                % Declare matrix for current results and iteration
                % variables for parfor loop
                cur_res=NaN(iter,8);
                iter_tnoise = noiselevel.thermal(thermalind);
                iter_pnoise = noiselevel.physio(noiseind);
                iter_bias = superficial_bias(biasind);
                iter_att = attentional_modulation(attind);
                parfor iterind=1:iter
                    cur_res(iterind,:)=attention_simulation_iteration(iterind,nvoxel,GLM_var,voxelsignal,iter_tnoise,iter_pnoise,iter_bias,iter_att);
                end
                res_table(cur_row,:)=[{noiselevel.physio(noiseind),noiselevel.thermal(thermalind),superficial_bias(biasind),attentional_modulation(attind)},num2cell(nanmean(cur_res,1))];
                var_table(cur_row,:)=[{noiselevel.physio(noiseind),noiselevel.thermal(thermalind),superficial_bias(biasind),attentional_modulation(attind)},num2cell(nanstd(cur_res,0,1))];
                med_table(cur_row,:)=[{noiselevel.physio(noiseind),noiselevel.thermal(thermalind),superficial_bias(biasind),attentional_modulation(attind)},num2cell(nanmedian(cur_res,1))];
                p25_table(cur_row,:)=[{noiselevel.physio(noiseind),noiselevel.thermal(thermalind),superficial_bias(biasind),attentional_modulation(attind)},num2cell(prctile(cur_res,25))];
                p75_table(cur_row,:)=[{noiselevel.physio(noiseind),noiselevel.thermal(thermalind),superficial_bias(biasind),attentional_modulation(attind)},num2cell(prctile(cur_res,75))];
                cur_row=cur_row+1;
            end
        end
    end
end
save('att_sim_results_same_pref.mat','res_table','var_table','med_table','p25_table','p75_table');

end

function cur_est = attention_simulation_iteration(count,nvoxel,GLM_var,voxelsignal,tnoise,pnoise,cur_bias,cur_att)
    
    final_design=cell(2,1);
    for j=1:2 %2 design matrices for TaskD+ and TaskD- 
        % Construct GLM matrix 
        for k = 1:4     % 4 runs per condition
            %randomize order of blocks within subruns (1 is faces, 2 is houses)
            total_length = 2*GLM_var.restdur+2*GLM_var.blocks*(GLM_var.blockdur+GLM_var.restdur);
            design = zeros(2,total_length);
            SVM_design = zeros(2*GLM_var.blocks,total_length);
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
        hrf_distri = spm_hrf(GLM_var.TR)';
        conv_design = conv2(design_mat,hrf_distri,'same');
        
        final_design{j} = conv_design;
        
        SVM_design = blkdiag(SVM_subrun{1:4});
        conv_SVM = conv2(SVM_design,hrf_distri,'same');
        
        SVM_design_mat{j} = conv_SVM;
    end


    % let's suppose this example ROI has more face cells than house cells (let's call it
    % FFA)
    % Calculating this here so that each iteration can have different
    % frequencies
    frequency.face_cells = abs(normrnd(0,2,[1, nvoxel]));
    frequency.house_cells = abs(normrnd(0,1,[1, nvoxel]));

    
    %======================================================================
    % Task with distractor - TaskD+
    %======================================================================
    
    % face stimulus, attention task, with distractor
    attention.face = cur_att;
    attention.house = 1;
    firingrate.face = 1;
    firingrate.house = 1;  
    face_task_dplus_response = voxelsignal(firingrate, frequency, attention);

    % house stimulus, attention task, with distractor
    attention.face = 1;
    attention.house = cur_att;
    firingrate.face = 1;
    firingrate.house = 1;
    house_task_dplus_response = voxelsignal(firingrate, frequency, attention);
    
    raw_response = [face_task_dplus_response', house_task_dplus_response']*final_design{1};
    thermal_noise = normrnd(0,2,size(raw_response))*tnoise;
    physio_noise = normrnd(0,2,size(raw_response))*pnoise;
    measured_response = thermal_noise+cur_bias*(raw_response+physio_noise);
    
    % Linear, first order sinosodial and mean detrending
    nvols = size(measured_response,2);
    linear_trend=[1:nvols];
    sin_trend=sin((linear_trend)*2*pi/(nvols-1));
    cos_trend=cos((linear_trend)*2*pi/(nvols-1));
    mean_trend = ones(1,nvols);
    
    dt_design = [linear_trend;sin_trend;cos_trend;mean_trend];
    
    measured_trend = measured_response/dt_design;
    measured_est = measured_trend*dt_design;
    detrended_response = measured_response - measured_est;
    
    design_trend = final_design{1}/dt_design;
    design_est = design_trend*dt_design;
    detrended_design = final_design{1} - design_est;
    
    SVM_design_trend = SVM_design_mat{1}/dt_design;
    SVM_design_est = SVM_design_trend*dt_design;
    detrended_SVM = SVM_design_mat{1} - SVM_design_est;
    
    %Method 1: Z-scoring
    zbeta_estimates = zscore(detrended_response')'/detrended_design;
    zscore_taskD = mean(zbeta_estimates(:,1)-zbeta_estimates(:,2));
    
    %Fitting the measured response to the GLM for the other methods
    beta_estimates = detrended_response/detrended_design;
    contrast_estimates = beta_estimates(:,1)-beta_estimates(:,2);
    
    % Split data into sub runs for SVM and LDC cross-validation
    for i=1:4
        subrun_data{i}=detrended_response(:,(i-1)*total_length+1:i*total_length);
        SVM_mat{i}=detrended_SVM((i-1)*20+1:i*20,(i-1)*total_length+1:i*total_length);
        LDC_mat{i}=detrended_design(:,(i-1)*total_length+1:i*total_length);
    end    
    
    %Method 2: SVM classifier
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
    SVM_acc = correct_class/8*GLM_var.blocks;
    
    %Method 3: LDC
    for i=1:4
        train_ind = setdiff(1:4,i);
        test_ind = i;
        train_design = horzcat(LDC_mat{train_ind});
        test_design = LDC_mat{test_ind};
        conv_train_design = conv2(train_design,hrf_distri,'same');
        conv_test_design = conv2(test_design,hrf_distri,'same');
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
    
    LDC_res = mean(LDC_dist);
    
    %Store contrast estimate for Methods 4,5 and 6
    measured_bold_task_dplus = contrast_estimates;
    
    %Store real contrast for reference
    real_bold_task_dplus = face_task_dplus_response'-house_task_dplus_response';
    

    %======================================================================
    % next, attention task without distractor - TaskD-
    %======================================================================
    
    % face stimulus, attention task, no distractor
    attention.face = cur_att;
    attention.house = 1;
    firingrate.face = 1;
    firingrate.house = 0;
    face_task_dminus_response = voxelsignal(firingrate, frequency, attention);

    % house stimulus, attention task, no distractor
    attention.face = 1;
    attention.house = cur_att;
    firingrate.face = 0;
    firingrate.house = 1;
    house_task_dminus_response = voxelsignal(firingrate, frequency, attention);
    
    raw_response = [face_task_dminus_response', house_task_dminus_response']*final_design{2};
    thermal_noise = normrnd(0,2,size(raw_response))*tnoise;
    physio_noise = normrnd(0,2,size(raw_response))*pnoise;
    measured_response = thermal_noise+cur_bias*(raw_response+physio_noise);
    
    
    measured_trend = measured_response/dt_design;
    measured_est = measured_trend*dt_design;
    detrended_response = measured_response - measured_est;
    
    design_trend = final_design{2}/dt_design;
    design_est = design_trend*dt_design;
    detrended_design = final_design{2} - design_est;
    
    %Fitting the measured response to the GLM 
    beta_estimates = detrended_response/detrended_design;
    contrast_estimates = beta_estimates(:,1)-beta_estimates(:,2);
    
    %Store contrast estimate for Methods 4,5 and 6
    measured_bold_task_dminus = contrast_estimates;
    
    %======================================================================
    % Evaluating Methods 4,5 and 6
    %======================================================================
    
    % Method 4: Estimating Attention Modulation with mean(TaskD/TaskND)
    taskD_taskND_ratio = mean(measured_bold_task_dplus./measured_bold_task_dminus);
    taskD_taskND_ratio_est =  1/(1-taskD_taskND_ratio);
    
    % Method 5: Estimating Attention Modulation with mean(TaskD)/mean(TaskND)
    taskD_taskND_ROI_ratio = mean(measured_bold_task_dplus)/mean(measured_bold_task_dminus);
    taskD_taskND_ROI_ratio_est =  1/(1-taskD_taskND_ROI_ratio);
   
    % Method 6: Estimating Attention Modulation with TaskD/TaskND Deming Regression
    slope_taskD_taskND = deming(measured_bold_task_dminus,measured_bold_task_dplus);
    taskD_taskND = slope_taskD_taskND(2);
    taskD_taskND_att_est = 1/(1-taskD_taskND);
    
    cur_est = [zscore_taskD, SVM_acc,LDC_res,taskD_taskND_ratio_est,taskD_taskND_ROI_ratio_est, taskD_taskND_att_est,mean(real_bold_task_dplus),mean(measured_bold_task_dplus)];
    
    %Generate deming plots for comparison with real data
    fname = sprintf('sim_scatter/TaskD_TaskND_pnoise_%g_tnoise_%g_salient_%g_bias_%g_att_%g.png',pnoise,tnoise,1,cur_bias,cur_att);
    if count==1
        figure
        scatter(measured_bold_task_dminus,measured_bold_task_dplus)
        axis_limits = floor(max([measured_bold_task_dminus;measured_bold_task_dplus]))+1;
        ylim([-axis_limits axis_limits]);
        xlim([-axis_limits axis_limits]);
        ylabel('TaskD+')
        xlabel('TaskD-')
        xFit = linspace(-axis_limits ,axis_limits, 1000);
        yFit = polyval([taskD_taskND,slope_taskD_taskND(1)] , xFit);
        hold on;
        plot(xFit, yFit, 'r-')
        saveas(gcf,fname,'png');
        close all
    end
    
    
    
end
