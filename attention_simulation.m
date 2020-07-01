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
sim_par.TR = 2.5; %TR
sim_par.blocks = 10; %no of blocks of house/face for each subrun
sim_par.blockdur = 6; %block dur in terms of TR
sim_par.restdur = 1; %rest dur between blocks in terms of TR
sim_par.n_subruns = 4; %no of subruns per attention condition (TaskD+/TaskD-)
% For Section 4.4 results, set both to 0.7 
sim_par.n_neurons_face_sigma = 1.45; % sigma of half-Gaussian controlling face cell frequency per voxel
sim_par.n_neurons_house_sigma = 0.7; % same for house cell frequency

% two noise sources at voxel level 
% Physiological noise that scales with laminar bias
physio_sigma_list = [11];
% Dimensionality of physiological noise
physio_vects = 20; % why not 10 as in HBM sim?
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
    glm = create_glm(sim_par);
    for physio_sigma = physio_sigma_list
        for thermal_sigma = thermal_sigma_list
            for superficial_bias = superficial_bias_list
                for attentional_modulation = attention_list
                    iter_par = struct('physio_sigma',physio_sigma,'physio_vects',physio_vects,'thermal_sigma',thermal_sigma,'superficial_bias',superficial_bias,'attentional_modulation',attentional_modulation);
                    [estimates, plot_vars] = attention_simulation_iteration(n_voxels,sim_par,iter_par,glm);
                    if iterind==1 && pflag            % only plot first iteration
                        scatter_plot(iter_par,plot_vars.dplus,plot_vars.dminus,plot_vars.deming_regression);
                    end
                    results(iterind).params(count) = iter_par;
                    results(iterind).estimates(count) = estimates; 
                    count = count+1;
                end
            end
        end
    end
end
save('att_sim_results.mat','results');

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
    ax = gca;
    ax.FontSize = 10;
    saveas(gcf,fname,'png');
    close all
end
