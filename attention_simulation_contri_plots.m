function attention_simulation_contri_plots(res_file)
% Plot results from the output of attention_simulation.m
% This function queries the output for the parameters we are interested in
% and generate plots comparing the ground truth against the various
% estimates using the different metrics.
%
% If the input res is undefined, we load the sample result from
% sample_results/att_sim_results.mat
%
% attention_simulation_contri_plots(res_file)

if ~exist('res_file', 'var') || isempty(res_file)
    res_file = fullfile(fileparts(mfilename('fullpath')),'sample_results', 'att_sim_results.mat');
end

% define which values we are interested in
% Change the values here to query different results
physio_sigma_list = repmat(11,1,3);
thermal_sigma_list = repmat(15,1,3);
superficial_bias = [2, 1.5, 1];
attentional_modulation = [3 2 3];  

cmap = [251,180,174;
    179,205,227];
cmap = cmap/255;

res = load(res_file);  
res = res.results;

%Find the indices for the 3 rows
match = NaN(3,size(res(1).params,2));
for i=1:3
    match(i,:) = [res(1).params.physio_sigma]==physio_sigma_list(i) & [res(1).params.thermal_sigma]==thermal_sigma_list(i) & [res(1).params.superficial_bias]==superficial_bias(i) & [res(1).params.attentional_modulation]==attentional_modulation(i);
end

var_names=fieldnames(res(1).estimates);
res_initial = struct();
for iter=1:size(res,2)
    estimates = res(iter).estimates;
    for i=1:3
        cur_est = estimates(match(i,:)==1);
        for j=1:size(var_names,1)
            if size([estimates.(var_names{j})],1)==1 %skipping real bold response and measured bold response because they are nvox*iter matrices
                res_initial.(var_names{j})(iter,i) = cur_est.(var_names{j});
            end
        end
    end
end

niter = size(res_initial.deming_regression,1);
%Order res_table in plotting order
res_table = struct('ground_truth',attentional_modulation,'deming_regression',res_initial.deming_regression,'raw_ratio',res_initial.raw_ratio,'ROI_ratio',res_initial.ROI_ratio,'zscore',res_initial.zscore,'SVM',res_initial.SVM,'LDC',res_initial.LDC, 'l2_dplus', res_initial.l2_dplus);
var_names=fieldnames(res_table);

contri_vects =[1 0 -1; 0.5 -1 0.5]';
assert(corr(contri_vects(:,1), contri_vects(:,2)) == 0, ...
    'the multiple regression correlation algorithm below is only valid for orthogonal regressors');
contri_tiles = NaN([size(var_names, 1), 2, 3]);

for i=1:size(var_names,1)
    % beta for zscored X and Y == pearson r, if regressors are orthogonal
    contri = zscore(contri_vects) \ zscore(res_table.(var_names{i})');
    contri_tiles(i,:,:) = prctile(contri, [25, 50, 75], 2);
end

contri_median = contri_tiles(:,:,2);
contri_lower = contri_median - contri_tiles(:,:,1);
contri_upper = contri_tiles(:,:,3) - contri_median;

%Plot the data
figure
% medians
att_plot=bar(contri_median);
hold on
for b=1:2
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(contri_median, 1);
nbars = size(contri_median, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, contri_median(:,i), contri_lower(:,i), contri_upper(:,i), 'k.');
end

%Tidying up the plot and adding labels
set(gca, 'XTickLabel', {'Ground Truth','Deming Regression','Voxel Ratio', 'ROI Ratio', 'Z-scoring', 'SVM classification', 'LDC', 'L2 norm'});
set(gca,'XTickLabelRotation',20);
ylabel({'Laminar contributions', '(median r \pm 25 percentiles)'});
x0=10;
y0=10;
width=950;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend('Superficial bias', 'Attentional modulation','location','southeast')

%save figure
fname = sprintf('att_%g_%g_%g_contri_plot.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off
