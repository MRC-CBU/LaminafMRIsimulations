function attention_simulation_plots(res_file)
% Plot results from the output of attention_simulation.m
% This function queries the output for the parameters we are interested in
% and generate plots comparing the ground truth against the various
% estimates using the different metrics.
% 
% If the input res is undefined, we load the sample result from
% sample_results/att_sim_results.mat
%
% attention_simulation_plots(res_file)

if ~exist('res_file', 'var') || isempty(res_file)
    res_file = fullfile(fileparts(mfilename('fullpath')),'sample_results', 'att_sim_results.mat');
end

% define which values we are interested in
% Change the values here to query different results
physio_sigma_list = repmat(11,1,3);
thermal_sigma_list = repmat(15,1,3);
superficial_bias = [2, 1.5, 1];
attentional_modulation = [3 2 3];  

cmap = [8 81 156;
    107 174 214; 
    189 215 231];
cmap = cmap/255;
    
res = load(res_file);  
res = res.results;

%Extract the 3 rows
params = [res.params];
estimates = [res.estimates];
var_names=fieldnames(estimates);
for i = 1:3
    match = [params.physio_sigma]==physio_sigma_list(i) & [params.thermal_sigma]==thermal_sigma_list(i) & [params.superficial_bias]==superficial_bias(i) & [params.attentional_modulation]==attentional_modulation(i);
    
    cur_estimates = estimates(match==1);
    
    for j=1:size(var_names,1)
        median_table.(var_names{j})(i)=median([cur_estimates.(var_names{j})]);
        p25_table.(var_names{j})(i)=prctile([cur_estimates.(var_names{j})],25);
        p75_table.(var_names{j})(i)=prctile([cur_estimates.(var_names{j})],75);
    end
end
   
% Split into multiple plots:
% Att Modulation
% 1. Deming/Mean/MeanROI
% 2. Z-score
% 3. SVM
% 4. LDC
% 5. L2 normalization


%--------------------------------------------------------------------------------------------
%Plotting attentional modulation
figure
att_plot=bar(attentional_modulation);
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'Ground Truth'});
set(gca,'XTickLabelRotation',20);
ylim([0 4]);
ylabel({'Attentional modulation', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=185;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_att.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off


%--------------------------------------------------------------------------------------------
%Plotting measured response
figure
att_plot=bar(median_table.dplus_measured_response);
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end

%Add error bars
ngroups = 3;
nbars = 1;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, median_table.dplus_measured_response(i,:), median_table.dplus_measured_response(i,:)-p25_table.dplus_measured_response(i,:),p75_table.dplus_measured_response(i,:)-median_table.dplus_measured_response(i,:), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'Simulated Response'});
set(gca,'XTickLabelRotation',20);
% ylim([0 2.5]);
ylabel({'Attentional modulation', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=185;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_measured.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------
%Plotting ratio approaches
med_mat=[median_table.deming_regression; median_table.raw_ratio; median_table.ROI_ratio];
p25_mat=[p25_table.deming_regression; p25_table.raw_ratio; p25_table.ROI_ratio];
p75_mat=[p75_table.deming_regression; p75_table.raw_ratio; p75_table.ROI_ratio];

figure
att_plot=bar(med_mat);
hold on
for b=1:3
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(med_mat, 1);
nbars = size(med_mat, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, med_mat(:,i), med_mat(:,i)-p25_mat(:,i),p75_mat(:,i)-med_mat(:,i), 'k.');
end

%Tidying up the plot and adding labels
xticks([1 2 3])
set(gca, 'XTickLabel', {'Deming Regression','Voxel Ratio', 'ROI Ratio'});
set(gca,'XTickLabelRotation',20);
% ylim([-0.1 1.2]);
ylabel({'Selectivity (A.U.)', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=470;
height=550;
set(gcf,'position',[x0,y0,width,height])
legend('Superficial','Middle','Deep','location','northwest')

%save figure
fname = sprintf('att_%g_%g_%g_plot_1.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off


%--------------------------------------------------------------------------------------------
%Plotting z-score
figure
att_plot=bar(median_table.zscore);
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
ngroups = 3;
nbars = 1;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, median_table.zscore(i,:), median_table.zscore(i,:)-p25_table.zscore(i,:),p75_table.zscore(i,:)-median_table.zscore(i,:), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'Z-scoring'});
set(gca,'XTickLabelRotation',20);
%ylim([0 10]);
ylabel({'Region-mean contrast estimate (A.U.)', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=190;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_2.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------
%Plotting SVM
figure
att_plot=bar(median_table.SVM);
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
ngroups = 3;
nbars = 1;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, median_table.SVM(i,:), median_table.SVM(i,:)-p25_table.SVM(i,:),p75_table.SVM(i,:)-median_table.SVM(i,:), 'k.');
end



%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'SVM Classification'});
set(gca,'XTickLabelRotation',20);
ylim([50 100]);
ylabel({'Classification accuracy (%)', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=185;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_3.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------
%Plotting LDC
figure
att_plot=bar(median_table.LDC);
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
ngroups = 3;
nbars = 1;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, median_table.LDC(i,:), median_table.LDC(i,:)-p25_table.LDC(i,:),p75_table.LDC(i,:)-median_table.LDC(i,:), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'LDC'});
set(gca,'XTickLabelRotation',20);
%ylim([60 100]);
ylabel({'Linear discriminant contrast (A.U.)', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=180;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_4.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------
%Plotting l2_dplus
figure
att_plot=bar(median_table.l2_dplus);
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
ngroups = 3;
nbars = 1;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, median_table.l2_dplus(i,:), median_table.l2_dplus(i,:)-p25_table.l2_dplus(i,:),p75_table.l2_dplus(i,:)-median_table.l2_dplus(i,:), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'L2 norm'});
set(gca,'XTickLabelRotation',20);
%ylim([60 100]);
ylabel({'Region-mean contrast estimate (A.U.)', '(median \pm25 percentiles)'});
x0=10;
y0=10;
width=180;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_5.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------


%--------------------------------------------------------------------------------------------

close all
end
