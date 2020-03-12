function attention_simulation_contri_plots()
% Plotting results from the output of attention_simulation.m
% This function queries the output for the parameters we are interested in
% and generate plots comparing the ground truth against the various
% estimates using the different metrics.

% define which values we are interested in
% Change the values here to query different results, e.g.
% int_attention=[1.5 1 1.5] for retinotopic attention plots
int_physio_noise = repmat(4,1,3);              
int_thermal_noise = repmat(6,1,3);     
int_attention = [3 2 3];                     
int_bias = [2, 1.5, 1];               





cmap = [251,180,174;
    179,205,227];
cmap = cmap/255;


%load res_table and var_table
res = load('att_sim_results.mat'); 
res_table = res.res_table;
var_table = res.var_table;
p25_table = res.p25_table;
p75_table = res.p75_table;

%Pre-define interested results+variance table
int_res = res_table(1:3,:);
int_res{:,:}=NaN;
int_var = int_res;
int_p25 = int_res;
int_p75 = int_res;

%Extract the 3 rows
for i=1:3
    selected_row=all([res_table.Physio_Noise==int_physio_noise(i),res_table.Thermal_Noise==int_thermal_noise(i),res_table.Superficial_Bias==int_bias(i),res_table.Attention_Modulation==int_attention(i)],2);
    int_res(i,:)=res_table(selected_row,:);
    int_var(i,:)=var_table(selected_row,:);
    int_p25(i,:)=p25_table(selected_row,:);
    int_p75(i,:)=p75_table(selected_row,:);    
end
   
% Split into multiple plots:
% 1. Att Modulation/Deming/Mean/MeanROI
% 2. SVM
% 3. LDC
% 4. Z-score

%Combine the data into a results and error matrix
res_mat=[int_res.Attention_Modulation, int_res.Deming_TaskD_TaskND, int_res.Mean_TaskD_TaskND, int_res.Mean_ROI_TaskD_TaskND, int_res.Zscore_TaskD, int_res.SVM_TaskD, int_res.LDC_TaskD]';
var_mat=[zeros(3,1), int_var.Deming_TaskD_TaskND, int_var.Mean_TaskD_TaskND, int_var.Mean_ROI_TaskD_TaskND, int_var.Zscore_TaskD, int_var.SVM_TaskD, int_var.LDC_TaskD]';
p25_mat=[int_res.Attention_Modulation, int_p25.Deming_TaskD_TaskND, int_p25.Mean_TaskD_TaskND, int_p25.Mean_ROI_TaskD_TaskND, int_p25.Zscore_TaskD, int_p25.SVM_TaskD, int_p25.LDC_TaskD]';
p75_mat=[int_res.Attention_Modulation, int_p75.Deming_TaskD_TaskND, int_p75.Mean_TaskD_TaskND, int_p75.Mean_ROI_TaskD_TaskND,int_p75.Zscore_TaskD, int_p75.SVM_TaskD, int_p75.LDC_TaskD]';

mean_res = mean(res_mat,2);
var_mat_new = mean((p75_mat - p25_mat)/2,2)./mean_res;
for i=1:3
    res_mat(:,i) = res_mat(:,i)./mean_res;
end
contri_vects =[1 0 -1; 0.5 -1 0.5];
contri_mat = contri_vects*res_mat';
contri_mat = contri_mat';



%Plot the data
figure
att_plot=bar(contri_mat);
hold on
for b=1:2
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(contri_mat, 1);
nbars = size(contri_mat, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, contri_mat(:,i), var_mat_new, 'k.');
end

%Tidying up the plot and adding labels
set(gca, 'XTickLabel', {'Ground Truth','Deming Regression','Ratio of individual voxels', 'Ratio of entire ROI', 'Z-scoring', 'SVM classification', 'LDC'});
set(gca,'XTickLabelRotation',20);
%ylim([0 4]);
ylabel('Laminar Contributions')
x0=10;
y0=10;
width=950;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend('Superficial Bias', 'Attentional Modulation','location','northwest')

%save figure
fname = sprintf('att_%g_%g_%g_contri_plot.png',int_attention);
saveas(gcf,fname,'png');
hold off
