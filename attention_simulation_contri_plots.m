function attention_simulation_contri_plots(filename)
% Plotting results from the output of attention_simulation.m
% This function queries the output for the parameters we are interested in
% and generate plots comparing the ground truth against the various
% estimates using the different metrics.

% define which values we are interested in
% Change the values here to query different results
physio_sigma_list = repmat(8,1,3);
thermal_sigma_list = repmat(12,1,3);
superficial_bias = [2, 1.5, 1];
attentional_modulation = [3 2 3];             

cmap = [251,180,174;
    179,205,227];
cmap = cmap/255;

res = load(filename); 
res = cat(1,res.all_res{:});


sz = [3,12];
var_types = repmat({'double'},1,12);
var_names = {'Physio_Noise','Thermal_Noise','Superficial_Bias','attentional_modulation','Zscore_TaskD','SVM_TaskD','LDC_TaskD','Mean_TaskD_TaskND','Mean_ROI_TaskD_TaskND','Deming_TaskD_TaskND','Real_BOLD_TaskD','Measured_BOLD_TaskD'};
res_table = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names);
var_table = res_table;
med_table = res_table;
p25_table = res_table;
p75_table = res_table;


%Extract the 3 rows
for i = 1:3
    ind = find((res(:,1)==physio_sigma_list(i) & res(:,2)==thermal_sigma_list(i) & res(:,3)==superficial_bias(i) & res(:,4)==attentional_modulation(i)));
    res_table(i,:)=[{physio_sigma_list(i),thermal_sigma_list(i),superficial_bias(i),attentional_modulation(i)},num2cell(nanmean(res(ind,5:end)))];
    var_table(i,:)=[{physio_sigma_list(i),thermal_sigma_list(i),superficial_bias(i),attentional_modulation(i)},num2cell(nanstd(res(ind,5:end),0,1))];
    med_table(i,:)=[{physio_sigma_list(i),thermal_sigma_list(i),superficial_bias(i),attentional_modulation(i)},num2cell(nanmedian(res(ind,5:end),1))];
    p25_table(i,:)=[{physio_sigma_list(i),thermal_sigma_list(i),superficial_bias(i),attentional_modulation(i)},num2cell(prctile(res(ind,5:end),25))];
    p75_table(i,:)=[{physio_sigma_list(i),thermal_sigma_list(i),superficial_bias(i),attentional_modulation(i)},num2cell(prctile(res(ind,5:end),75))];
end


%Combine the data into a results and error matrix
res_mat=[res_table.attentional_modulation, res_table.Deming_TaskD_TaskND, res_table.Mean_TaskD_TaskND, res_table.Mean_ROI_TaskD_TaskND, res_table.Zscore_TaskD, res_table.SVM_TaskD, res_table.LDC_TaskD]';
var_mat=[zeros(3,1), var_table.Deming_TaskD_TaskND, var_table.Mean_TaskD_TaskND, var_table.Mean_ROI_TaskD_TaskND, var_table.Zscore_TaskD, var_table.SVM_TaskD, var_table.LDC_TaskD]';
p25_mat=[res_table.attentional_modulation, p25_table.Deming_TaskD_TaskND, p25_table.Mean_TaskD_TaskND, p25_table.Mean_ROI_TaskD_TaskND, p25_table.Zscore_TaskD, p25_table.SVM_TaskD, p25_table.LDC_TaskD]';
p75_mat=[res_table.attentional_modulation, p75_table.Deming_TaskD_TaskND, p75_table.Mean_TaskD_TaskND, p75_table.Mean_ROI_TaskD_TaskND,p75_table.Zscore_TaskD, p75_table.SVM_TaskD, p75_table.LDC_TaskD]';

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
ylim([-0.2 1]);
ylabel('Laminar Contributions')
x0=10;
y0=10;
width=950;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend('Superficial Bias', 'Attentional Modulation','location','northwest')

%save figure
fname = sprintf('att_%g_%g_%g_contri_plot.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off
