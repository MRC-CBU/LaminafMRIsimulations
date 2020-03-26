function attention_simulation_plots(filename)
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

cmap = [8 81 156;
    107 174 214;
    189 215 231];
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
   
% Split into multiple plots:
% 1. Att Modulation/Deming/Mean/MeanROI
% 2. SVM
% 3. LDC
% 4. Z-score

%Combine the data into a results and error matrix
res_mat=[res_table.attentional_modulation, res_table.Deming_TaskD_TaskND, res_table.Mean_TaskD_TaskND, res_table.Mean_ROI_TaskD_TaskND]';
var_mat=[zeros(3,1), var_table.Deming_TaskD_TaskND, var_table.Mean_TaskD_TaskND, var_table.Mean_ROI_TaskD_TaskND]';
p25_mat=[res_table.attentional_modulation, p25_table.Deming_TaskD_TaskND, p25_table.Mean_TaskD_TaskND, p25_table.Mean_ROI_TaskD_TaskND]';
p75_mat=[res_table.attentional_modulation, p75_table.Deming_TaskD_TaskND, p75_table.Mean_TaskD_TaskND, p75_table.Mean_ROI_TaskD_TaskND]';

%Plot the data
figure
att_plot=bar(res_mat);
hold on
for b=1:3
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(res_mat, 1);
nbars = size(res_mat, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, res_mat(:,i), res_mat(:,i)-p25_mat(:,i),p75_mat(:,i)-res_mat(:,i), 'k.');
end

%Tidying up the plot and adding labels
set(gca, 'XTickLabel', {'Ground Truth','Deming Regression','Ratio of individual voxels', 'Ratio of entire ROI'});
set(gca,'XTickLabelRotation',20);
ylim([0 4]);
ylabel('Attentional Modulation (A.U.)')
x0=10;
y0=10;
width=680;
height=520;
set(gcf,'position',[x0,y0,width,height])
legend('Superficial','Mid','Deep','location','northwest')

%save figure
fname = sprintf('att_%g_%g_%g_plot_1.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off
% 

%Plot the data
figure
att_plot=bar(res_table.Zscore_TaskD);
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
    errorbar(x, res_table.Zscore_TaskD(:,i), res_table.Zscore_TaskD(:,i)-p25_table.Zscore_TaskD(:,i),p75_table.Zscore_TaskD(:,i)-res_table.Zscore_TaskD(:,i), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'Z-scoring'});
set(gca,'XTickLabelRotation',20);
%ylim([0 10]);
ylabel('Z-scoring (A.U.)')
x0=10;
y0=10;
width=190;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_2.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off


%Plot the data
figure
att_plot=bar(res_table.SVM_TaskD);
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
    errorbar(x, res_table.SVM_TaskD(:,i), res_table.SVM_TaskD(:,i)-p25_table.SVM_TaskD(:,i),p75_table.SVM_TaskD(:,i)-res_table.SVM_TaskD(:,i), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'SVM Classification'});
set(gca,'XTickLabelRotation',20);
ylim([60 100]);
ylabel('Classification Accuracy (%)')
x0=10;
y0=10;
width=185;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_3.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off


%Plot the data
figure
att_plot=bar(res_table.LDC_TaskD);
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
    errorbar(x, res_table.LDC_TaskD(:,i), res_table.LDC_TaskD(:,i)-p25_table.LDC_TaskD(:,i),p75_table.LDC_TaskD(:,i)-res_table.LDC_TaskD(:,i), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'LDC'});
set(gca,'XTickLabelRotation',20);
%ylim([60 100]);
ylabel('Linear Discriminant Contrast (A.U.)')
x0=10;
y0=10;
width=180;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = sprintf('att_%g_%g_%g_plot_4.png',attentional_modulation);
saveas(gcf,fname,'png');
hold off
close all
end