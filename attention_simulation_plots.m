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
res = res.results;

%Extract the 3 rows
params = [res.params];
estimates = [res.estimates];
var_names=fieldnames(estimates);
for i = 1:3
    match = [params.physio_sigma]==physio_sigma_list(i) & [params.thermal_sigma]==thermal_sigma_list(i) & [params.superficial_bias]==superficial_bias(i) & [params.attentional_modulation]==attentional_modulation(i);
    
    cur_estimates = estimates(match==1);
    
    for j=1:size(var_names,1)
        if size([cur_estimates.(var_names{j})],1)==1 %skipping real bold response and measured bold response because they are nvox*iter matrices
            res_table.(var_names{j})(i)=mean([cur_estimates.(var_names{j})]);
            var_table.(var_names{j})(i)=std([cur_estimates.(var_names{j})]);
            p25_table.(var_names{j})(i)=prctile([cur_estimates.(var_names{j})],25);
            p75_table.(var_names{j})(i)=prctile([cur_estimates.(var_names{j})],75);
        end
    end
end
   
% Split into multiple plots:
% 1. Att Modulation/Deming/Mean/MeanROI
% 2. SVM
% 3. LDC
% 4. Z-score

%Combine the data into a results and error matrix
res_mat=[attentional_modulation; res_table.deming_est; res_table.raw_ratio_est; res_table.ROI_ratio_est];
var_mat=[zeros(1,3); var_table.deming_est; var_table.raw_ratio_est; var_table.ROI_ratio_est];
p25_mat=[attentional_modulation; p25_table.deming_est; p25_table.raw_ratio_est; p25_table.ROI_ratio_est];
p75_mat=[attentional_modulation; p75_table.deming_est; p75_table.raw_ratio_est; p75_table.ROI_ratio_est];

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
att_plot=bar(res_table.zscore);
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
    errorbar(x, res_table.zscore(i,:), res_table.zscore(i,:)-p25_table.zscore(i,:),p75_table.zscore(i,:)-res_table.zscore(i,:), 'k.');
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
att_plot=bar(res_table.SVM);
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
    errorbar(x, res_table.SVM(i,:), res_table.SVM(i,:)-p25_table.SVM(i,:),p75_table.SVM(i,:)-res_table.SVM(i,:), 'k.');
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'SVM Classification'});
set(gca,'XTickLabelRotation',20);
ylim([80 100]);
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
att_plot=bar(res_table.LDC);
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
    errorbar(x, res_table.LDC(i,:), res_table.LDC(i,:)-p25_table.LDC(i,:),p75_table.LDC(i,:)-res_table.LDC(i,:), 'k.');
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