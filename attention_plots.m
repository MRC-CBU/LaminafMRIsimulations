function attention_plots(res_file)
% Plot the real fMRI results in a similar manner to the simulated data (see
% attention_simulation_plots.m).
%
% If the input real_results is undefined, we load the sample result from
% sample_results/real_fMRI_results.mat
%
% attention_plots(filename)

if ~exist('res_file', 'var') || isempty(res_file)
    % default to example pre-computed results
    res_file = fullfile(fileparts(mfilename('fullpath')),'sample_results', 'real_fMRI_results.mat');
end

cmap = [8 81 156;
    49 130 189;
    107 174 214;
    189 215 231];



lmap = [31 120 180; %blue
     51 160 44; %green
     227 26 28; %pink
     200 180 0; %gold
    106 61 154; %purple
    0 0 0]; %black

cmap = cmap/255;
lmap = lmap/255;

res = load(res_file);  
res = res.results;

%convert results into matrix form
res_table = struct();
var_names=fieldnames(res(1).estimates);
for j=1:size(var_names,1)
    temp_data = zeros(6,3);
    for i=1:6
        for k=1:3
            if strcmp(var_names{j},'full')==0 
                temp_data(i,4-k)=res(i).estimates(k).(var_names{j});
            end
        end
    end
    res_table.(var_names{j})=temp_data;
end

% Split into multiple plots:
% 1. Deming/Mean/MeanROI
% 2. SVM
% 3. LDC
% 4. Z-score
% 5. L2 norm

%--------------------------------------------------------------------------------------------
%Plotting selectivity
data(:,1,:)=res_table.deming_ratio;
data(:,2,:)=res_table.raw_ratio; 
data(:,3,:)=res_table.ROI_ratio;

figure
hold on
task_sum = bar(squeeze(nanmedian(data,1)));
for b=1:3
    task_sum(b).FaceColor = cmap(b,:);
end

%Tidying up the plot and adding labels
xticks([1 2 3])
set(gca, 'XTickLabel', {'Deming Regression','Voxel Ratio', 'ROI Ratio'});
set(gca,'XTickLabelRotation',20);
% ylim([-0.5 2.5]);
ylabel('Median selectivity (A.U.)')
x0=10;
y0=10;
width=570;
height=520;
set(gcf,'position',[x0,y0,width,height])
box on

pause(0.1); %pause allows the figure to be created
for layer = 1:3
    for t = 0:2:0
        clear xData
        for b=1:3
            xData(b) = task_sum(t+b).XData(1,layer)+task_sum(t+b).XOffset;
        end
        for sub=1:6
            sp = plot(xData,[squeeze(data(sub,layer,:))],'-o','MarkerSize',5);
            sp.Color = lmap(sub,:);
        end
    end
end
legend('Superficial','Middle','Deep','location','northwest')

%save figure
fname = 'plot_1.png';
saveas(gcf,fname,'png');
hold off


%--------------------------------------------------------------------------------------------
%Plotting Z scoring
data = res_table.zscore;
figure
att_plot=bar(squeeze(nanmedian(data,1)));
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
clear xData
for b=1:3
    xData(b) = att_plot.XData(b);
end
for sub=1:6
    sp = plot(xData,[squeeze(data(sub,:))],'-o','MarkerSize',5);
    sp.Color = lmap(sub,:);
end

%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'Z-scoring'});
set(gca,'XTickLabelRotation',20);
%ylim([0 10]);
ylabel('Median region-mean contrast estimate (A.U.)')
x0=10;
y0=10;
width=190;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = 'plot_2.png';
saveas(gcf,fname,'png');
hold off


%--------------------------------------------------------------------------------------------
%Plotting SVM
data = res_table.SVM;
figure
att_plot=bar(squeeze(nanmedian(data,1)));
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
clear xData
for b=1:3
    xData(b) = att_plot.XData(b);
end
for sub=1:6
    sp = plot(xData,[squeeze(data(sub,:))],'-o','MarkerSize',5);
    sp.Color = lmap(sub,:);
end


%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'SVM Classification'});
set(gca,'XTickLabelRotation',20);
ylim([50 100]);
ylabel('Median classification accuracy (%)')
x0=10;
y0=10;
width=185;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = 'plot_3.png';
saveas(gcf,fname,'png');
hold off


%--------------------------------------------------------------------------------------------
%Plotting LDC
data = res_table.LDC;
figure
att_plot=bar(squeeze(nanmedian(data,1)));
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
clear xData
for b=1:3
    xData(b) = att_plot.XData(b);
end
for sub=1:6
    sp = plot(xData,[squeeze(data(sub,:))],'-o','MarkerSize',5);
    sp.Color = lmap(sub,:);
end


%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'LDC'});
set(gca,'XTickLabelRotation',20);
%ylim([60 100]);
ylabel('Median linear discriminant contrast (A.U.)')
x0=10;
y0=10;
width=180;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = 'plot_4.png';
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------
%Plotting l2_dplus
data = res_table.l2_dplus;
figure
att_plot=bar(squeeze(nanmedian(data,1)));
hold on
att_plot.FaceColor = 'flat';
for b=1:3
    att_plot.CData(b,:) = cmap(b,:);
end


%Add error bars
clear xData
for b=1:3
    xData(b) = att_plot.XData(b);
end
for sub=1:6
    sp = plot(xData,[squeeze(data(sub,:))],'-o','MarkerSize',5);
    sp.Color = lmap(sub,:);
end


%Tidying up the plot and adding labels
xticks([2])
set(gca, 'XTickLabel', {'L2 norm'});
set(gca,'XTickLabelRotation',20);
%ylim([60 100]);
ylabel('Region-mean contrast estimate (A.U.)');
x0=10;
y0=10;
width=180;
height=500;
set(gcf,'position',[x0,y0,width,height])

%save figure
fname = 'plot_5.png';
saveas(gcf,fname,'png');
hold off

%--------------------------------------------------------------------------------------------


close all
end
