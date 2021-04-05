function tstat_comparison_plots(real_results, simulated_results)
% visualise matching between real and simulated fMRI data in terms of regional-mean T
% statistics and regional standard deviations.
%
% The inputs are paths to relevant .mat files. If they are undefined, the examples from
% the sample_results folder are used.
%
% tstat_comparison_plots(real_results, simulated_results)

if ~exist('real_results', 'var') || isempty(real_results)
    % default to example pre-computed results
    real_results = load(fullfile(fileparts(mfilename('fullpath')), ...
        'sample_results', 'real_fMRI_results.mat'));
end

if ~exist('simulated_results', 'var') || isempty(simulated_results)
    simulated_results = load(fullfile(fileparts(mfilename('fullpath')), ...
        'sample_results', 'att_sim_results.mat'));
end

tstat_mean = NaN([7, 2]);
tstat_std = NaN([7, 2]);
tstat_mean_std = NaN([7, 2]);
tstat_std_std = NaN([7, 2]);

for layer = 1:3 % 1=deep, 2=middle, 3=superficial
    % Real Data
    layer_names = {'Deep','Middle','Superficial'};
    fprintf('%s Layer \n',layer_names{layer})
    tmean_dplus_real = zeros(6,1);
    tmean_dminus_real = zeros(6,1);
    tstd_dplus_real = zeros(6,1);
    tstd_dminus_real = zeros(6,1);
    for i=1:6
        tmean_dplus_real(i) = real_results.results(i).estimates(layer).tstat_dplus_mean;
        tstd_dplus_real(i) = real_results.results(i).estimates(layer).tstat_dplus_std;

        tmean_dminus_real(i) = real_results.results(i).estimates(layer).tstat_dminus_mean;
        tstd_dminus_real(i) = real_results.results(i).estimates(layer).tstat_dminus_std;
    end

    fprintf('Real fMRI data \n')
    fprintf('Dplus Tstat Mean: %f%s%f \n', mean(tmean_dplus_real), char(177), std(tmean_dplus_real))
    fprintf('Dplus Tstat Std: %f%s%f \n', mean(tstd_dplus_real), char(177), std(tstd_dplus_real))
    fprintf('Dminus Tstat Mean: %f%s%f \n', mean(tmean_dminus_real), char(177), std(tmean_dminus_real))
    fprintf('Dminus Tstat Std: %f%s%f \n', mean(tstd_dminus_real), char(177), std(tstd_dminus_real))
    
    tstat_mean(4-layer,1) =  mean(tmean_dplus_real);
    tstat_std(4-layer,1) =  mean(tstd_dplus_real);
    tstat_mean(8-layer,1) =  mean(tmean_dminus_real);
    tstat_std(8-layer,1) =  mean(tstd_dminus_real);
    
    tstat_mean_std(4-layer,1) =  std(tmean_dplus_real);
    tstat_std_std(4-layer,1) =  std(tstd_dplus_real);
    tstat_mean_std(8-layer,1) =  std(tmean_dminus_real);
    tstat_std_std(8-layer,1) =  std(tstd_dminus_real);
    
    % Simulation Data
    sim_size = numel(simulated_results.results);
    tmean_dplus_sim = zeros(sim_size,1);
    tmean_dminus_sim = zeros(sim_size,1);
    tstd_dplus_sim = zeros(sim_size,1);
    tstd_dminus_sim = zeros(sim_size,1);
    for i=1:sim_size
        tmean_dplus_sim(i) = simulated_results.results(i).estimates(layer*2).tstat_dplus_mean;
        tstd_dplus_sim(i) = simulated_results.results(i).estimates(layer*2).tstat_dplus_std;

        tmean_dminus_sim(i) = simulated_results.results(i).estimates(layer*2).tstat_dminus_mean;
        tstd_dminus_sim(i) = simulated_results.results(i).estimates(layer*2).tstat_dminus_std;
    end

    fprintf('Sim fMRI data \n')
    fprintf('Dplus Tstat Mean: %f%s%f \n', mean(tmean_dplus_sim), char(177), std(tmean_dplus_sim))
    fprintf('Dplus Tstat Std: %f%s%f \n', mean(tstd_dplus_sim), char(177), std(tstd_dplus_sim))
    fprintf('Dminus Tstat Mean: %f%s%f \n', mean(tmean_dminus_sim), char(177), std(tmean_dminus_sim))
    fprintf('Dminus Tstat Std: %f%s%f \n', mean(tstd_dminus_sim), char(177), std(tstd_dminus_sim))
    fprintf('\n')
    
    tstat_mean(4-layer,2) =  mean(tmean_dplus_sim);
    tstat_std(4-layer,2) =  mean(tstd_dplus_sim);
    tstat_mean(8-layer,2) =  mean(tmean_dminus_sim);
    tstat_std(8-layer,2) =  mean(tstd_dminus_sim);
    
    tstat_mean_std(4-layer,2) =  std(tmean_dplus_sim);
    tstat_std_std(4-layer,2) =  std(tstd_dplus_sim);
    tstat_mean_std(8-layer,2) =  std(tmean_dminus_sim);
    tstat_std_std(8-layer,2) =  std(tstd_dminus_sim);
end

cmap = [251,180,174;
    179,205,227];
cmap = cmap/255;


figure
att_plot=bar(tstat_mean);
hold on
for b=1:2
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(tstat_mean, 1);
nbars = size(tstat_mean, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, tstat_mean(:,i), tstat_mean_std(:,i), 'k.');
end

%Tidying up the plot and adding labels
xticks([1 2 3 5 6 7])
set(gca, 'XTickLabel', {'Superficial','Middle', 'Deep'});
set(gca,'XTickLabelRotation',20);
ylabel({'Region-mean T statistic', '(mean \pm1 standard deviation)'});
x0=10;
y0=10;
width=700;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend('Real','Simulated','location','northwest')

%save figure
fname = 'tstat.png';
saveas(gcf,fname,'png');
hold off

figure
att_plot=bar(tstat_std);
hold on
for b=1:2
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(tstat_std, 1);
nbars = size(tstat_std, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, tstat_std(:,i), tstat_std_std(:,i), 'k.');
end

%Tidying up the plot and adding labels
xticks([1 2 3 5 6 7])
set(gca, 'XTickLabel', {'Superficial','Middle', 'Deep'});
set(gca,'XTickLabelRotation',20);
ylabel('T-statistic (Standard Deviation)')
ylabel({'Standard deviation of region-mean T statistic', '(mean \pm1 standard deviation)'});
x0=10;
y0=10;
width=700;
height=500;
set(gcf,'position',[x0,y0,width,height])
legend('Real','Simulated','location','northwest')

%save figure
fname = 'tstat_std.png';
saveas(gcf,fname,'png');
hold off
