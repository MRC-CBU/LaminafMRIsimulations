
for layer = 1:3 % 1=deep, 2=middle, 3=superficial
    % Real Data
    load('real_fMRI_results.mat')
    layer_names = {'Deep','Middle','Superficial'};
    fprintf('%s Layer \n',layer_names{layer})
    tmean_dplus_real = zeros(6,1);
    tmean_dminus_real = zeros(6,1);
    tstd_dplus_real = zeros(6,1);
    tstd_dminus_real = zeros(6,1);
    for i=1:6
        tmean_dplus_real(i) = results_exclude(i).estimates(layer).tstat_dplus_mean;
        tstd_dplus_real(i) = results_exclude(i).estimates(layer).tstat_dplus_std;

        tmean_dminus_real(i) = results_exclude(i).estimates(layer).tstat_dminus_mean;
        tstd_dminus_real(i) = results_exclude(i).estimates(layer).tstat_dminus_std;

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
    load('att_sim_results.mat')
    sim_size = 20;
    tmean_dplus_sim = zeros(sim_size,1);
    tmean_dminus_sim = zeros(sim_size,1);
    tstd_dplus_sim = zeros(sim_size,1);
    tstd_dminus_sim = zeros(sim_size,1);
    for i=1:sim_size
        tmean_dplus_sim(i) = results(i).estimates(layer*2).tstat_dplus_mean;
        tstd_dplus_sim(i) = results(i).estimates(layer*2).tstat_dplus_std;

        tmean_dminus_sim(i) = results(i).estimates(layer*2).tstat_dminus_mean;
        tstd_dminus_sim(i) = results(i).estimates(layer*2).tstat_dminus_std;
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
xticks([1 2 3 4 5 6 7])
set(gca, 'XTickLabel', {'Superficial','Middle', 'Deep',' '});
set(gca,'XTickLabelRotation',20);
ylabel('T-statistic (Mean)')
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
xticks([1 2 3 4 5 6 7])
set(gca, 'XTickLabel', {'Superficial','Middle', 'Deep',' '});
set(gca,'XTickLabelRotation',20);
ylabel('T-statistic (Standard Deviation)')
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
