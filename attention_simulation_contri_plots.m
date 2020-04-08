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
res = res.results;

%Find the indices for the 3 rows
match = NaN(3,size(res(1).params,2));
for i=1:3
    match(i,:) = [res(1).params.physio_sigma]==physio_sigma_list(i) & [res(1).params.thermal_sigma]==thermal_sigma_list(i) & [res(1).params.superficial_bias]==superficial_bias(i) & [res(1).params.attentional_modulation]==attentional_modulation(i);
end

var_names=fieldnames(res(1).estimates);
res_table = struct();
for iter=1:size(res,2)
    estimates = res(iter).estimates;
    for i=1:3
        cur_est = estimates(match(i,:)==1);
        for j=1:size(var_names,1)
            if size([estimates.(var_names{j})],1)==1 %skipping real bold response and measured bold response because they are nvox*iter matrices
                res_table.(var_names{j})(iter,i) = cur_est.(var_names{j});
            end
        end
    end
end

%Adding in Ground Truth
res_table.ground_truth(1,:) = attentional_modulation;

%Order res_table in plotting order
order = {'ground_truth','deming_est','raw_ratio_est','ROI_ratio_est','zscore','SVM','LDC'};
res_table = orderfields(res_table,order);
var_names=fieldnames(res_table);

%demean the data
for i=1:size(var_names,1)
    for j=1:size(res_table.(var_names{i}),1)
        res_table.(var_names{i})(j,:) = res_table.(var_names{i})(j,:)/mean(res_table.(var_names{i})(j,:));
    end
end

contri_vects =[1 0 -1; 0.5 -1 0.5]';
mean_contri=NaN(size(var_names,2),2);
p25_contri=mean_contri;
p75_contri=mean_contri;

for i=1:size(var_names,1)
    contri=res_table.(var_names{i}) * contri_vects;
    mean_contri(i,:)=mean(contri,1);
    p25_contri(i,:)=prctile(contri,25,1);
    p75_contri(i,:)=prctile(contri,75,1);
end

%Plot the data
figure
att_plot=bar(mean_contri);
hold on
for b=1:2
    att_plot(b).FaceColor = cmap(b,:);
end

%Add error bars
ngroups = size(mean_contri, 1);
nbars = size(mean_contri, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, mean_contri(:,i), mean_contri(:,i)-p25_contri(:,i),p75_contri(:,i)-mean_contri(:,i), 'k.');
end

%Tidying up the plot and adding labels
set(gca, 'XTickLabel', {'Ground Truth','Deming Regression','Ratio of individual voxels', 'Ratio of entire ROI', 'Z-scoring', 'SVM classification', 'LDC'});
set(gca,'XTickLabelRotation',20);
ylim([-0.2 1.2]);
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
