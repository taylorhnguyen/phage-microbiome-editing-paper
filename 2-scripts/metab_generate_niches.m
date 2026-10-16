%% generate niche overlaps using untargeted metabolomics data

clear all;
close all;
clc;

%% import untargeted metabolomics data

t_neg = readtable('../1-data/metabolomics_datasheets/TaylorNguyen_KCHU02_CombineSubmit_NoAnnotations.xlsx', ...
    'ReadVariableNames', true, 'Sheet', 'Neg');
t_pos = readtable('../1-data/metabolomics_datasheets/TaylorNguyen_KCHU02_CombineSubmit_NoAnnotations.xlsx', ...
    'ReadVariableNames', true, 'Sheet', 'Pos');

% combine pos and negative into one dataframe
t_neg = t_neg(8:end, 7:end);
t_pos = t_pos(22:end, 7:end); 

t = [t_neg; t_pos];
snames = t.Properties.VariableNames;
metabnames = t.Properties.VariableNames;

%% organize data
metab_all = [];
for i = 1:72
    metab_all(:,i) = table2array(t(:,i));
end

%% take log10 then the mean

% take log10 of each value
metab_all_log10 = metab_all;
metab_all_log10(metab_all_log10==0)=1;
metab_all_log10 = log10(metab_all_log10);

% get mean log10 BHI blank values
idx_BHI = find(contains(snames,"BHI_blank"));
metab_BHI = metab_all(:,idx_BHI);
metab_BHI_mean = mean(metab_BHI, 2);
metab_BHI_log10 = metab_all_log10(:,idx_BHI);
metab_BHI_log10_mean = mean(metab_BHI_log10, 2);

%% take the mean for each species
metab_log10_mean = [];
metab_mean = [];

names = {"TT1_", "TT2_", "TT3_", "TT4_", "TT5_", "TT6_", ...
         "TT7_", "TT8_", "TT9_", "TT10_", "TT11_", "TT12", ...
         "TT13_", "TT14_", "TT15_", "Klebsiella", "Salmonella", "MG1655"};

n = length(names);
m = length(metab_all);

for i = 1:n
    idx_isol = find(contains(snames, names{i}))
    metab_log10_mean(:,i) = mean(metab_all_log10(:, idx_isol), 2);
    metab_mean(:,i) = mean(metab_all(:, idx_isol), 2);
end

% get fold-change from BHI / normalize to BHI 
metab_log10_mean_norm = metab_log10_mean - metab_BHI_log10_mean;

%% get p-values (FDR) between isolate vs. BHI
pvals = zeros(m,n);
qvals = zeros(m,n);

for i = 1:n
    for j = 1:m
        idx_isol = find(contains(snames, names{i}));
        x = metab_all_log10(j, idx_isol);
        y = metab_BHI_log10(j,:);
        [~,p] = ttest2(x, y);
        pvals(j,i) = p;
    end
end

% apply benjamini-hochberg (BH) correction
for i = 1:n
    qvals(:,n) = mafdr(pvals(:,n), 'BHFDR', true);
end

%% find unique features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set threshold values
% depleted: decreased >100-fold compared to BHI & difference is significant 
% (p-value < 0.05) from two-sample t-test
% produced: increaed >10-fold compared to BHI & difference is significant 
% (p-value < 0.05) from two-sample t-test
% and feature >5,000 in BHI (initial value)

fc_thold = -2; % 100-fold cut-off
pval_thold = 0.05;
qval_thold = 0.05;
count_thold = 5000;

% make consumption matrix
metab_consumed = zeros(m,2);

for i = 1:18
    idx_isol = find(contains(snames, names{i}));
    isol_consumed(:,1) = metab_log10_mean_norm(:,i) < fc_thold; % check consumption
    isol_consumed(:,2) = qvals(:,i) < qval_thold; % check qval (after FDR correction)
    isol_consumed(:,3) = metab_BHI_mean > count_thold; % check starting counts in BHI

    % find features that pass all tests
    isol_consumed(:,4) = isol_consumed(:,1)==1 & ...
                         isol_consumed(:,2)==1 & ...
                         isol_consumed(:,3)==1;
    metab_consumed(:,i) = isol_consumed(:,4);
end

metab_consumed_total = sum(metab_consumed)

% save for downstream analysis
save("../1-data/metabolomics_isolates_consumed_all.mat", "metab_consumed")

%%