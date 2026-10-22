%% script to generate figure 5 and supp figures 9
% author: taylor nguyen

clear all;
close all;
clc;

cd('/Users/taynguyen/Dropbox/_phage/github_upload/2-scripts/')

%%

%% figure 5d 
% compare ef and stm metabolomics

% load untargeted metabolomics data
% run "/2-scripts/metab_generate_niches.m" to generate dataframe
load('../1-data/metabolomics_isolates_consumed_all.mat', 'metab_consumed')

names = {"TT1_", "TT2_", "TT3_", "TT4_", "TT5_", "TT6_", ...
         "TT7_", "TT8_", "TT9_", "TT10_", "TT11_", "TT12", ...
         "TT13_", "TT14_", "TT15_", "Klebsiella", "Salmonella", "MG1655"};

% make a venn diagram between ef and stm
metab_consumed_filt = [];
metab_consumed_filt = metab_consumed(:,[2, 17]);
metab_consumed_filt(:,3) = metab_consumed_filt(:,1)==1 & metab_consumed_filt(:,2)==1;
t = sum(metab_consumed_filt)

% pct. unique to stm
(t(2)-t(3))/t(2)

ef_unique = t(1)-t(3)
stm_unique = t(2)-t(3)
ef_stm_shared = t(3)

figure
subplot(1,2,1)
vals = [ef_unique, stm_unique, ef_stm_shared];
venn(vals)
xticks([])
yticks([])
xticklabels([])
yticklabels([])

% make venn diagram between kp and stm
metab_consumed_filt = [];
metab_consumed_filt = metab_consumed(:,[16, 17]);
metab_consumed_filt(:,3) = metab_consumed_filt(:,1)==1 & metab_consumed_filt(:,2)==1;
t = sum(metab_consumed_filt)

% pct. unique to stm
(t(2)-t(3))/t(2)

kp_unique = t(1)-t(3)
stm_unique = t(2)-t(3)
kp_stm_shared = t(3)

subplot(1,2,2)
vals = [kp_unique, stm_unique, kp_stm_shared];
venn(vals)
xticks([])
yticks([])
xticklabels([])
yticklabels([])

% save figure
w=6;
h=1.5;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = "../3-figs/fig5/venn_efvsstm_kpvsstm";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%%