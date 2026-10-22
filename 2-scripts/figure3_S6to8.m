%% script to generate figure 3 and supp figures 6 to 8
% author: taylor nguyen

clear all;
close all;
clc;

cd('/Users/taynguyen/Dropbox/_phage/github_upload/2-scripts/')

%% figure 3a
% plot niche overlap of each species with Ef

% load untargeted metabolomics data
% run "/2-scripts/metab_generate_niches.m" to generate dataframe
load('../1-data/metabolomics_isolates_consumed_all.mat', 'metab_consumed')

names = {"TT1_", "TT2_", "TT3_", "TT4_", "TT5_", "TT6_", ...
         "TT7_", "TT8_", "TT9_", "TT10_", "TT11_", "TT12", ...
         "TT13_", "TT14_", "TT15_", "Klebsiella", "Salmonella", "MG1655"};

ef_shared_niche = [];
ef_shared_niche_pct = [];

% put in phylogenetic order
% klebsiella pneumoniae
% e hirae, e faecalis, e faecium
% e fergusonii (omit)
% b theta, fragilis, uniformis, distasonis
% f plautii, c scindens, c hylemonae, c symbiosum, c clostridioforme
% b. producta, c hathewayi

isols = [18, 16, 1, 10, 8, 4, 12, 13, 15, 13, 7, 9, 3, 5, 6, 11];
names(isols)

colors2plot = [84,39,143; 84,39,143; ...
               140,81,10; 140,81,10; 140,81,10; ...
               253,141,60; 253,141,60; 253,141,60; 107,174,214;
               33,113,181; 65,171,93; 65,171,93; 65,171,93; 65,171,93;
               65,171,93; 65,171,93] ./ 256;

for i = 1:length(isols)

    isol = isols(i);
    metab_consumed_filt = [];
    metab_consumed_filt = metab_consumed(:,[2, isol]);
    metab_consumed_filt(:,3) = metab_consumed_filt(:,1)==1 & metab_consumed_filt(:,2)==1;
    t = sum(metab_consumed_filt);
    
    ef_shared_niche(1,i) = t(1); % ef all
    ef_shared_niche(2,i) = t(2); % other isolate all
    ef_shared_niche(3,i) = t(3); % shared
    ef_shared_niche(4,i) = t(3) ./ t(1); % pct of ef niche that is shared

end

ef_shared_niche

figure
b = bar(ef_shared_niche(4,:), 'EdgeColor', 'none');
b.FaceColor = "flat";

for i = 1:length(colors2plot)
    b.CData(i,:) = colors2plot(i,:);
end

ylim([0 1.1])
yline(1, '--')

xticklabels(["Kp", "Eh", "Efs", "Efm", "Bt", "Bf", "Bu", "Pd", "Fp", "Csc", ...
                "Chy","Csy","Cc","Bp","Cha"]);

ylabel("Pct. metabolites shared with Ef (100X cutoff)")

% save figure
w=3;
h=2;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', ...
    'PaperSize', [w, h], 'DefaultFigureRenderer', 'painters')

fname = "../3-figs/fig3/15mb_shared_metabolites_wEf";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%%