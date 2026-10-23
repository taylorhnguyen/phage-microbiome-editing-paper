%% script to generate figure 2 and supp figures 4 to 5
% author: taylor nguyen

clear all;
close all;
clc;

cd('/Users/taynguyen/Dropbox/_phage/github_upload/2-scripts/')

%%  load master variables

tay4_color = [201, 148, 199] / 255;
kp_color = [223,101,176] ./ 255;
cocktail_color = [233,14,139] ./ 255;

load ColorsKat % loads colors_kat, a cell array with color codes for each L and OTUs_kat, the corresponfding OTUs for each color in colors_kat

rows = ["A","B","C","D","E","F","G","H"];
cols = ["1","2","3","4","5","6","7","8","9","10","11","12"];

%% supp: figure s4a
% solid isolation of Tay4-treated Ef from monoculture

% load 96-well plate OD600 growth curves
load('../1-data/fig2_250325_res_Ef_solid_isolation.mat', 'res_Ef_solid_isolation')

figure
subplot(1,2,1)
AUC_all = [];
count = 1;

for c = 3:4 
    for r = 1:8

        % plot growth curves / plot control first
        well = strcat(rows(r), cols(c));
        plot(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank, 'color', [0, 0, 0]); hold all;
        AUC_all(count,1) = trapz(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank);

        % plot +Tay4
        well = strcat(rows(r), cols(c+6));
        plot(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank, 'color', tay4_color);
        AUC_all(count,2) = trapz(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank);

        count = count + 1;
    end
end

% figure details
xlim([0 48])
set(gca, 'yscale', 'log')
ylim([5*10^-4 2])
yticks([10^-2 10^0])
yticklabels({'10^{-2}', '10^0'})

% plot AUC
subplot(1,2,2)
AUC_norm = AUC_all ./ (mean(AUC_all(:,1)));
AUC_mean = mean(AUC_norm);
AUC_std = std(AUC_norm);

% plot bar graph
b = bar(AUC_mean, 'edgecolor', 'none'); hold all;
b.FaceColor = 'flat';
b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = tay4_color;

% draw rectangle for WT because can't change edges for bars
% get x position and height
x = b.XData(1);
h = AUC_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, 'FaceColor', 'none');

% plot errorbars
er = errorbar(1:2, AUC_mean, AUC_std, 'linestyle', 'none', 'color', [0,0,0]);

% add scatter plot
for k = 1:2
    scatter(repelem(k,16), AUC_norm(:,k), [],'black', 'HandleVisibility','off');
end

[h,p] = ttest2(AUC_all(:,1), AUC_all(:,2), "Tail", "right");
ylim([0 1.5])
xticklabels([])

% save figure
w=4.5;
h=1.75;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/supp/phage_treated_mono_res");
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s4b
% solid isolation of Tay4-treated Ef from community

% load 96-well plate OD600 growth curves
load('../1-data/fig2_250325_res_Ef_solid_isolation.mat', 'res_Ef_solid_isolation')

figure
subplot(1,2,1)
AUC_all = [];
count = 1;

for c = 5:6
    for r = 1:8

        % plot growth curves / plot control first
        well = strcat(rows(r), cols(c));
        plot(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank, 'color', [0, 0, 0]); hold all;
        AUC_all(count,1) = trapz(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank);

        % plot +Tay4
        well = strcat(rows(r), cols(c+6));
        plot(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank, 'color', tay4_color);
        AUC_all(count,2) = trapz(res_Ef_solid_isolation.(well).time, res_Ef_solid_isolation.(well).OD600blank);

        count = count + 1;
    end
end

% figure details
xlim([0 48])
set(gca, 'yscale', 'log')
ylim([5*10^-4 2])
yticks([10^-2 10^0])
yticklabels({'10^{-2}', '10^0'})

% plot AUC
subplot(1,2,2)
AUC_norm = AUC_all ./ (mean(AUC_all(:,1)));
AUC_mean = mean(AUC_norm);
AUC_std = std(AUC_norm);

% plot bar graph
b = bar(AUC_mean, 'edgecolor', 'none'); hold all;
b.FaceColor = 'flat';
b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = tay4_color;

% draw rectangle for WT because can't change edges for bars
% get x position and height
x = b.XData(1);
h = AUC_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, 'FaceColor', 'none');

% plot errorbars
er = errorbar(1:2, AUC_mean, AUC_std, 'linestyle', 'none', 'color', [0,0,0]);

% add scatter plot
for k = 1:2
    scatter(repelem(k,16), AUC_norm(:,k), [],'black', 'HandleVisibility','off');
end

[h,p] = ttest2(AUC_all(:,1), AUC_all(:,2), "Tail", "right");
ylim([0 1.5])
xticklabels([])

% save figure
w=4.5;
h=1.75;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/supp/phage_treated_SIC_res");
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% figure 2a
% Tay4 PFUS/mL during phage treatment in monoculture and SIC

% load PFU data
t = readtable('../1-data/CFUs_PFUs_data/fig2_250604_tay4_treatment_PFU_quantification.xlsx');

fig = figure;
color_exp = [150,150,150] ./ 256;
lw = 1.;

% plot PFUs when Ef monoculture is treated with Tay4
d_mono = [];
idx = strcmp(t.strain, 'TT2+tay4'); % ef monoculture: +tay4, 4 recovery
temp = t(idx,:);

count = 1;
for pass = 1:5
    for rep = 1:3
        d_mono(pass,rep) = temp.tay4PFU_mL(count);
        count = count + 1;
    end
end

% get ++tay4 passage 5s
idx = strcmp(t.strain, 'TT2++tay4'); % ef monoculture: +tay4, 3 recovery, +tay4
temp_tay4 = t(idx,:);
d_mono_tay4 = d_mono;
d_mono_tay4(5,:) = temp_tay4.tay4PFU_mL;

% get expected pfus
d_mono_exp = d_mono(1,:);
for p = 2:5
    d_mono_exp(p,:) = d_mono_exp(p-1,:) ./ 200;
end

% get average across 3 replicates
d_mono_avg = mean(d_mono, 2);
d_mono_std = std(d_mono, [], 2);
d_mono_tay4_avg = mean(d_mono_tay4, 2);
d_mono_tay4_std = std(d_mono_tay4, [], 2);
d_mono_exp_avg = mean(d_mono_exp, 2);
d_mono_exp_std = std(d_mono_exp, [], 2);

% plot monoculture
subplot(1,2,1)

% plot monoculture - expected
errorbar(1:5, d_mono_exp_avg, d_mono_exp_std, 'o-', 'color', color_exp, ...
    'MarkerFaceColor', "white", 'linewidth', lw); hold all;

% plot monoculture - data
errorbar(1:5, d_mono_avg, d_mono_std, 'o-', 'color', tay4_color, ...
    'MarkerFaceColor', "white", 'linewidth', lw); hold all;

plot(1, d_mono_avg(1), 'o-', 'color', tay4_color, ...
    'MarkerFaceColor', tay4_color, 'linewidth', lw); % fill passage 1 tay4 treatment

yscale('log')
ylim([10^0 10^10])
xlim([0.7, 5.3])
yline(10^3, '--')

% plot PFUs when MD2 is treated with Tay4
d_preMDm2 = [];
idx = strcmp(t.strain, 'preMDm2+tay4'); % ef monoculture: +tay4, 4 recovery
temp = t(idx,:);

count = 1;
for pass = 1:5
    for rep = 1:3
        d_preMDm2(pass,rep) = temp.tay4PFU_mL(count);

        count = count + 1;
    end
end

% get ++tay4 passage 5s
idx = strcmp(t.strain, 'preMDm2++tay4'); % ef monoculture: +tay4, 3 recovery, +tay4
temp_tay4 = t(idx,:);
d_preMDm2_tay4 = d_preMDm2;
d_preMDm2_tay4(5,:) = temp_tay4.tay4PFU_mL;

% get expected pfus
d_preMDm2_exp = d_preMDm2(1,:);
for p = 2:5
    d_preMDm2_exp(p,:) = d_preMDm2_exp(p-1,:) ./ 200;
end

% get average across 3 replicates
d_preMDm2_avg = mean(d_preMDm2, 2);
d_preMDm2_std = std(d_preMDm2, [], 2);
d_preMDm2_tay4_avg = mean(d_preMDm2_tay4, 2);
d_preMDm2_tay4_std = std(d_preMDm2_tay4, [], 2);
d_preMDm2_exp_avg = mean(d_preMDm2_exp, 2);
d_preMDm2_exp_std = std(d_preMDm2_exp, [], 2);

subplot(1,2,2)
% plot monoculture - expected
errorbar(1:5, d_preMDm2_exp_avg, d_preMDm2_exp_std, 'o-', 'color', color_exp, ...
    'MarkerFaceColor', "white", 'linewidth', lw); hold all;

% plot monoculture - data
errorbar(1:5, d_preMDm2_avg, d_preMDm2_std, 'o-', 'color', tay4_color, ...
    'MarkerFaceColor', "white", 'linewidth', lw); hold all;

plot(1, d_preMDm2_avg(1), 'o-', 'color', tay4_color, ...
    'MarkerFaceColor', tay4_color, 'linewidth', lw); % fill passage 1 tay4 treatment

yscale('log')
ylim([10^0 10^10])
xlim([0.7, 5.3])
yline(10^3, '--')

% save figure
w=4;
h=3;
set(fig, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h]);

fname = "../3-figs/supp/Tay4_PFUs_during_growth";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s4c
% plaque images

%% supp: figure s4d
% Tay4 PFUs/mL when treated exponential vs. stationary

t = readtable('../1-data/CFUs_PFUs_data/fig2_250716_demonstrating_inability_to_infect_stationary_phase.xlsx');

fig = figure;
color_exp = [150,150,150] ./ 256;
lw = 1.;

d_input = [];
idx = strcmp(t.strain, 'media_only+tay4'); % ef monoculture: +tay4, 4 recovery
d_input = t(idx,:);
d_input = d_input(1:3,:);

d_0h = [];
idx = strcmp(t.strain, 'TT2+tay4_t=0'); % ef monoculture: +tay4, 4 recovery
d_0h = t(idx,:);

d_24h = [];
idx = strcmp(t.strain, 'TT2+tay4_t=24'); % ef monoculture: +tay4, 4 recovery
d_24h = t(idx,:);

data2plot = [];
data2plot = [d_input.tay4PFU_mL, d_0h.tay4PFU_mL, d_24h.tay4PFU_mL];
b = bar(mean(data2plot), 'FaceColor', 'white', 'EdgeColor', 'black', ...
    'LineWidth', 1.3);
hold all;
er = errorbar(mean(data2plot), std(data2plot), 'linestyle', 'none', ...
    'color', 'black', 'LineWidth', 1.3);
hold all;

for i = 1:3
    scatter(repelem(i,3), data2plot(:,i), 'black');
end

yscale('log')
ylim([10^0, 10^12])
xlim([0.25, 3.75])
xticklabels([])

% get p-values
[h, p] = ttest2(d_input.tay4PFU_mL, d_0h.tay4PFU_mL, "tail", "left")
[h, p] = ttest2(d_input.tay4PFU_mL, d_24h.tay4PFU_mL, "tail", "left")

% save figure
w=1.75;
h=2;
set(fig, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h]);

fname = "../3-figs/supp/pfus_stationary";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% fig. 2c
% plot resistance phenotype when isolated with or without phage pressure

% load 96-well plate OD600 growth curves
% p1: isolation without Tay4 pressure
% p2: isolation with Tay4 pressure
% p3: second passage  - test if isolates from p1 and p2 are still resistant
% to Tay4
load('../1-data/fig2_250410_res_Ef_liquid_isolation.mat', 'res_Ef_liquid_isolation')

figure
cols2plot = [1,4]; % plot +-+ first, then +++
reps2plot = [8,7];
AUC_all = [];
p3 = res_Ef_liquid_isolation(3);

for i = 1:2
    subplot(1,2,i)
    c = cols2plot(i);
    n = reps2plot(i);
    for r = 1:n

        % plot control first
        well = strcat(rows(r), cols(c));
        plot(p3.(well).time, p3.(well).OD600blank, 'color', [0,0,0], ...
            'linewidth', 1.);
        hold all;

        % plot tay4-treated
        well = strcat(rows(r), cols(c+6))
        plot(p3.(well).time, p3.(well).OD600blank, 'color', tay4_color, ...
            'linewidth', 1.3);
    end

    % figure details
    set(gca, 'yscale', 'log')
    ylim([2*10^-3 2])
    xlim([0 10])
    yticks([10^-2 10^0])
    yticklabels({'10^{-2}', '10^0'})
end

% save figure
w=4;
h=1.5;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = "../3-figs/fig2/liquid_isol_phage_pressure";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s5a
% Ef genomic point mutations
% see code: 

%% figure 2d
% plot hyx inverted relative abundance

% load DELLY data
d = readtable('../1-data/DELLY/fig2_DELLY_variant_freq_hyx.txt');

% 10 monoculture, -tay4
% 19 monoculture, +tay4
% 14 SIC, -tay4
% 26 SIC, +tay4

% just for hyx data
d_mono = d(strcmp(d.SUBSET, 'mono –Tay4 pooled'), :);
d_mono_tay4 = d(strcmp(d.SUBSET, 'mono +Tay4 pooled'), :);
d_SIC = d(strcmp(d.SUBSET, 'SIC –Tay4'), :);
d_SIC_tay4 = d(strcmp(d.SUBSET, 'SIC +Tay4'), :);

% get data for hyx inversion
xdata = (repelem([1, 2, 3, 4], ...
    [height(d_mono), height(d_mono_tay4), height(d_SIC), height(d_SIC_tay4)]));
ydata = d.AF;

figure
sc(1) = swarmchart(repelem(1, height(d_mono)), d_mono.AF, 'filled', 'MarkerFaceColor', 'black');
hold on
sc(2) = swarmchart(repelem(1.75, height(d_mono_tay4)), d_mono_tay4.AF, 'filled', 'MarkerFaceColor', tay4_color);
sc(3) = swarmchart(repelem(3, height(d_SIC)), d_SIC.AF, 'filled', 'MarkerFaceColor', 'black');
sc(4) = swarmchart(repelem(3.75, height(d_SIC_tay4)), d_SIC_tay4.AF, 'filled', 'MarkerFaceColor', tay4_color);

% plot details for swarmchart
alpha_vals = [0.5, 0.6, 0.5, 0.6];
for j = 1:4
    sc(j).MarkerFaceAlpha = alpha_vals(j);
    sc(j).XJitter = 'density';
    sc(j).XJitterWidth = 0.4;
    sc(j).SizeData = 20;
    sc(j).MarkerEdgeColor = 'none';
end

% boxplot
bc(1) = boxchart(repelem(1, height(d_mono)), d_mono.AF, 'BoxFaceColor', 'black');
bc(2) = boxchart(repelem(1.75, height(d_mono_tay4)), d_mono_tay4.AF, 'BoxFaceColor', tay4_color);
bc(3) = boxchart(repelem(3, height(d_SIC)), d_SIC.AF, 'BoxFaceColor', 'black');
bc(4) = boxchart(repelem(3.75, height(d_SIC_tay4)), d_SIC_tay4.AF, 'BoxFaceColor', tay4_color);

alpha_vals = [0.2, 0.2, 0.2, 0.2];
for k = 1:4
    bc(k).MarkerStyle = 'none';
    bc(k).BoxEdgeColor = 'black';
    bc(k).BoxFaceAlpha = alpha_vals(k);
end

% other figure details
set(gca, 'FontName', 'Arial')
ylim([-0.05 1.15])
xticklabels([])
xticks([])
box on 

p = [];
p(1) = ranksum(d_mono.AF, d_mono_tay4.AF);
p(2) = ranksum(d_SIC.AF, d_SIC_tay4.AF);

% save figure
w=2.5;
h=2.5;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], ... 
    'PaperUnits', 'Inches', 'PaperSize', [w, h], 'Renderer', 'painters')

fname = "../3-figs/fig3/DELLY_hyx";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s5b
% plot ALL structural variations in Ef isolates

% load DELLY data
d_all = readtable('../1-data/DELLY/fig2_DELLY_variant_freq_all.txt');

% 10 monoculture, -tay4
% 19 monoculture, +tay4
% 14 SIC, -tay4
% 26 SIC, +tay4

% all the variant frequencies
d_all_mono = d_all(strcmp(d_all.SUBSET, 'mono –Tay4'), :);
d_all_mono_tay4 = d_all(strcmp(d_all.SUBSET, 'mono +Tay4'), :);
d_all_SIC = d_all(strcmp(d_all.SUBSET, 'SIC –Tay4'), :);
d_all_SIC_tay4 = d_all(strcmp(d_all.SUBSET, 'SIC +Tay4'), :);

genes = unique(d_all.TRA_COORD);

figure
p_all = [];
for i = 1:8
    gene2plot = genes{i};
    dd = d_all(strcmp(d_all.TRA_COORD, gene2plot), :); % get all data for specific coordinate

    dd_mono = dd(strcmp(dd.SUBSET, 'mono –Tay4'), :);
    dd_mono_tay4 = dd(strcmp(dd.SUBSET, 'mono +Tay4'), :);
    dd_SIC = dd(strcmp(dd.SUBSET, 'SIC –Tay4'), :);
    dd_SIC_tay4 = dd(strcmp(dd.SUBSET, 'SIC +Tay4'), :);

    subplot(2,4,i)
    % get data for hyx inversion
    xdata = (repelem([1, 2, 3, 4], ...
        [height(dd_mono), height(dd_mono_tay4), height(dd_SIC), height(dd_SIC_tay4)]));
    ydata = dd.AF;

    sc(1) = swarmchart(repelem(1, height(dd_mono)), dd_mono.AF, 'filled', 'MarkerFaceColor', 'black');
    hold on
    sc(2) = swarmchart(repelem(1.75, height(dd_mono_tay4)), dd_mono_tay4.AF, 'filled', 'MarkerFaceColor', tay4_color);
    sc(3) = swarmchart(repelem(3, height(dd_SIC)), dd_SIC.AF, 'filled', 'MarkerFaceColor', 'black');
    sc(4) = swarmchart(repelem(3.75, height(dd_SIC_tay4)), dd_SIC_tay4.AF, 'filled', 'MarkerFaceColor', tay4_color);

    % plot details for swarmchart
    alpha_vals = [0.5, 0.6, 0.5, 0.6];
    for j = 1:4
        sc(j).MarkerFaceAlpha = alpha_vals(j);
        sc(j).XJitter = 'density';
        sc(j).XJitterWidth = 0.4;
        sc(j).SizeData = 20;
        sc(j).MarkerEdgeColor = 'none';
    end

    % boxplot
    bc(1) = boxchart(repelem(1, height(dd_mono)), dd_mono.AF, 'BoxFaceColor', 'black'); 
    bc(2) = boxchart(repelem(1.75, height(dd_mono_tay4)), dd_mono_tay4.AF, 'BoxFaceColor', tay4_color); 
    bc(3) = boxchart(repelem(3, height(dd_SIC)), dd_SIC.AF, 'BoxFaceColor', 'black'); 
    bc(4) = boxchart(repelem(3.75, height(dd_SIC_tay4)), dd_SIC_tay4.AF, 'BoxFaceColor', tay4_color); 
    
    alpha_vals = [0.2, 0.2, 0.2, 0.2];
    for k = 1:4
        bc(k).MarkerStyle = 'none';
        bc(k).BoxEdgeColor = 'black';
        bc(k).BoxFaceAlpha = alpha_vals(k);
    end

    p_all(i,1) = ranksum(dd_mono.AF, dd_mono_tay4.AF);
    p_all(i,2) = ranksum(dd_SIC.AF, dd_SIC_tay4.AF);

    % other figure details
    set(gca, 'FontName', 'Arial')
    ylim([-0.05 1.15])
    xticklabels([])
    title(gene2plot)
    xticks([])
    box on 
end

% save figure
w=10;
h=4.5;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h]);

fname = "../3-figs/fig3/DELLY_all";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%%