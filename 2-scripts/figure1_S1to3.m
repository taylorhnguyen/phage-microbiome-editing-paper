%% script to generate figure 1 and supp figures 1-3
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

%% figure 1a
% imaging

%% figure 1b
% Ef monoculture + Tay4
% pass 1 (+Tay4), pass 2 (recovery), pass 3 (+Tay4)

% load 96-well plate OD600 growth curves
load('../1-data/fig1_230524_OD_monoculture_1X_MOI.mat', 'mono_1X_MOI')

% EDIT HERE TO PLOT DIFFERENT PASSAGE
%pass = 1;
pass = 3;

% plot growth curves first
figure
subplot(1,2,1)

AUC2plot = [];
for r = 1:8
        % plot controls
        well = strcat(rows(r),cols(1));
        plot(mono_1X_MOI(pass).(well).time, mono_1X_MOI(pass).(well).OD600blank, 'Linewidth', 1, 'Color', [0,0,0]./255); hold all;
        AUC2plot(r,1) = trapz(mono_1X_MOI(pass).(well).time, mono_1X_MOI(pass).(well).OD600blank);
        
        % plot +phage
        well = strcat(rows(r),cols(3)); % column 3: tay4
        plot(mono_1X_MOI(pass).(well).time, mono_1X_MOI(pass).(well).OD600blank, 'Linewidth', 1, 'Color', tay4_color); hold all;
        AUC2plot(r,2) = trapz(mono_1X_MOI(pass).(well).time, mono_1X_MOI(pass).(well).OD600blank);
end

set(gca,'yscale','log')
xlim([0 48])
ylim([10^-3 1])
xlabel('time (h)');
ylabel('blanked OD600')
axis square

% plot AUC
subplot(1,2,2)
AUC2plot = AUC2plot ./ (mean(AUC2plot(:,1)));
AUC_mean = mean(AUC2plot);
AUC_std = std(AUC2plot);
b = bar(AUC_mean, 'FaceColor', 'flat', 'LineWidth', 1.3, 'linestyle', 'none'); hold all;

b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = tay4_color;

% draw rectangle for WT because can't change edges for bars
% get x position and height
x = b.XData(1);
h = AUC_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, ...
    'FaceColor', 'none');

% add errorbars
er = errorbar(AUC_mean, AUC_std); hold all;
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% add scatter plot
for k = 1:2
    scatter(repelem(k,8),AUC2plot(:,k), [],'black', 'HandleVisibility','off'); hold all;
end

axis square
set(gca, 'XTickLabel', {'WT', 'Tay4'})
ylabel('normalized AUC');
ylim([0 1.3])
yticks([0 0.2 0.4 0.6 0.8 1.0])
xlim([0.5 2.5])

p = ranksum(AUC2plot(:,1), AUC2plot(:,2), "tail", "right")

% save figure
w=4;
h=3;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/fig1/ef_mono_tay4_pass", string(pass));
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s1a
% Ef monoculture + Tay4
% pass 1 (+Tay4), pass 2 (recovery), pass 3 (+Tay4)

% load 96-well plate OD600 growth curves
load('../1-data/fig1_230524_OD_monoculture_1X_MOI.mat', 'mono_1X_MOI')

% EDIT HERE TO PLOT DIFFERENT PASSAGE
%pass = 1;
pass = 3;

% params:
% column (3 for Tay4), passage, dataframe
p_vals = plot_phage_growth_comparison(3, pass, mono_1X_MOI)

% save figure
w=10;
h=10;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat(strcat("../3-figs/supp/ef_mono_tay4_growth_pass", string(pass)));
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s1b
% plot heat-killed controls

% load 96-well plate OD600 growth curves
load('../1-data/fig1_230524_OD_monoculture_heat_killed.mat', 'mono_heat_killed')

figure
% plot control wells
ctrl_wells = ["B1", "C1", "D1", "E1", "B12", "C12", "D12", "E12"];
OD2keep = [];
for i = 1:8
    well = ctrl_wells(i);
    OD2keep = [OD2keep, mono_heat_killed(1).(well).OD600blank];
end

OD_avg = mean(OD2keep,2);
OD_std = std(OD2keep,1,2);

t = mono_heat_killed(1).("A1").time;
plot(t, OD_avg, 'color', 'black', 'linewidth', 1); 
hold all;

% plot standard deviation shaded region
y_upper = OD_avg + OD_std;
y_lower = OD_avg - OD_std;

patch([t', fliplr(t')], [y_lower', fliplr(y_upper')], 'black', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');

% plot 1 uL highest amount of phage
c = 10;
OD2keep = [];
for r = 1:8
    well = strcat(rows(r),cols(c));
    OD2keep = [OD2keep, mono_heat_killed(1).(well).OD600blank];
end
OD_avg = mean(OD2keep,2);
OD_std = std(OD2keep,1,2);

plot(t, OD_avg, 'Color', tay4_color); 
hold all;

% plot standard deviation shaded region
y_upper = OD_avg + OD_std;
y_lower = OD_avg - OD_std;

p = patch([t', fliplr(t')], [y_lower', fliplr(y_upper')], 'black',...
    'FaceAlpha', 0.4, 'EdgeColor', 'none');
p.FaceColor = 'flat';
p.FaceColor = tay4_color;

% figure details
ylim([0 0.8])
xlim([0 48])
xlabel('Time (h)')
ylabel('Blanked OD600')
set(gca,'FontSize',14)

% save figure
w=2.5;
h=2;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/supp/tay4_heatkilled");
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s1c
% compare Tay4 treatment MOI 1 vs. 10

% load 96-well plate OD600 growth curves
load('../1-data/fig1_230524_OD_monoculture_1X_MOI.mat', 'mono_1X_MOI')
load('../1-data/fig1_230524_OD_monoculture_10X_MOI.mat', 'mono_10X_MOI')

% EDIT HERE TO PLOT DIFFERENT PASSAGE 
%pass = 1;
pass = 3;

% plot passage 1 for MOI 1 and 10 together 
figure

for r = 1:8 % get all 8 biological replicates

        % plot controls
        well = strcat(rows(r),cols(1));
        AUC2plot(r,1) = trapz(mono_1X_MOI(pass).(well).time, mono_1X_MOI(pass).(well).OD600blank);
        
        % plot +phage MOI~1
        well = strcat(rows(r),cols(3));
        AUC2plot(r,2) = trapz(mono_1X_MOI(pass).(well).time, mono_1X_MOI(pass).(well).OD600blank);

        % plot +phage MOI~10
        well = strcat(rows(r),cols(3));
        AUC2plot(r,3) = trapz(mono_10X_MOI(pass).(well).time, mono_10X_MOI(pass).(well).OD600blank);
end

AUC2plot = AUC2plot ./ (mean(AUC2plot(:,1)));
AUC_mean = mean(AUC2plot);
AUC_std = std(AUC2plot);
b = bar(AUC_mean, 'FaceColor', 'flat', 'LineWidth', 1.3, 'linestyle', 'none');
hold on

b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = tay4_color;
b.CData(3,:) = tay4_color;

% draw rectangle for WT because can't change edges for bars
% get x position and height
x = b.XData(1);
h = AUC_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, ...
    'FaceColor', 'none');

% add errorbars
er = errorbar(AUC_mean, AUC_std); hold all;
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% add scatter plot
for k = 1:3
    scatter(repelem(k,8),AUC2plot(:,k), [],'black', 'HandleVisibility','off'); hold all;
end

ylim([0 1.5])
xlim([0.2 3.8])
yticks([0 0.2 0.4 0.6 0.8 1.0])

% get p-values
p1 = [];
p1(1) = ranksum(AUC2plot(:,1), AUC2plot(:,2), "tail", "right");
p1(2) = ranksum(AUC2plot(:,1), AUC2plot(:,3), "tail", "right");
p1(3) = ranksum(AUC2plot(:,2), AUC2plot(:,3));
p1

% save figure
w=2.5;
h=1.75;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/supp/tay4_MOI_1_vs_10_pass", string(pass));
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s1d
% plaque photos

%% supp: figure s1e
% testing resistance at passage 1 and 2

% load 96-well plate OD600 growth curves
% p1: passage 1
% p2: passage 2
load('../1-data/fig1_250402_Ef_mono_resistance.mat', 'mono_resistance')

% EDIT HERE TO CHANGE PASSAGE
%pass = 1;
pass = 2;

figure
AUC2plot = [];

for r = 1:8
    % control 
    well = strcat(rows(r), cols(2));
    plot(mono_resistance(pass).(well).time, mono_resistance(pass).(well).OD600blank, 'black', 'LineWidth', 1.3); hold all;
    AUC2plot(r,1) = trapz(mono_resistance(pass).(well).time, mono_resistance(pass).(well).OD600blank);

    % plot Tay4-treated
    well = strcat(rows(r), cols(8));
    plot(mono_resistance(pass).(well).time, mono_resistance(pass).(well).OD600blank, 'color', tay4_color, 'LineWidth', 1.3); hold all;
    AUC2plot(r,2) = trapz(mono_resistance(pass).(well).time, mono_resistance(pass).(well).OD600blank);
end

set(gca, 'YScale', 'log')
ylim([10^-3 10^0])
xlim([0 48])

w=2.5;
h=1.7;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = "../3-figs/supp/growth";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% figure 1c
% imaging

%% supp: figure s2a
% plot isolates treated with Tay4 & cocktail

% load 96-well plate OD600 growth curves
% p1: isols + Tay4
% p2: isols + cocktail
load('../1-data/fig1_250611_15mb_isols_phage.mat', 'isols_phage')
p1 = isols_phage(1);
p2 = isols_phage(2);

% order to plot
species_names_ordered = ...
    ["E. fergusonii", "K. pneumoniae", "E. hirae", "E. faecium", ...
     "E. faecalis", "B. theta", "B. fragilis", "B. uniformis", ...
     "P. distasonis", "F. plautii", "C. symbiosum", "C. clostridioforme", ...
     "B. producta", "C. scindens", "C. hylemonae", "C. hathewayi"];

sp_ordered = [3, 1, 11, 6, ...
              12, 13, 14, 4, ...
              15, 5, 16, 7, ...
              10, 8, 9, 2];

fig = figure;
count = 1;
t = p1.('A1').time;

AUC_ctrl = [];
AUC_tay4 = [];
AUC_cocktail = [];

for i = [1, 4] % go through the first 8 species, then second 8 species

    for r = 1:8

        % +TAY4 / get all data in one matrix
        OD_all = [];
        for c = 1:12 % get all columns
            well = strcat(rows(r), cols(c));
            OD_all(:,c) = p1.(well).OD600blank;
        end

        % get mean of replicates
        OD_ctrl_avg = mean(OD_all(:, i:i+2), 2);
        OD_ctrl_std = std(OD_all(:, i:i+2), [], 2);

        OD_tay4_avg = mean(OD_all(:, i+6:i+6+2), 2);
        OD_tay4_std = std(OD_all(:, i+6:i+6+2), [], 2);

        subplot(4,4, sp_ordered(count))

        % plot controls & standard deviation shaded region
        plot(t, OD_ctrl_avg, 'Color', 'black');
        hold all;
        y_upper = OD_ctrl_avg + OD_ctrl_std;
        y_lower = OD_ctrl_avg - OD_ctrl_std;

        p = patch([t', fliplr(t')], [y_lower', fliplr(y_upper')], 'black',...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');

        % plot +Tay4 & standard deviation shaded region
        plot(t, OD_tay4_avg, 'Color', tay4_color);
        y_upper = OD_tay4_avg + OD_tay4_std;
        y_lower = OD_tay4_avg - OD_tay4_std;

        p = patch([t', fliplr(t')], [y_lower', fliplr(y_upper')], 'black',...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
        p.FaceColor = 'flat';
        p.FaceColor = tay4_color;

        % +COCKTAIL / get all data in one matrix
        OD_all = [];
        for c = 1:12 % get all columns
            well = strcat(rows(r), cols(c));
            OD_all(:,c) = p2.(well).OD600blank;
        end

        % get mean of replicates
        OD_cocktail_avg = mean(OD_all(:, i+6:i+6+2), 2);
        OD_cocktail_std = std(OD_all(:, i+6:i+6+2), [], 2);

        % plot +cocktail & standard deviation shaded region
        plot(t, OD_cocktail_avg, 'Color', cocktail_color);
        y_upper = OD_cocktail_avg + OD_cocktail_std;
        y_lower = OD_cocktail_avg - OD_cocktail_std;

        p = patch([t', fliplr(t')], [y_lower', fliplr(y_upper')], 'black',...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
        p.FaceColor = 'flat';
        p.FaceColor = cocktail_color;

        count = count+1;
    end
end

for i = 1:16
    subplot(4,4,i)

    % figure details
    title(species_names_ordered(i));
    ylim([0 1.2])
    xlim([0 48])
end

w=8;
h=6;
set(fig, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h], ...
    'Renderer', 'painters');

fname = "../3-figs/isols_tay4_cocktail_ordered";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% supp: figure s2b
% plot lag time differences in isols + Tay4/cocktail

% load 96-well plate OD600 growth curves
% p1: isols + Tay4
% p2: isols + cocktail
load('../1-data/fig1_250611_15mb_isols_phage.mat', 'isols_phage')

% edit here to plot lag time comparison between tay4 or cocktail
p1 = isols_phage(1); % use for tay4
%p1 = isols_phage(2); % use for cocktail

% order to plot
species_names_ordered = ...
    ["E. fergusonii", "K. pneumoniae", "E. hirae", "E. faecium", ...
     "E. faecalis", "B. theta", "B. fragilis", "B. uniformis", ...
     "P. distasonis", "F. plautii", "C. symbiosum", "C. clostridioforme", ...
     "B. producta", "C. scindens", "C. hylemonae", "C. hathewayi"];

sp_ordered = [3, 1, 11, 6, ...
              12, 13, 14, 4, ...
              15, 5, 16, 7, ...
              10, 8, 9, 2];

lag_ctrl = [];
lag_phage = [];
count = 1;

% get control first
% lag time: time to half max OD
for c = [1, 4] % go through the first 8 species, then second 8 species
    for r = 1:8
        for rep = 1:3

            % 1: get control lag time
            well = strcat(rows(r), cols(c + rep - 1));
            lag_ctrl(sp_ordered(count), rep) = p1.(well).lag_time;

            % 2: get +phage lag time
            well = strcat(rows(r), cols(c + rep + 5));
            lag_phage(sp_ordered(count), rep) = p1.(well).lag_time;
        end

        count = count + 1;
    end
end

% get stats
mean_ctrl = mean(lag_ctrl, 2);
mean_phage = mean(lag_phage, 2);
sd_ctrl = std(lag_ctrl, [], 2);
sd_phage = std(lag_phage, [], 2);

% plot bar graph
fig = figure
b = bar([mean_ctrl, mean_phage], 'grouped', 'groupwidth', 0.7);
hold all;
b(1).FaceColor = "white";
b(1).EdgeColor = "black";
b(2).FaceColor = tay4_color;
%b(2).FaceColor = cocktail_color;

% plot errorbar
x_ctrl = b(1).XEndPoints;
x_phage = b(2).XEndPoints;

er(1) = errorbar(x_ctrl, mean_ctrl, sd_ctrl, ...
    'linestyle', 'none', 'color', 'black', 'linewidth', 1.);
er(1) = errorbar(x_phage, mean_phage, sd_phage, ...
    'linestyle', 'none', 'color', 'black', 'linewidth', 1.);

% add individual datapoints
jitterWidth = 0.15;
alphaVal = 0.4;

for i = 1:16
    % control points
    jitter_x = x_ctrl(i) + (rand(1, size(lag_ctrl, 2)) - 0.5) * jitterWidth;
    scatter(jitter_x, lag_ctrl(i, :), 20, 'filled', 'MarkerFaceAlpha', alphaVal, ...
        'markerfacecolor', 'black');
    
    % +Tay4 points
    jitter_x = x_phage(i) + (rand(1, size(lag_phage, 2)) - 0.5) * jitterWidth;
    scatter(jitter_x, lag_phage(i, :), 20, 'k', 'filled', 'MarkerFaceAlpha', alphaVal, ...
        'markerfacecolor', 'black');
end

% get p_vals
p_vals = [];
for i = 1:16
    [h, p_vals(i)] = ttest2(lag_ctrl(i,:), lag_phage(i,:));
end

% figure details
ylim([0 40])
xticks([])

% save figure
w=9;
h=2;
set(fig, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h]);
fname = "../3-figs/isols_tay4_cocktail_ordered_lagtime_tay4";
%fname = "../3-figs/isols_tay4_cocktail_ordered_lagtime_cocktail";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% figure 1d
% SICs + Tay4 (MOI~1) for passages 1-6
% plot Ef and Kp relative abundance

% load 16s ASV relative abundance table
load("../1-data/fig1_240108_SIC_Ef_resistance_ASV_table.mat", "SIC_res_ASV_relab_table")

% plot average of Ef (top panel) and Kp (bottom panel) 
% plot as 1 2 3 4 5 6/6P (two lines) 
ef_idx = 1;
kp_idx = 6;

% figure parameters
SIC_shift = [0:12:84];
lw = 1.;

figure
for s = 4:7
    
    Ef_relab_ctrl = [];
    Kp_relab_ctrl = [];
    Ef_relab_phage = [];
    Kp_relab_phage = [];
    Ef_relab_phage_phage = [];
    Kp_relab_phage_phage = [];

    for r = 1:3
        idx2plot = [r+SIC_shift(s):96:576];
        Ef_relab_ctrl(:,r) = SIC_res_ASV_relab_table(idx2plot, ef_idx);
        Kp_relab_ctrl(:,r) = SIC_res_ASV_relab_table(idx2plot, kp_idx);
    
        Ef_relab_phage(:,r) = SIC_res_ASV_relab_table(idx2plot+3, ef_idx);
        Kp_relab_phage(:,r) = SIC_res_ASV_relab_table(idx2plot+3, kp_idx);

        Ef_relab_phage_phage(:,r) = SIC_res_ASV_relab_table(idx2plot+3, ef_idx);
        Kp_relab_phage_phage(:,r) = SIC_res_ASV_relab_table(idx2plot+3, kp_idx);
    end

    % rewrite the 6th passage for 6th passage + Tay4
    for r = 1:3
        idx2plot = [r+SIC_shift(s)+96*6];
        Ef_relab_phage_phage(6,r) = SIC_res_ASV_relab_table(idx2plot+3, ef_idx);
        Kp_relab_phage_phage(6,r) = SIC_res_ASV_relab_table(idx2plot+3, kp_idx);
    end

    % plot Ef (first row)
    subplot(2,4,s-3)

    % plot control
    y = mean(Ef_relab_ctrl, 2);
    sdev = std(Ef_relab_ctrl, [], 2);
    plot(1:6, y, 'o-', 'Color', 'black', 'MarkerFaceColor', 'white', 'linewidth', lw);
    hold all;
    errorbar(1:6, y, sdev, 'Color', 'black', "LineStyle","none", 'linewidth', lw);

    % plot phage
    y = mean(Ef_relab_phage, 2);
    sdev = std(Ef_relab_phage, [], 2);

    % plot untreated passages with white fill
    plot(2:6, y(2:6), 'o-', 'Color', tay4_color, 'MarkerFaceColor', "white", 'linewidth', lw); 

    % plot phage treated passages with fill
    plot(1, y(1), 'o-', 'Color', tay4_color, 'MarkerFaceColor', tay4_color, 'linewidth', lw);
    errorbar(1:6, y, sdev, 'Color', tay4_color, "LineStyle","none", 'linewidth', lw);

    % plot phage-phage
    y = mean(Ef_relab_phage_phage, 2);
    sdev = std(Ef_relab_phage, [], 2);

    % plot untreated passages with white fill
    plot(1:5, y(1:end-1), 'o-', 'Color', tay4_color, 'MarkerFaceColor', "white", 'linewidth', lw);
    plot(1, y(1), 'o', 'Color', tay4_color, 'MarkerFaceColor', tay4_color, 'linewidth', lw);
    plot(5:6, y(5:6), 'o-', 'Color', tay4_color, 'MarkerFaceColor', "white", 'linewidth', lw); % change linestyle in illustrator
    plot(6, y(6), 'o', 'Color', tay4_color, 'MarkerFaceColor', tay4_color, 'linewidth', lw);

    errorbar(1:6, y, sdev, 'Color', tay4_color, "LineStyle","none", 'linewidth', lw);

    xticklabels([[]])
    ylim([0 0.5]);
    xlim([0.5 6.5])

    % get p-values
    [h, p_vals{s}(1)] = ttest2(Ef_relab_ctrl(1,:), Ef_relab_phage(1,:)); % WT vs. +Tay4 pass 1
    [h, p_vals{s}(2)] = ttest2(Ef_relab_phage(1,:), Ef_relab_phage(6,:)); % WT vs. +Tay4 pass 6
    [h, p_vals{s}(3)] = ttest2(Ef_relab_phage(6,:), Ef_relab_phage_phage(6,:)); % +Tay4 vs. ++Tay4 pass 6

    % plot Kp
    subplot(2,4,s-3+4)

    % plot control
    y = mean(Kp_relab_ctrl, 2);
    sdev = std(Kp_relab_ctrl, [], 2);
    plot(1:6, y, 'o-', 'Color', 'black', 'MarkerFaceColor', 'white', 'linewidth', lw);
    hold all;
    errorbar(1:6, y, sdev, 'Color', 'black', "LineStyle","none", 'linewidth', lw);
    
    % plot phage
    y = mean(Kp_relab_phage, 2);
    sdev = std(Kp_relab_phage, [], 2);

    % plot untreated passages with white fill
    plot(2:6, y(2:end), 'o-', 'Color', tay4_color, 'MarkerFaceColor', "white",'linewidth', lw);

    % plot phage treated passages with fill
    plot(1, y(1), 'o-', 'Color', tay4_color, 'MarkerFaceColor', tay4_color,'linewidth', lw);   
    errorbar(1:6, y, sdev, 'Color', tay4_color, "LineStyle","none",'linewidth', lw);

    % plot phage-phage
    y = mean(Kp_relab_phage_phage, 2);
    sdev = std(Kp_relab_phage_phage, [], 2);

    % plot phage treated passages with fill
    plot(1:5, y(1:5), 'o-', 'Color', tay4_color, 'MarkerFaceColor', "white",'linewidth', lw);
    plot(1, y(1), 'o-', 'Color', tay4_color, 'MarkerFaceColor', tay4_color,'linewidth', lw);
    plot(5:6, y(5:6), 'o-', 'Color', tay4_color, 'MarkerFaceColor', "white",'linewidth', lw); % change linestyle in illustrator
    plot(6, y(6), 'o-', 'Color', tay4_color, 'MarkerFaceColor', tay4_color,'linewidth', lw);

    errorbar(1:6, y, sdev, 'Color', tay4_color, "LineStyle","none",'linewidth', lw);
    
    xticklabels([])
    xlim([0.5 6.5])
    yline([0 0.1])

end

p_vals

w=7;
h=3;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = "../3-figs/fig1/Entero_ASV_relab_tay4";
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% figure 1e
% compare ef and kp metabolomics

% load untargeted metabolomics data
% run "/2-scripts/metab_generate_niches.m" to generate dataframe
load('../1-data/metabolomics_isolates_consumed_all.mat', 'metab_consumed')

names = {"TT1_", "TT2_", "TT3_", "TT4_", "TT5_", "TT6_", ...
         "TT7_", "TT8_", "TT9_", "TT10_", "TT11_", "TT12", ...
         "TT13_", "TT14_", "TT15_", "Klebsiella", "Salmonella", "MG1655"};

% make a venn diagram between ef and kp
metab_consumed_filt = [];
metab_consumed_filt = metab_consumed(:,[2, 16]);
metab_consumed_filt(:,3) = metab_consumed_filt(:,1)==1 & metab_consumed_filt(:,2)==1;
t = sum(metab_consumed_filt);

% pct. unique to kp
(t(2)-t(3))/t(2)

ef_unique = t(1)-t(3)
kp_unique = t(2)-t(3)
ef_kp_shared = t(3)

figure
vals = [ef_unique, kp_unique, ef_kp_shared];

venn(vals)

% save figure
w=4;
h=2;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/fig1/venn_efvskp")
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% figure 1f
% plot Ef+Kp pairwise +/- Tay4

% load 16S data
load("../1-data/fig1_240701_Ef_Kp_pairwise_ASV_relab.mat", "Ef_Kp_pairwise_ASV_relab"); % relative abundance data
load("../1-data/fig1_240701_Ef_Kp_pairwise_ASV_raw_reads.mat", "Ef_Kp_pairwise_raw_reads"); % all raw reads (organized by sample number listed in NCBI sample name)
load("../1-data/fig1_240701_Ef_Kp_ASV_table.mat", "Ef_Kp_pairwise_ASV_taxa") % DADA2 ASV names

% 1:10 ef:kp +/- tay4
% line plots

SICs2plot = [1:12:48];  % +Kp (1:100) rep1, rep2, +Kp (1:100) +Tay4, rep2
SICS2plot_phage = [7:12:48]; % MOI=1

color2plot_Ef = [241,125,177] ./ 255;
color2plot_Kp = [0, 174, 239] ./ 255;
colors2plot = [color2plot_Ef; color2plot_Kp];

ef_idx = 1;
kp_idx = 2;
   
figure
for r = 1:1 % iterate through the different ratios
    ef_data2plot = [];
    kp_data2plot = [];
    ef_phage_data2plot = [];
    kp_phage_data2plot = [];

    for rep = 1:3 % iterate to get all the reps
        idx2plot = [SICs2plot(r)+rep-1:48:48*6]; % get all 6 passages
        ef_data2plot(:,rep) = Ef_Kp_pairwise_ASV_relab(idx2plot, ef_idx);
        kp_data2plot(:,rep) = Ef_Kp_pairwise_ASV_relab(idx2plot, kp_idx);
        
        idx2plot = [SICS2plot_phage(r)+rep-1:48:48*6]; % get all 6 passages
        ef_phage_data2plot(:,rep) = Ef_Kp_pairwise_ASV_relab(idx2plot, ef_idx);
        kp_phage_data2plot(:,rep) = Ef_Kp_pairwise_ASV_relab(idx2plot, kp_idx);
    end

    % plot relative abundance (- phage)
    ef_avg = mean(ef_data2plot,2);
    ef_std = std(ef_data2plot,[],2);
    kp_avg = mean(kp_data2plot,2);
    kp_std = std(kp_data2plot,[],2);

    subplot(2,1,1)
    errorbar([1:6], ef_avg, ef_std, 'o-', 'color', colors2plot(1,:), ...
        'LineWidth', 1.);
    hold all;
    errorbar([1:6], kp_avg, kp_std, 'o-', 'color', colors2plot(2,:), ...
        'linewidth', 1.);

    title('1:10 Ef:Kp')
    xlabel('passages')
    ylim([0 1.1])
    yline(1,'--')
    xlim([0.5 6.5])

    % plot relative abundance (+ phage)
    ef_avg_phage = mean(ef_phage_data2plot,2);
    ef_std_phage = std(ef_phage_data2plot,[],2);
    kp_avg_phage = mean(kp_phage_data2plot,2);
    kp_std_phage = std(kp_phage_data2plot,[],2);

    data2plot_phage = [ef_avg_phage, kp_avg_phage];

    subplot(2,1,2)
    errorbar([1:6], ef_avg_phage, ef_std_phage, 'o-', 'color', colors2plot(1,:), ...
        'linewidth', 1);
    hold all;
    errorbar([1:6], kp_avg_phage, kp_std_phage, 'o-', 'color', colors2plot(2,:), ...
        'linewidth', 1);

    title('1:10 Ef:Kp + phage')
    xlabel('passages')
    ylim([0 1.1])
    yline(1,'--')
    xlim([0.5 6.5])

end

% save figure
w=2.5;
h=4;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

fname = strcat("../3-figs/fig1/Ef_Kp_pairwise_MOI1");
print(gcf, strcat(fname, '.png'),'-dpng','-r300');
print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

%% figure 1g
% niche test: 14mb SIC + Ef titration

% load 16S data
load("../1-data/fig1_210114_niche_test_15mb_ASV_relab.mat", "niche_test_15mb_ASV_relab"); % relative abundance data
load("../1-data/fig1_210114_niche_test_15mb_ASV_raw_reads.mat", "niche_test_15mb_raw_reads"); % all raw reads (organized by sample number listed in NCBI sample name)
load("../1-data/fig1_210114_niche_test_15mb_ASV_DADA2_taxa.mat", "niche_test_15mb_ASV_taxa") % DADA2 ASV names

pass_shift = [0, 96, 192];

for p = 1:3

    idx2plot = [25, 34:-1:26] + pass_shift(p);

    data2plot = niche_test_15mb_ASV_relab(idx2plot,:); % rel abundance
    asvs2plot = find(any(data2plot>0,1)); % finds the families that have relative abundance >0 to plot only these

    % order by family
    t = niche_test_15mb_ASV_taxa(asvs2plot)';
    [asvs, asvs2plot_sort] = sort(t);
    asvs = flip(asvs);
    asvs2plot_sort = flip(asvs2plot_sort);
    data2plot_filt = data2plot(:,asvs2plot);
    data2plot_sort = data2plot_filt(:,asvs2plot_sort);

    for i = 1:length(asvs)
        asv = asvs(i);
        fam = split(asv, ";g");
        color2plot = colors_kat{5}(strcmp(fam(1),OTUs_kat{5}),:); % finds the color code in kat's list of OTUs

        if isempty(color2plot)==0
            colors(i,:) = color2plot;

        elseif strcmp(asv, "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;NA;NA;NA")==1
            x = "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae";
            color2plot = colors_kat{5}(strcmp(x, OTUs_kat{5}),:); % finds the color code in kat's list of OTUs
            colors(i,:) = color2plot;

        else
            fam=split(asv, ";NA;NA");
            color2plot = colors_kat{5}(strcmp(fam(1),OTUs_kat{5}),:); % finds the color code in kat's list of OTUs
            colors(i,:) = color2plot;

        end
    end

    barplotconnect(data2plot_sort(:,:),1:length(idx2plot),colors,0,-1,0);

    %bar(data2plot_sort, 'stacked');
    %axis square
    xticklabels([])
    ylabel('Relative abundance')
    %set(gca,'FontSize',18);

    % save figure
    w=2.5;
    h=1.75;
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

    fname = strcat("../3-figs/fig1/niche_test_15mb_Eftitration_pass", string(p));
    print(gcf, strcat(fname, '.png'),'-dpng','-r300');
    print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');

end

%%