function [p_vals] = plot_phage_growth_comparison(c, pass, plate)

% master variables
p_vals = [];
rows = ["A","B","C","D","E","F","G","H"];
cols = ["1","2","3","4","5","6","7","8","9","10","11","12"];

phage_names = ["WT","hk","tay4","tay62","tay119","tay120","tay121","tay122","tay126","tay129","cocktail"];

figure
subplot(1,5,1)

lag2plot = [];
maxOD2plot = [];
maxGR2plot = [];
AUC2plot = [];
colors2plot = [0,0,0; % WT grey
                189,189,189; % hk
                201,148,199; % tay 4 hnl purple
                255,127,0; % tay62 richmond dark green
                253,191,111; % tay119 richmond light green
                251,128,114; %tay120 codiga pink
                51,160,44; %tay121 richmond day 2 orange
                178,223,138; % tay122 richmond day 2 light orange
                33,113,181; % tay126 dark blue
                107,174,214; % tay129 light blue
                165,15,21] ./ 255; % cocktail red
              
% make dataframes & plot OD600 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:8
        % plot controls
        well = strcat(rows(r),cols(1));
        plot(plate(pass).(well).time, plate(pass).(well).OD600blank, 'Linewidth', 1, 'Color', colors2plot(1,:)); hold all;
        lag2plot(r,1) = plate(pass).(well).lag_time;
        maxOD2plot(r,1) = max(plate(pass).(well).OD600blank);
        maxGR2plot(r,1) = plate(pass).(well).growth_rate_max;
        AUC2plot(r,1) = trapz(plate(pass).(well).time, plate(pass).(well).OD600blank);
        
        % plot +phage
        well = strcat(rows(r),cols(c));
        plot(plate(pass).(well).time, plate(pass).(well).OD600blank, 'Linewidth', 1, 'Color', colors2plot(c,:)); hold all;
        lag2plot(r,2) = plate(pass).(well).lag_time;
        maxOD2plot(r,2) = max(plate(pass).(well).OD600blank);
        maxGR2plot(r,2) = plate(pass).(well).growth_rate_max;
        AUC2plot(r,2) = trapz(plate(pass).(well).time, plate(pass).(well).OD600blank);
end

% other figure things
set(gca,'yscale','log')
xlim([0 48])
ylim([10^-3 1])
xlabel('time (h)');
ylabel('blanked OD600')
axis square

% plot lag time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,5,2)
lag_mean = mean(lag2plot);
lag_std = std(lag2plot);
b = bar(lag_mean, 'EdgeColor', 'none'); hold all;
b.FaceColor = 'flat';

b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = colors2plot(c,:);

% draw rectangle for WT because can't change edges for bars
% Get x position and height
x = b.XData(1);
h = lag_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, ...
    'FaceColor', 'none');

% add errorbars
er = errorbar(lag_mean,lag_std);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% add scatter plot
for k = 1:2
    scatter(repelem(k,8),lag2plot(:,k), [],'black', 'HandleVisibility','off'); hold all;
end

% add other figure things
axis square
set(gca, 'XTickLabel', {'WT' phage_names(c)})
ylabel('lag time (h)');
ylim([0 max(lag2plot,[],'all')*1.2])
xlim([0.5 2.5])

% add p-values
[p,h] = ranksum(lag2plot(:,1), lag2plot(:,2));
p_vals(1) = p;
x = 1.4; y = max(lag2plot,[],'all')*1.1;
if p > .05 | isnan(p)==1
    text(x, y, 'NS');
elseif p < 0.05 & p > 0.005
    text(x, y, '*');
elseif p < 0.005 & p > 0.0005
    text(x, y, '**');
elseif p < 0.0005 & p > 0.00005
    text(x, y, '***');
else
    text(x, y, '****');
end

% plot max OD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,5,3)
maxOD_mean = mean(maxOD2plot);
maxOD_std = std(maxOD2plot);
b = bar(maxOD_mean, 'EdgeColor', 'none'); hold all;
b.FaceColor = 'flat';

b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = colors2plot(c,:);

% draw rectangle for WT because can't change edges for bars
% Get x position and height
x = b.XData(1);
h = maxOD_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, ...
    'FaceColor', 'none');

% add error bar
er = errorbar(maxOD_mean, maxOD_std);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% add scatter plot
for k = 1:2
    scatter(repelem(k,8),maxOD2plot(:,k), [],'black', 'HandleVisibility','off'); hold all;
end

% other figure things
axis square
set(gca, 'XTickLabel', {'WT' phage_names(c)})
ylabel('max yield (OD600)');
ylim([0 max(maxOD2plot,[],'all')*1.2])
xlim([0.5 2.5])

% add p-values
[p,h] = ranksum(maxOD2plot(:,1), maxOD2plot(:,2));
p_vals(2) = p;
x = 1.4; y = max(maxOD2plot,[],'all')*1.1;
if p > .05 | isnan(p)==1
    text(x, y, 'NS');
elseif p < 0.05 & p > 0.005
    text(x, y, '*');
elseif p < 0.005 & p > 0.0005
    text(x, y, '**');
elseif p < 0.0005 & p > 0.00005
    text(x, y, '***');
else
    text(x, y, '****');
end

% plot max growth rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,5,4)
maxGR_mean = mean(maxGR2plot);
maxGR_std = std(maxGR2plot);
b = bar(maxGR_mean, 'EdgeColor', 'none'); hold all;
b.FaceColor = 'flat';

b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = colors2plot(c,:);

% draw rectangle for WT because can't change edges for bars
% Get x position and height
x = b.XData(1);
h = maxGR_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, ...
    'FaceColor', 'none');

% add errorbar
er = errorbar(maxGR_mean, maxGR_std);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% add scatter plot
for k = 1:2
    scatter(repelem(k,8),maxGR2plot(:,k), [],'black', 'HandleVisibility','off'); hold all;
end

% other figure things
axis square
set(gca, 'XTickLabel', {'WT' phage_names(c)})
ylabel('max growth rate (hr^{-1})');
ylim([0 max(maxGR2plot,[],'all')*1.2])
xlim([0.5 2.5])

% add p-values
[p,h] = ranksum(maxGR2plot(:,1), maxGR2plot(:,2));
p_vals(3) = p;
x = 1.4; y = max(maxGR2plot,[],'all')*1.1;
if p > .05 | isnan(p)==1
    text(x, y, 'NS');
elseif p < 0.05 & p > 0.005
    text(x, y, '*');
elseif p < 0.005 & p > 0.0005
    text(x, y, '**');
elseif p < 0.0005 & p > 0.00005
    text(x, y, '***');
else
    text(x, y, '****');
end

% plot AUC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,5,5)
AUC2plot = AUC2plot ./ (mean(AUC2plot(:,1)));
AUC_mean = mean(AUC2plot);
AUC_std = std(AUC2plot);
b = bar(AUC_mean, 'EdgeColor', 'none'); hold all;
b.FaceColor = 'flat';
b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = colors2plot(c,:);

% draw rectangle for WT because can't change edges for bars
% Get x position and height
x = b.XData(1);
h = AUC_mean(1);
rectangle('Position', [x - 0.4, 0, 0.8, h], ...
    'EdgeColor', '0,0,0', ...  % set to black
    'LineWidth', 1, ...
    'FaceColor', 'none');

% add errorbars
er = errorbar(AUC_mean, AUC_std);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

% add scatter plot
for k = 1:2
    scatter(repelem(k,8),AUC2plot(:,k), [],'black', 'HandleVisibility','off'); hold all;
end

% other figure things
axis square
set(gca, 'XTickLabel', {'WT' phage_names(c)})
ylabel('normalized AUC');
ylim([0 max(AUC2plot,[],'all')*1.2])
yticks([0 0.2 0.4 0.6 0.8 1.0])
xlim([0.5 2.5])

% add p-values
[p,h] = ranksum(AUC2plot(:,1), AUC2plot(:,2), "Tail", "right");
%[p,h] = ranksum(AUC2plot(:,1), AUC2plot(:,2));
p_vals(4) = p;
x = 1.4; y = max(AUC2plot,[],'all')*1.1;
if p > .05 | isnan(p)==1
    text(x, y, 'NS');
elseif p < 0.05 & p > 0.005
    text(x, y, '*');
elseif p < 0.005 & p > 0.0005
    text(x, y, '**');
elseif p < 0.0005 & p > 0.00005
    text(x, y, '***');
else
    text(x, y, '****');
end

% save figure in main script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

