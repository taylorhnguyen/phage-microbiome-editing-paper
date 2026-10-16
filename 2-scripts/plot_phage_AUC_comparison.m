function [fig, p] = plot_phage_AUC_comparison_for230524(phage_name, pass, c, plate, color2plot, plottype)

% master variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows = ["A","B","C","D","E","F","G","H"];
cols = ["1","2","3","4","5","6","7","8","9","10","11","12"];
phage_names = ["WT","hk","tay4","tay62","tay119","tay120","tay121","tay122","tay126","tay129","cocktail"];

lag2plot = [];
maxOD2plot = [];
maxGR2plot = [];
AUC2plot = [];

ctrl_color = [0,0,0]./255; % black

figsdir = '/Users/taynguyen/Dropbox/_phage/analysis_forphagepaper/monoculture_MOI/figs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure

if plottype == "vertical"
    subplot(2,1,1)
elseif plottype == "horizontal"
    subplot(1,2,1)
end
          
for r = 1:8
        % plot controls
        well = strcat(rows(r),cols(1));
        plot(plate(pass).(well).time, plate(pass).(well).OD600blank, 'Linewidth', 1, 'Color', ctrl_color); hold all;
        %lag2plot(r,1) = plate(pass).(well).lag_time;
        %maxOD2plot(r,1) = max(plate(pass).(well).OD600blank);
        %maxGR2plot(r,1) = max(plate(pass).(well).growth_rate);
        AUC2plot(r,1) = trapz(plate(pass).(well).time, plate(pass).(well).OD600blank);
        
        % plot +phage
        well = strcat(rows(r),cols(c));
        plot(plate(pass).(well).time, plate(pass).(well).OD600blank, 'Linewidth', 1, 'Color', color2plot); hold all;
        %lag2plot(r,2) = plate(pass).(well).lag_time;
        %maxOD2plot(r,2) = max(plate(pass).(well).OD600blank);
        %maxGR2plot(r,2) = max(plate(pass).(well).growth_rate);
        AUC2plot(r,2) = trapz(plate(pass).(well).time, plate(pass).(well).OD600blank);
end

set(gca,'yscale','log')
xlim([0 48])
ylim([10^-3 1])
xlabel('time (h)');
ylabel('blanked OD600')
axis square


% add AUC
if plottype == "vertical"
    subplot(2,1,2)
elseif plottype == "horizontal"
    subplot(1,2,2)
end

AUC2plot = AUC2plot ./ (mean(AUC2plot(:,1)));
AUC_mean = mean(AUC2plot);
AUC_std = std(AUC2plot);
b = bar(AUC_mean, 'FaceColor', 'flat', 'LineWidth', 1.3, 'linestyle', 'none');
hold on

b.CData(1,:) = [255, 255, 255] ./ 255;
b.CData(2,:) = color2plot;

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
set(gca, 'XTickLabel', {'WT' phage_names(c)})
ylabel('normalized AUC');
ylim([0 1.3])
yticks([0 0.2 0.4 0.6 0.8 1.0])

% test for normal distribution
%[h,p] = kstest(AUC2plot(:,1)); % 1 / not normal
%[h,p] = kstest(AUC2plot(:,2)); % 1 / not normal

% test for equal variance
%[h,p] = vartest2(AUC2plot(:,1), AUC2plot(:,2)); % equal variance

%[h,p] = ttest2(AUC2plot(:,1), AUC2plot(:,2), "Tail", "right");

p = ranksum(AUC2plot(:,1), AUC2plot(:,2), "tail", "right")

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

xlim([0.5 2.5])

% save 
%w=2.2;
%h=4;
%set(gcf, 'Units', 'Inches', 'Position', [0, 0, w, h], 'PaperUnits', 'Inches', 'PaperSize', [w, h])

%fname = strcat(figsdir, '/', phage_name, '/ef_phage_', phage_name, '_pass', string(pass), '_MOI', MOI);
%print(gcf, strcat(fname, '.png'),'-dpng','-r300');
%print(gcf, strcat(fname, '.pdf'),'-dpdf','-r300');
%print(gcf, strcat(fname, '.svg'),'-dsvg','-r300');
%savefig(gcf, strcat(fname, '.fig'))
