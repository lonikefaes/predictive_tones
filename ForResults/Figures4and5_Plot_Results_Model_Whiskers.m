clear; clc
load('output_t2_model_Mispred_ptcutoffpt01_notflipped_depths791011.mat')
outfolder = ['\PATHtoData\'];

cmap = colormap(jet);cmap = cmap(25:52:end,:);cmap = [cmap;zeros(2,3);cmap;zeros(2,3);cmap];

%%

% ======================== Mispred versus Pred
figure(1002); clf;

subplot(1,5,1);
boxplot(resp_HG_oddball,'Colors',cmap(5,:),'Whisker',1, 'Widths',[1,1],'Symbol', '','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-3 7]; x.XTick = [-2 0 2 4 6]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on;
xline(0);title('HG');   
ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted

subplot(1,5,2);
boxplot(resp_PP_oddball,'Colors',cmap(4,:),'Whisker',1, 'Widths',1,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-3 7]; x.XTick = [-2 0 2 4 6]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('PP');

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted

subplot(1,5,3);
boxplot(resp_PT_oddball,'Colors',cmap(3,:),'Whisker',1, 'Widths',1,'Symbol','-','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-3 7]; x.XTick = [-2 0 2 4 6]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('PT');

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted

subplot(1,5,4);
boxplot(resp_aSTG_oddball,'Colors',cmap(2,:),'Whisker',1, 'Widths',1,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-3 7]; x.XTick = [-2 0 2 4 6]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('aSTG');

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted


subplot(1,5,5);
boxplot(resp_pSTG_oddball,'Colors',cmap(1,:),'Whisker',1, 'Widths',1,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-3 7]; x.XTick = [-2 0 2 4 6]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('pSTG');

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted

c = gcf; 

set(findobj(gcf,'type','Axes'),'FontName','Arial','FontSize',10);
sgtitle('Mispredicted versus Predictable', 'FontName', 'Arial', 'FontSize', 16, 'fontweight', 'bold')



% saveas(c, [outfolder, 'Results_MispredvsPred_Model_Whiskers.svg'])
% saveas(c, [outfolder, 'Results_MispredvsPred_Model_Whiskers.png'])

%% Plot the Unpredictable versus Predictable conditions
% ======================== unexp
figure(1003)

load('output_t2_model_Unpred_ptcutoffpt01_notflipped_depths791011.mat')

subplot(1,5,1);
boxplot(resp_HG_unexp,'Colors',cmap(5,:),'Whisker',1, 'Widths',0.3,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-2 4]; x.XTick = [-2 0 2 4]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on;
xline(0);title('HG', 'fontname', 'Arial', 'fontsize', 12);   

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted


subplot(1,5,2);
boxplot(resp_PP_unexp,'Colors',cmap(4,:),'Whisker',1, 'Widths',0.3,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-2 4]; x.XTick = [-2 0 2 4]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('PP', 'fontname', 'Arial', 'fontsize', 12);

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted


subplot(1,5,3);
boxplot(resp_PT_unexp,'Colors',cmap(3,:),'Whisker',1, 'Widths',0.3,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-2 4]; x.XTick = [-2 0 2 4]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('PT', 'fontname', 'Arial', 'fontsize', 12);

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted


subplot(1,5,4);
boxplot(resp_aSTG_unexp,'Colors',cmap(2,:),'Whisker',1, 'Widths',0.3,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-2 4]; x.XTick = [-2 0 2 4]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('aSTG', 'fontname', 'Arial', 'fontsize', 12);

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted


subplot(1,5,5);
boxplot(resp_pSTG_unexp,'Colors',cmap(1,:),'Whisker',1, 'Widths',0.3,'Symbol','','orientation','horizontal','PlotStyle','compact'); hold on
x = gca; x.YLim = [0 4];x.XLim = [-2 4]; x.XTick = [-2 0 2 4]; x.Box = 'off'; x.YAxis.Visible = 'off'; x.XAxisLocation = 'top';grid on
xline(0);title('pSTG', 'fontname', 'Arial', 'fontsize', 12);

ax = gca;
% do all the plotting
ax.Position(2) = ax.Position(2)*0.05; % This may be adjusted

c = gcf; 

set(findobj(gcf,'type','Axes'),'FontName','Arial','FontSize',10);
sgtitle('Unpredictable versus Predictable', 'FontName', 'Arial', 'FontSize', 16, 'fontweight', 'bold')


b = gca;
b.XAxis.FontSize = 12;
b.XAxis.FontName = 'Arial';
b.YAxis.FontSize = 12;
b.YAxis.FontName = 'Arial';

saveas(c, [outfolder, 'Results_UnpredvsPred_Model_Whiskers.svg'])
saveas(c, [outfolder, 'Results_UnpredvsPred_Model.png'])
