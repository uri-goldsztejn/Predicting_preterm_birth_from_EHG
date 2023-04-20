% This script generates the ROC curves shown in Fig. 2
clear all;  clc; close all;
rng('default')
addpath('results/roc_aucs')
addpath(genpath('utils'))

% Load clinical information data
results_class_meta = load('classification_metadata.mat');
results_reg_meta = load('regression_metadata.mat');

% EHG
results_class_EHG = load('EHG_classification.mat');
results_reg_EHG = load('EHG_regression.mat');

% combined
results_class_combined = load('classification_combined.mat');
results_reg_combined = load('regression_combined.mat');


%% Performance bound

for i = 1:20
    y_days{i} = results_reg_meta.y(:,2);
    y_days_class{i} = y_days{i} > 258;
    y_with_noise{i} = y_days{i} + 7*randn(size(y_days{i}));
end
[X_ultrasound,Y_ultrasound,T_auc,AUC_ultrasound] = perfcurve(y_days_class,y_with_noise,1);
[X_ultrasound,Y_ultrasound,T_auc,AUC_ultrasound]  = perfcurve(y_days_class,y_with_noise,1,'XVals',X_ultrasound(:,1));

% References for ultrasound error
% https://onlinelibrary.wiley.com/doi/full/10.1016/j.jmwh.2008.11.003
% https://journals.lww.com/greenjournal/Fulltext/2017/05000/Committee_Opinion_No_700__Methods_for_Estimating.50.aspx


%%
linewidth = 2;
figure('defaultAxesFontName','Arial','DefaultAxesFontSize',14)

subplot(1,2,1)

legend_fontsize = 10;


% plot bound
plot(X_ultrasound(:,1),Y_ultrasound(:,1),'LineWidth',linewidth,'Color','k');hold on
xconf = [X_ultrasound(:,1)' flip(X_ultrasound(:,1))'] ;
yconf = [Y_ultrasound(:,2)' flip(Y_ultrasound(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','k','EdgeColor','none', 'FaceAlpha',0.6);


xconf = [X_ultrasound(:,1)'  0] ;
yconf = [(Y_ultrasound(:,3))' 1] ;
p = fill(xconf,yconf,'red','FaceColor','k','EdgeColor','none', 'FaceAlpha',0.2);


% combined

[X_auc_c_comb,Y_auc_c_comb,T_auc,AUC_total_comb_c] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1);
[X_auc_c_comb,Y_auc_c_comb,T_auc,AUC_total] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1,'XVals',X_auc_c_comb(:,1));
plot(X_auc_c_comb(:,1),Y_auc_c_comb(:,1),'LineWidth',linewidth,'Color','#0072BD');hold on
xconf = [X_auc_c_comb(:,1)' flip(X_auc_c_comb(:,1))'] ;
yconf = [Y_auc_c_comb(:,2)' flip(Y_auc_c_comb(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','#0072BD','EdgeColor','none', 'FaceAlpha',0.25);



% EHG
[X_auc_c_EHG,Y_auc_c_EHG,T_auc,AUC_total_ehg_c] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1);
[X_auc_c_EHG,Y_auc_c_EHG,T_auc,AUC_total] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1,'XVals',X_auc_c_EHG(:,1));
plot(X_auc_c_EHG(:,1),Y_auc_c_EHG(:,1),'LineWidth',linewidth,'Color','#D95319');hold on
xconf = [X_auc_c_EHG(:,1)' flip(X_auc_c_EHG(:,1))'] ;
yconf = [Y_auc_c_EHG(:,2)' flip(Y_auc_c_EHG(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','#D95319','EdgeColor','none', 'FaceAlpha',0.25);


%metadata

[X_auc_c_meta,Y_auc_c_meta,T_auc,AUC_total_meta_c] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1);
[X_auc_c_meta,Y_auc_c_meta,T_auc,AUC_total] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1,'XVals',X_auc_c_meta(:,1));
plot(X_auc_c_meta(:,1),Y_auc_c_meta(:,1),'LineWidth',linewidth,'Color','#EDB120');hold on
xconf = [X_auc_c_meta(:,1)' flip(X_auc_c_meta(:,1))'] ;
yconf = [Y_auc_c_meta(:,2)' flip(Y_auc_c_meta(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','#EDB120','EdgeColor','none', 'FaceAlpha',0.25);



title('Classification models')

set(gca,'xtick', [0:0.2:1])
set(gca,'xticklabel',  [0:0.2:1])

set(gca,'ytick', [0:0.2:1])
set(gca,'yticklabel', [0:0.2:1])
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

xlabel({'False positive rate','(1 – specificity)','(a)'},'FontSize',14 );
ylabel('True positive rate','FontSize',15 )

xlim([0 1])
ylim([0 1])

% Values are obtained from the variables "AUC_total"
lh = legend('','','',...
    [' Combined' newline ' AUC = 0.78' newline ' (0.76 - 0.80)'],'',...
    [' EHG' newline ' AUC = 0.74' newline ' (0.73 - 0.76)'],'',...
    [ ' Clinical info.' newline ' AUC = 0.65' newline ' (0.63 - 0.67)'],'','Location','southeast','FontSize',legend_fontsize );
legend('boxoff')



%% Regression

subplot(1,2,2)

% bound
plot(X_ultrasound(:,1),Y_ultrasound(:,1),'LineWidth',linewidth,'Color','k');hold on
xconf = [X_ultrasound(:,1)' flip(X_ultrasound(:,1))'] ;
yconf = [Y_ultrasound(:,2)' flip(Y_ultrasound(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','k','EdgeColor','none', 'FaceAlpha',0.6);


xconf = [X_ultrasound(:,1)'  0] ;
yconf = [(Y_ultrasound(:,3))' 1] ;
p = fill(xconf,yconf,'red','FaceColor','k','EdgeColor','none', 'FaceAlpha',0.2);


for i=1:length(results_reg_combined.y_test_stored)
    
    results_reg_combined.y_test_stored_binary{i} =results_reg_combined.y_test_stored{i}>258;
end

[X_auc_r_comb,Y_auc_r_comb,T_auc,AUC_total_comb_r] = perfcurve(results_reg_combined.y_test_stored_binary,results_reg_combined.y_hat_stored,1);
[X_auc_r_comb,Y_auc_r_comb,T_auc,AUC_total] = perfcurve(results_reg_combined.y_test_stored_binary,results_reg_combined.y_hat_stored,1,'XVals',X_auc_r_comb(:,1));
plot(X_auc_r_comb(:,1),Y_auc_r_comb(:,1),'LineWidth',linewidth,'Color','#0072BD');hold on
xconf = [X_auc_r_comb(:,1)' flip(X_auc_r_comb(:,1))'] ;
yconf = [Y_auc_r_comb(:,2)' flip(Y_auc_r_comb(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','#0072BD','EdgeColor','none', 'FaceAlpha',0.25);



for i=1:length(results_reg_EHG.y_test_stored)
    
    results_reg_EHG.y_test_stored_binary{i} = results_reg_EHG.y_test_stored{i}>258;
end

[X_auc_r_EHG,Y_auc_r_EHG,T_auc,AUC_total_ehg_r] = perfcurve(results_reg_EHG.y_test_stored_binary,results_reg_EHG.y_hat,1);
[X_auc_r_EHG,Y_auc_r_EHG,T_auc,AUC_total] = perfcurve(results_reg_EHG.y_test_stored_binary,results_reg_EHG.y_hat,1,'XVals',X_auc_r_EHG(:,1));
plot(X_auc_r_EHG(:,1),Y_auc_r_EHG(:,1),'LineWidth',linewidth,'Color','#D95319');hold on
xconf = [X_auc_r_EHG(:,1)' flip(X_auc_r_EHG(:,1))'] ;
yconf = [Y_auc_r_EHG(:,2)' flip(Y_auc_r_EHG(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','#D95319','EdgeColor','none', 'FaceAlpha',0.25);



for i=1:length(results_reg_meta.y_test_cell)
    
    results_reg_meta.y_test_cell_binary{i} = results_reg_meta.y_test_cell{i}>258;
end

[X_auc_r_meta,Y_auc_r_meta,T_auc,AUC_total_meta_r] = perfcurve(results_reg_meta.y_test_cell_binary,results_reg_meta.y_hat_cell,1);
[X_auc_r_meta,Y_auc_r_meta,T_auc,AUC_total] = perfcurve(results_reg_meta.y_test_cell_binary,results_reg_meta.y_hat_cell,1,'XVals',X_auc_r_meta(:,1));
plot(X_auc_r_meta(:,1),Y_auc_r_meta(:,1),'LineWidth',linewidth,'Color','#EDB120');hold on
xconf = [X_auc_r_meta(:,1)' flip(X_auc_r_meta(:,1))'] ;
yconf = [Y_auc_r_meta(:,2)' flip(Y_auc_r_meta(:,3))'] ;
p = fill(xconf,yconf,'red','FaceColor','#EDB120','EdgeColor','none', 'FaceAlpha',0.25);


% Figure axes
title('Regression models')
set(gca,'xtick', [0:0.2:1])
set(gca,'xticklabel',  [0:0.2:1])

set(gca,'ytick', [0:0.2:1])
set(gca,'yticklabel', [0:0.2:1])
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'False positive rate','(1 – specificity)', '(b)'},'FontSize',14 );



xlim([0 1])
ylim([0 1])
% Values are obtained from the variables "AUC_total"

legend([' Bound' newline ' AUC = 0.98' newline ' (0.98 - 0.98)'],'','',...
    [' Combined' newline ' AUC = 0.75' newline ' (0.73 - 0.77)'],'',...
    [' EHG' newline ' AUC = 0.70' newline ' (0.68 - 0.73)'],'',...
    [ ' Clinical info.' newline ' AUC = 0.67' newline ' (0.65 - 0.70)'],'','Location','southeast','FontSize',legend_fontsize )
legend('boxoff')

set(gcf,'color','w');

