% Calcualte values in Table 2 of the manuscript
clear all; close all; clc
addpath('roc_aucs_5')
addpath(genpath('utils'))

% Load results from clinical information
results_class_meta = load('classification_metadata.mat');

%EHG
results_class_EHG = load('EHG_classification.mat');

%combined
results_class_combined = load('classification_combined.mat');

%% Ssensitivity


specificities = [50 70 90]/100;

% metadata
[X_auc_c_meta,Y_auc_c_meta,T_auc,AUC_total_meta_c] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1);
[X_auc_c_meta,Y_auc_c_meta,T_auc,AUC_total] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1,'XVals',X_auc_c_meta(:,1));


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_meta));
sensitivity_meta(i,:) = Y_auc_c_meta(idx(i),:);
end

% The columns indicate the mean, lower bound and upper bound of the 95% CI.
% The rows correspond to different specificity values
round(sensitivity_meta,3)*100

% EHG

[X_auc_c_EHG,Y_auc_c_EHG,T_auc,AUC_total_ehg_c] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1);
[X_auc_c_EHG,Y_auc_c_EHG,T_auc,AUC_total] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1,'XVals',X_auc_c_EHG(:,1));


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_EHG));
sensitivity_EHG(i,:) = Y_auc_c_EHG(idx(i),:);
end

round(sensitivity_EHG,3)*100


% Combined

[X_auc_c_comb,Y_auc_c_comb,T_auc,AUC_total_comb_c] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1);
[X_auc_c_comb,Y_auc_c_comb,T_auc,AUC_total] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1,'XVals',X_auc_c_comb(:,1));


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_comb));
sensitivity_comb(i,:) = Y_auc_c_comb(idx(i),:);
end

round(sensitivity_comb,3)*100

%%  PPV

% Clinical information
[X_auc_c_meta,Y_auc_c_meta_ppv,T_auc,AUC_total_meta_c] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1,'YCrit','ppv');
[X_auc_c_meta,Y_auc_c_meta_ppv,T_auc,AUC_total] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1,'XVals',X_auc_c_meta(:,1),'YCrit','ppv');


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_meta));
ppv_meta(i,:) = Y_auc_c_meta_ppv(idx(i),:);
end

round(ppv_meta,3)*100

% EHG

[X_auc_c_EHG,Y_auc_c_EHG_ppv,T_auc,AUC_total_EHG_c] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1,'YCrit','ppv');
[X_auc_c_EHG,Y_auc_c_EHG_ppv,T_auc,AUC_total] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1,'XVals',X_auc_c_EHG(:,1),'YCrit','ppv');


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_EHG));
ppv_EHG(i,:) = Y_auc_c_EHG_ppv(idx(i),:);
end

round(ppv_EHG,3)*100

% combined


[X_auc_c_combined,Y_auc_c_combined_ppv,T_auc,AUC_total_combined_c] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1,'YCrit','ppv');
[X_auc_c_combined,Y_auc_c_combined_ppv,T_auc,AUC_total] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1,'XVals',X_auc_c_combined(:,1),'YCrit','ppv');

for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_combined));
ppv_combined(i,:) = Y_auc_c_combined_ppv(idx(i),:);
end

round(ppv_combined,3)*100

%% NPV

% Clinical information

[X_auc_c_meta,Y_auc_c_meta_npv,T_auc,AUC_total_meta_c] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1,'YCrit','npv');
[X_auc_c_meta,Y_auc_c_meta_npv,T_auc,AUC_total] = perfcurve(results_class_meta.y_test_cell,results_class_meta.y_hat_cell,1,'XVals',X_auc_c_meta(:,1),'YCrit','npv');


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_meta));
npv_meta(i,:) = Y_auc_c_meta_npv(idx(i),:);
end

round(npv_meta,3)*100

% EHG

[X_auc_c_EHG,Y_auc_c_EHG_npv,T_auc,AUC_total_EHG_c] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1,'YCrit','npv');
[X_auc_c_EHG,Y_auc_c_EHG_npv,T_auc,AUC_total] = perfcurve(results_class_EHG.y_test_stored,results_class_EHG.y_hat_2,1,'XVals',X_auc_c_EHG(:,1),'YCrit','npv');


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_EHG));
npv_EHG(i,:) = Y_auc_c_EHG_npv(idx(i),:);
end

round(npv_EHG,3)*100


% combined

[X_auc_c_combined,Y_auc_c_combined_npv,T_auc,AUC_total_combined_c] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1,'YCrit','npv');
[X_auc_c_combined,Y_auc_c_combined_npv,T_auc,AUC_total] = perfcurve(results_class_combined.y_test_stored,results_class_combined.y_hat_stored,1,'XVals',X_auc_c_combined(:,1),'YCrit','npv');


for i=1:length(specificities)
[~,idx(i)] = min(abs((1-specificities(i))-X_auc_c_combined));
npv_combined(i,:) = Y_auc_c_combined_npv(idx(i),:);
end


round(npv_combined,3)*100

