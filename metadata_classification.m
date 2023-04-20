% Clinical information classification model
clear all; close all; clc
% Initialize random seed for reproducibility
rng('default')
addpath(genpath('utils'))

% Select folder with data
myDir = uigetdir(pwd,'Open database directory');

signalFiles = dir(fullfile(myDir,'*m.mat'));
headerFiles = dir(fullfile(myDir,'*.hea'));
LabelsVector = {};

count =1;
for k = 1:length(headerFiles)
    [X(count,:),y(count,:),names{count} ] = read_data(k,myDir);    
    % Exclude recordings done before 26 weeks of gestation
    if  X(count,1) < 26*7
        k;
        continue
    end
        
    count = count+1;
end

% Remove diagnoses of hypertension, diabetes, placental_position, and
% funneling, as explained in the paper.
X(:,[6,7,8,11])=[];

count = 1;
%% Repeat cross-validation 20 times with different sampling.
for j = 1:20
    % classification with crossval
    k=5;
    c = cvpartition(y(:,1),'KFold',k);
    y_test_cat = [];
    y_hat_cat = [];
    
    for i = 1:k
        idxTrain = training(c,i);
        idxTest = test(c,i);   
        XTrain = X(idxTrain,:);
        XTest = X(idxTest,:);
        X_train_no_nan = XTrain;
        X_test_no_nan  = XTest;
        
        % Complete missing entries
        for j = 1:size(XTrain,2)            
            mode_i = mode(XTrain(:,j));  
            is_nan_in_column = isnan(XTrain(:,j));
            X_train_no_nan(is_nan_in_column,j) = mode_i;
            is_nan_in_column_test = isnan(XTest(:,j));
            X_test_no_nan(is_nan_in_column_test,j) = mode_i;      
        end
        
        
        XTest = X_test_no_nan;        
        XTrain = X_train_no_nan;
        YTrain = y(idxTrain,1);
        mean_features = mean(XTrain);
        std_features = std(XTrain);
        
        try
            XTrain = (XTrain-mean_features)./std_features;
            XTest = (XTest-mean_features)./std_features;
            
            Mdl = fitclinear(XTrain,YTrain ,'learner','logistic','Regularization','lasso');
            
        catch
            % Cases where a feature may have the same value for all the
            % samples
            std_0 = find(std_features==0);
            XTrain = (XTrain-mean_features)./std_features;
            XTest = (XTest-mean_features)./std_features;
            for ii = 1:length(std_0)
                XTrain(:,std_0(ii)) = [];
                XTest(:,std_0(ii)) = [];
            end
            Mdl = fitclinear(XTrain,YTrain ,'learner','logistic','Regularization','lasso');
            
        end
                
        YTest = y(idxTest,1);
        y_hat = XTest*Mdl.Beta+ Mdl.Bias;
        y_hat_cell{count} = y_hat;
        y_test_cell{count} = YTest;
        
        count = count+1;
    end
end
%%

% Save data for analysis
% save('./results/roc_aucs/classification_metadata','y_test_cell','y_hat_cell');

