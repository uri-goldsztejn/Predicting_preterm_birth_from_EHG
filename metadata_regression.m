% Clinical information regression model
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
    % Exclude recordings done before 28 weeks of gestation
    if X(count,1) < 26*7
        k;
        continue
    end 
    count = count+1;  
end

% Remove diagnoses of hypertension, diabetes, placental_position,
% funneling, as explained in the paper.
X(:,[6,7,8,11])=[];


k=5;
count = 1;
%% Repeat cross-validation 20 times with different sampling.
for jj=1:20
    % Regression with cross-validation
    c = cvpartition(y(:,1),'KFold',k);
    y_test_cat=[];
    y_hat_cat =[];
    
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
        YTrain = y(idxTrain,2);
        Mdl = fitlm(XTrain,YTrain);
        YTest = y(idxTest,2);
        y_hat = predict(Mdl,XTest);
        y_hat_cell{count} = y_hat;
        y_test_cell{count} = YTest;
        

        count = count+1;
    end
    
end

%%


% Save data for analysis
save('./results/roc_aucs/regression_metadata','y_test_cell','y_hat_cell','y');
