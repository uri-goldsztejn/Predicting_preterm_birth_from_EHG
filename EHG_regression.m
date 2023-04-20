% EHG regression model

clear all; close all; clc
% Initialize random seed for reproducibility
rng('default')
addpath(genpath('utils'))

% Select folder with data
myDir = uigetdir(pwd,'Open database directory');
addpath(myDir)

signalFiles = dir(fullfile(myDir,'*m.mat'));
headerFiles = dir(fullfile(myDir,'*.hea'));
LabelsVector = {};

count =1;
for k = 1:length(headerFiles)
    % Preprocess data
    [downsample_signal,t_downsample,is_preterm(count),name,gestational_age(count),delivery_age(count)] = preprocessFile(k,myDir);
    % Exclude recordings done before 28 weeks of gestation
    if  gestational_age(count) < 26*7
        k;
        continue
    end
    
    % Truncate recordings longer than 30 minutes
    if t_downsample(end) > 30*60
        [~,idx] = min(abs(t_downsample-30*60));
        t_downsample = t_downsample(1:idx);
        downsample_signal = downsample_signal(1:idx,3);
    end
    
    
    names{count} = name;
    % Compute STFT
    fs = 1/(t_downsample(2)-t_downsample(1));
    f0 = 0.1;
    window_length = round(fs/f0*6);
    s = log(abs(stft(downsample_signal(:,1),fs,'Window',hamming(window_length),'Centered',0)));
    % Zero-pad a 26-minute long recording
    if size(s,2)<106
        s = [s,zeros(size(s,1),106-size(s,2))];
    end
    % STFT is symmetric in frequency, so we take only the positive
    % frequencies.
    s_half{count} = s(round(size(s,1)/2):end,1:106);
    
    count = count+1;
    
end


%% Repeat cross-validation 20 times with different sampling.
count=1;
for jj=1:20
    k=5;
    c = cvpartition(is_preterm,'KFold',k);
    y_test_cat=[];
    y_hat_cat =[];
    
    
    for i=1:k
        % Partition data
        idxTrain = training(c,i);
        XTrain = s_half(idxTrain);
        mean_train = mean(delivery_age(idxTrain));
        std_train = std(delivery_age(idxTrain));
        
        YTrain = (delivery_age(idxTrain) - mean_train)/std_train;
        idxTest = test(c,i);
        XTest = s_half(idxTest);
        YTest =(delivery_age(idxTest) - mean_train)/std_train;
        
        
        % Define network architecture
        numFeatures = size(s_half{1},1);
        numHiddenUnits = 100;
        numResponses = 1;
        numClasses = 2;
        
        
        layers = [ ...
            sequenceInputLayer(numFeatures)
            bilstmLayer(numHiddenUnits,'OutputMode','last')
            fullyConnectedLayer(2)
            fullyConnectedLayer(1)
            regressionLayer];
        
        maxEpochs = 25;
        miniBatchSize = 32;
        
        options = trainingOptions( 'adam',...
            'ExecutionEnvironment','cpu', ...
            'GradientThreshold',1, ...
            'MaxEpochs',maxEpochs, ...
            'ValidationData',{XTest,YTest'},...
            'ValidationFrequency', 8,...
            'InitialLearnRate',3e-3,...
            'L2Regularization',1e-3,...
            'MiniBatchSize',miniBatchSize, ...
            'SequenceLength','longest', ...
            'Shuffle','every-epoch', ...
            'Verbose',1);
        
        net = trainNetwork(XTrain,YTrain',layers,options);
        
        y_hat{count}   = (predict(net,XTest)*std_train)+ mean_train;
        YTest = (YTest*std_train)+ mean_train;
        y_test_stored{count} = YTest';
        
        count = count+1;
    end
end

%% Save data for analysis
for i=1:length(y_test_stored)

y_test_stored_binary{i} =y_test_stored{i}>258; 
end

[X_auc_r_EHG,Y_auc_r_EHG,T_auc,AUC_total_ehg_r] = perfcurve(y_test_stored_binary,y_hat,1);

% save('./results/roc_aucs/EHG_regression_69','y_hat','y_test_stored');
