% combined classification model

clear all; close all; clc
rng('default')
addpath(genpath('utils'))

myDir = uigetdir(pwd,'Open database directory');
addpath(myDir)

signalFiles = dir(fullfile(myDir,'*m.mat'));
headerFiles = dir(fullfile(myDir,'*.hea'));
LabelsVector = {};


%% load data


count =1;
for k = 1:length(headerFiles)
    
    % Read data and preprocess EHG
    [downsample_signal,t_downsample,is_preterm(count),name,gestational_age(count),delivery_age(count)] = preprocessFile(k,myDir);
    [X(count,:),y(count,:),names{count} ] = read_data(k,myDir);
    
    % Exclude recordings done before 26 weeks of gestation
    if gestational_age(count) < 26*7
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
    sizes(count) = size(s,2);
    % Zero-pad a 26-minute long recording
    if size(s,2)<106
        s = [s,zeros(size(s,1),106-size(s,2))];
    end
    % STFT is symmetric in frequency, so we take only the positive
    % frequencies.
    s_half{count} = s(round(size(s,1)/2):end,1:106);
    
    count = count+1;
    
    
end


% Remove diagnoses of hypertension, diabetes, placental_position,
% funneling, as explained in the paper.
X(:,[6,7,8,11])=[];


%% Repeat cross-validation 20 times with different sampling.
k=5;
count=1;
for jj=1:20
    c = cvpartition(y(:,1),'KFold',k);
    y_test_cat=[];
    y_hat_cat =[];
    for i =1:k
        
        idxTrain = training(c,i);
        idxTest = test(c,i);
        XTrain = X(idxTrain,:);
        XTest = X(idxTest,:);
        X_train_no_nan = XTrain;
        X_test_no_nan  = XTest;
        
        % Complete missing entries
        for j = 1:size(XTrain,2)
            mode_i = mean(XTrain(:,j),'omitnan'); 
            is_nan_in_column = isnan(XTrain(:,j));
            X_train_no_nan(is_nan_in_column,j) = mode_i;
            is_nan_in_column_test = isnan(XTest(:,j));
            X_test_no_nan(is_nan_in_column_test,j) = mode_i;
            
        end
        XTrain_metadata =X_train_no_nan;
        XTest_metadata = X_test_no_nan;
        YTrain_metadata = y(idxTrain,1);
        YTest_metadata = y(idxTest,1);
       
        XTrain_img = s_half(idxTrain);
        YTrain_img = categorical(is_preterm(idxTrain));       
        XTest_img = s_half(idxTest);
        YTest_img = categorical(is_preterm(idxTest));
        
        
        % Define network architecture
        numFeatures = size(s_half{1},1);
        numHiddenUnits = 100;
        numResponses = 1;
        numClasses = 2;
        
        classWeights = 1./countcats(YTrain_img);
        classWeights = classWeights'/mean(classWeights);
        
        layers = [ ...
            sequenceInputLayer(numFeatures)
            bilstmLayer(numHiddenUnits,'OutputMode','last','name','bilstm')
            fullyConnectedLayer(numClasses,'name','FC')
            softmaxLayer
            weightedClassificationLayer(classWeights')];
        
        maxEpochs = 30;
        miniBatchSize = 32; %12
        options = trainingOptions( 'adam',...
            'ExecutionEnvironment','cpu', ...
            'GradientThreshold',1, ...
            'MaxEpochs',maxEpochs, ...
            'ValidationData',{XTest_img,YTest_img'},...
            'ValidationFrequency', 8,...
            'InitialLearnRate',1e-4,... 
            'L2Regularization',1e-3,...
            'MiniBatchSize',miniBatchSize, ...
            'SequenceLength','longest', ...
            'Shuffle','every-epoch', ...
            'Verbose',1);
        

        
        net = trainNetwork(XTrain_img,YTrain_img',layers,options);
        
        % Get activation features to be combined with clinical information
        layer = 'FC';
        featuresTrain = activations(net,XTrain_img,layer);
        featuresTest = activations(net,XTest_img,layer);
        
        % Combine clinical info with features from EHG
        all_features_train = [XTrain_metadata, featuresTrain' ];
        mean_features = mean(all_features_train);
        std_features = std(all_features_train);
        all_features_train = (all_features_train-mean_features)./std_features;
        
        all_features_test = [XTest_metadata, featuresTest' ];
        all_features_test = (all_features_test-mean_features)./std_features;
        
        std_0 = find(std_features==0);
        for ii =1:length(std_0)
            all_features_train(:,std_0(ii)) = 0;
            all_features_test(:,std_0(ii)) = 0;
        end
        Mdl = fitclinear(all_features_train,YTrain_metadata ,'learner','logistic','Regularization','lasso');
        
        y_hat = all_features_test*Mdl.Beta+ Mdl.Bias;
        

        y_hat_stored{count} = y_hat;
        y_test_stored{count} = YTest_metadata;
        

        count=count+1;
    end
end


[X_auc_c_comb,Y_auc_c_comb,T_auc,AUC_total_comb_c] = perfcurve(y_test_stored,y_hat_stored,1);


%% Save data for analysis
%save('./results/roc_aucs/classification_combined','y_hat_stored','y_test_stored');
