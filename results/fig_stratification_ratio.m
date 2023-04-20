% Create figure illustrating data partitioning in the cross-validation
clear all; close all; clc
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
    [downsample_signal,t_downsample,is_preterm(count),name,gestational_age(count),delivery_age(count)] = preprocessFile(k,myDir);
    % Exclude recordings done before 26 weeks of gestation
    if gestational_age(count) < 26*7
        k;
        continue
    end
    
    count = count+1;
end
%%

% Create a five-fold cross-validation as in the other scripts.
k=5;
c = cvpartition(is_preterm,'KFold',k);
y_test_cat=[];
y_hat_cat =[];


for i=1:k
    idxTrain = training(c,i);
end

numFolds = 5;
numClasses=2;
nTestData = zeros(numFolds,numClasses);
for i = 1:numFolds
    testClasses = is_preterm(c.test(i));
    nCounts = countcats(categorical(testClasses)); % Number of test set observations in each class
    nTestData(i,:) = nCounts';
end

nTrainData = zeros(numFolds,numClasses);
for i = 1:numFolds
    trainClasses = is_preterm(c.training(i));
    nCounts = countcats(categorical(trainClasses)); % Number of test set observations in each class
    nTrainData(i,:) = nCounts';
end

figure('defaultAxesFontName','Arial','DefaultAxesFontSize',14)
% Show partitioning in the training set
subplot(1,2,1)
bar(nTrainData)
set(gca,'ytick', [0,25,50,75,100])
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'Training set (fold)', '(a)'},'FontSize',14)
ylabel('Number of observations','FontSize',14)
ylim([0 110])

% Show partitioning in the test set
subplot(1,2,2)
bar(nTestData)
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gca,'ytick', [0,10,20])
ylim([0 28])
xlabel({'Test set (fold)', '(b)'},'FontSize',14)
legend('Term','Preterm','FontSize',12)

