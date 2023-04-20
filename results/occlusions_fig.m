% figure 3 in manuscript - occlusions

clear all; close all; clc
rng('default')
addpath('occlusions/results')
addpath('results/roc_aucs')
addpath(genpath('utils'))
addpath(genpath('metadata_model'));

% Select folder with database
myDir = uigetdir(pwd,'Open database directory');

signalFiles = dir(fullfile(myDir,'*m.mat'));
headerFiles = dir(fullfile(myDir,'*.hea'));
LabelsVector = {};


count =1;
k = 1;

% Preprocess sample
[downsample_signal,t_downsample,is_preterm(count),name,gestational_age(count),delivery_age(count)] = preprocessFile(k,myDir);

names{count} = name;
fs = 1/(t_downsample(2)-t_downsample(1));

f0 = 0.1;
window_length = round(fs/f0*6);
% Compute STFT
s = log(abs(stft(downsample_signal(:,1),fs,'Window',hamming(window_length),'Centered',0)));
s_half = s(round(size(s,1)/2):end,1:106);
STFT_img = ind2rgb(uint8(rescale(s_half,1,128)),jet(128));


%% frequency bands

% Open figure
figure('defaultAxesFontName','Arial','DefaultAxesFontSize',13)
t = tiledlayout(3,3,'TileSpacing','loose');
% Plot entire STFT
ax = nexttile([1 2]);
stft(downsample_signal(:,1),fs,'Window',hamming(window_length),'Centered',0);
colorbar off

xl = xlim;
xticks = linspace(xl(1),xl(2),7);
set(gca,'xtick', xticks(1:end-1))
set(gca,'xticklabel', [])

set(gca,'ytick', [0:2:4])

xlabel('(a)')
ylabel('Frequency [Hz]','FontSize',12)

yline(1,'--','Linewidth',3,'Color','w')
yline(2.2,'--','Linewidth',3,'Color','w')
yline(3.5,'--','Linewidth',3,'Color','w')

fonts_size_bracket = 18;
fonts_size_letter = 11;

% Mark frequency bands
annotation('textbox', [0.605, 0.72, 0.05, 0.06], 'String', "}",'FontName','Arial','FontSize',fonts_size_bracket,'LineStyle','none')
annotation('textbox', [0.6+0.02, 0.72-0.01, 0.05, 0.06], 'String', "B0",'FontName','Arial','FontSize',fonts_size_letter,'LineStyle','none')
annotation('textbox', [0.605, 0.765, 0.05, 0.06], 'String', "}",'FontName','Arial','FontSize',fonts_size_bracket,'LineStyle','none')
annotation('textbox', [0.6+0.02, 0.765-0.01, 0.05, 0.06], 'String', "B1",'FontName','Arial','FontSize',fonts_size_letter,'LineStyle','none')
annotation('textbox', [0.605, 0.81, 0.05, 0.06], 'String', "}",'FontName','Arial','FontSize',fonts_size_bracket,'LineStyle','none')
annotation('textbox', [0.6+0.02, 0.81-0.01, 0.05, 0.06], 'String', "B2",'FontName','Arial','FontSize',fonts_size_letter,'LineStyle','none')
annotation('textbox', [0.605, 0.86, 0.05, 0.06], 'String', "}",'FontName','Arial','FontSize',fonts_size_bracket,'LineStyle','none')
annotation('textbox', [0.6+0.02, 0.86-0.01, 0.05, 0.06], 'String', "B3",'FontName','Arial','FontSize',fonts_size_letter,'LineStyle','none')

title("")
ylim([0 5])

nexttile

results_EHG_B0 = load('results_EHG_B0.mat');
results_EHG_B1 = load('results_EHG_B1.mat');
results_EHG_B2 = load('results_EHG_B2.mat');
results_EHG_B3 = load('results_EHG_B3.mat');

for i=1:length(results_EHG_B0.y_hat)
    results_EHG_B0.y_hat_2{i} = results_EHG_B0.y_hat{i}(:,2);
    results_EHG_B1.y_hat_2{i} = results_EHG_B1.y_hat{i}(:,2);
    results_EHG_B2.y_hat_2{i} = results_EHG_B2.y_hat{i}(:,2);
    results_EHG_B3.y_hat_2{i} = results_EHG_B3.y_hat{i}(:,2);
    
end


% Add results of EHG bands
[~,~,~,AUC_B0] = perfcurve(results_EHG_B0.y_test_stored,results_EHG_B0.y_hat_2,1);
[~,~,~,AUC_B1] = perfcurve(results_EHG_B1.y_test_stored,results_EHG_B1.y_hat_2,1);
[~,~,~,AUC_B2] = perfcurve(results_EHG_B2.y_test_stored,results_EHG_B2.y_hat_2,1);
[~,~,~,AUC_B3] = perfcurve(results_EHG_B3.y_test_stored,results_EHG_B3.y_hat_2,1);


% B0 B1 B2 B3
AUCs = [AUC_B0(1), AUC_B1(1), AUC_B2(1), AUC_B3(1)];


CI_neg = [(AUC_B0(1)-AUC_B0(2)), (AUC_B1(1)-AUC_B1(2)), (AUC_B2(1)-AUC_B2(2)), (AUC_B3(1)-AUC_B3(2)) ];
CI_pos = [(AUC_B0(3)-AUC_B0(1)), (AUC_B1(3)-AUC_B1(1)), (AUC_B2(3)-AUC_B2(1)), (AUC_B3(3)-AUC_B3(1)) ];

linewidth = 2;
capsize = 10;
marker_size = 3;
errorbar(1:4,AUCs,CI_neg,CI_neg,':ok','MarkerSize',marker_size,...
    'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',linewidth,'CapSize',capsize);hold on

plot(1:4,AUCs,':','LineWidth',linewidth,'Color','#0072BD');

xlabel({'Frequency bands', '(b)'});ylabel('AUC')
set(gca,'xtick', 1:4)
set(gca,'xticklabel',  {'B0','B1','B2', 'B3'});
set(gca,'ytick', [0.5:0.2:0.9])
xlim([0.5 4.5])
ylim([0.5 1])




ax2 = nexttile([1 2]);
stft(downsample_signal(:,1),fs,'Window',hamming(window_length),'Centered',0); hold on
colorbar off

xl = xlim;
xticks = linspace(xl(1),xl(2),7);
set(gca,'xtick', xticks(1:end-1))
set(gca,'xticklabel', [])
set(gca,'ytick', [0:2:4])

xlabel('(c)')
ylabel('Frequency [Hz]','FontSize',12)

ylim([0 5])


image_data = getimage(ax2);


child_image = get(ax2,'Children');
child_image.CData = child_image.CData(:,randperm(size(child_image.CData,2)));


rectangle('Position',[7,0,1,5],'LineWidth',2)
rectangle('Position',[15,0,1,5],'LineWidth',2)


% Random temporal rearrangement of STFT
title("")

axd = nexttile;
axd.Clipping = 'off';
results_EHG_rand_0 = load('EHG_classification.mat');
results_EHG_rand_25 = load('EHG_results_rand_25.mat');
results_EHG_rand_50 = load('EHG_results_rand_50.mat');
results_EHG_rand_75 = load('EHG_results_rand_75.mat');
results_EHG_rand_100 = load('EHG_results_rand_100.mat');


for i=1:length(results_EHG_rand_0.y_hat_2)
    results_EHG_rand_0.y_hat_2{i} = results_EHG_rand_0.y_hat_2{i};%(:,2);
    results_EHG_rand_25.y_hat_2{i} = results_EHG_rand_25.y_hat{i}(:,2);
    results_EHG_rand_50.y_hat_2{i} = results_EHG_rand_50.y_hat{i}(:,2);
    results_EHG_rand_75.y_hat_2{i} = results_EHG_rand_75.y_hat{i}(:,2);
    results_EHG_rand_100.y_hat_2{i} = results_EHG_rand_100.y_hat{i}(:,2);
    
end

% EHG
[~,~,~,AUC_rand_0] = perfcurve(results_EHG_rand_0.y_test_stored,results_EHG_rand_0.y_hat_2,1);
[~,~,~,AUC_rand_25] = perfcurve(results_EHG_rand_25.y_test_stored,results_EHG_rand_25.y_hat_2,1);
[~,~,~,AUC_rand_50] = perfcurve(results_EHG_rand_50.y_test_stored,results_EHG_rand_50.y_hat_2,1);
[~,~,~,AUC_rand_75] = perfcurve(results_EHG_rand_75.y_test_stored,results_EHG_rand_75.y_hat_2,1);
[~,~,~,AUC_rand_100] = perfcurve(results_EHG_rand_100.y_test_stored,results_EHG_rand_100.y_hat_2,1);


% rand_0 rand_25 rand_50 rand_75
AUCs_rand = [AUC_rand_0(1), AUC_rand_25(1), AUC_rand_50(1), AUC_rand_75(1), AUC_rand_100(1)];


CI_neg = [(AUC_rand_0(1)-AUC_rand_0(2)), (AUC_rand_25(1)-AUC_rand_25(2)), (AUC_rand_50(1)-AUC_rand_50(2)),...
    (AUC_rand_75(1)-AUC_rand_75(2)),(AUC_rand_100(1)-AUC_rand_100(2)) ];
CI_pos = [(AUC_rand_0(3)-AUC_rand_0(1)), (AUC_rand_25(3)-AUC_rand_25(1)), (AUC_rand_50(3)-AUC_rand_50(1)),...
    (AUC_rand_75(3)-AUC_rand_75(1)),(AUC_rand_100(3)-AUC_rand_100(1)) ];


errorbar(1:5,AUCs_rand,CI_neg,CI_neg,':ok','MarkerSize',marker_size,...
    'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',linewidth,'CapSize',capsize);hold on
plot(1:5,AUCs_rand,':','LineWidth',linewidth,'Color','#0072BD');

xlabel({'Columns rearranged (%)','(d)'});

ylabel('AUC')

set(gca,'xtick', 1:5)
set(gca,'xticklabel',  {'0','25','50', '75', '100'});

set(gca,'ytick', [0.5:0.2:0.9])

xlim([0.5 5.5])
ylim([0.5 1])

annotation('doublearrow',[0.245 0.40],[0.65 0.65])

% 20 min
% shortened sample
ax3 = nexttile([1 2]);
stft(downsample_signal(:,1),fs,'Window',hamming(window_length),'Centered',0)

child_image_3 = get(ax3,'Children');
child_image_3.CData(:,1:10) = 70;
child_image_3.CData(:,end-25:end) = 70;



h = colorbar;
set(get(h,'label'),'string','Magnitude [dB]','FontName','Arial','Fontsize',12);

xl = xlim;
xticks = linspace(xl(1),xl(2),7);
set(gca,'xtick', xticks(1:end-1))
set(gca,'xticklabel', [0:5:29])
set(gca,'ytick', [0:2:4])

xlabel({'Time [min]','(e)'})
ylabel('Frequency [Hz]','FontSize',12)
title("")
ylim([0 5])

nexttile
results_EHG_t1 = load('results_EHG_t1.mat');
results_EHG_t5 = load('results_EHG_t5.mat');
results_EHG_t10 = load('results_EHG_t10.mat');
results_EHG_t20 = load('results_EHG_t20.mat');
results_EHG_t30 = load('EHG_classification.mat');


for i=1:length(results_EHG_t1.y_hat)
    results_EHG_t1.y_hat_2{i} = results_EHG_t1.y_hat{i}(:,2);
    results_EHG_t5.y_hat_2{i} = results_EHG_t5.y_hat{i}(:,2);
    results_EHG_t10.y_hat_2{i} = results_EHG_t10.y_hat{i}(:,2);
    results_EHG_t20.y_hat_2{i} = results_EHG_t20.y_hat{i}(:,2);
    results_EHG_t30.y_hat_2{i} = results_EHG_t30.y_hat_2{i};%(:,2);
    
end

% EHG
[~,~,~,AUC_t1] = perfcurve(results_EHG_t1.y_test_stored,results_EHG_t1.y_hat_2,1);
[~,~,~,AUC_t5] = perfcurve(results_EHG_t5.y_test_stored,results_EHG_t5.y_hat_2,1);
[~,~,~,AUC_t10] = perfcurve(results_EHG_t10.y_test_stored,results_EHG_t10.y_hat_2,1);
[~,~,~,AUC_t20] = perfcurve(results_EHG_t20.y_test_stored,results_EHG_t20.y_hat_2,1);
[~,~,~,AUC_t30] = perfcurve(results_EHG_t30.y_test_stored,results_EHG_t30.y_hat_2,1);


% t1 t5 t10 t20
AUCs = [AUC_t1(1), AUC_t5(1), AUC_t10(1), AUC_t20(1), AUC_t30(1)];


CI_neg = [(AUC_t1(1)-AUC_t1(2)), (AUC_t5(1)-AUC_t5(2)), (AUC_t10(1)-AUC_t10(2)), (AUC_t20(1)-AUC_t20(2)), (AUC_t30(1)-AUC_t30(2))  ];
CI_pos = [(AUC_t1(3)-AUC_t1(1)), (AUC_t5(3)-AUC_t5(1)), (AUC_t10(3)-AUC_t10(1)), (AUC_t20(3)-AUC_t20(1)), (AUC_t30(3)-AUC_t30(1)) ];

errorbar(1:5,AUCs,CI_neg,CI_neg,'ok','MarkerSize',marker_size,...
    'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',linewidth,'CapSize',capsize);hold on

plot(1:5,AUCs,':','LineWidth',linewidth,'Color','#0072BD');

xlabel({'EHG duration [min]','(f)'});ylabel('AUC')
% title('')

set(gca,'xtick', 1:5)
set(gca,'xticklabel',  {'1','5','10', '20','30'});

set(gca,'ytick', [0.5:0.2:0.9])

xlim([0.5 5.5])
ylim([0.5 1])



set(gcf,'color','w');

