clear all; close all; clc
% S2 Fig
% This script generates 2 Fig and calculates skewness values of the
% distributions.

rng('default')
addpath(genpath('utils'))

myDir = uigetdir(pwd,'Open database directory');
addpath(myDir)

signalFiles = dir(fullfile(myDir,'*m.mat'));
headerFiles = dir(fullfile(myDir,'*.hea'));
LabelsVector = {};

count_tpehg = 1;
count = 1;

for k = 1:length(headerFiles)
    
    [x_tmp,~,tmp_name] = read_data(k,myDir);
    
    if x_tmp(1) < 26*7
        k;
        continue
    end
    
    [all.X(count,:),all.y(count,:),names{count} ] = read_data(k,myDir);
    
    if tmp_name(6) ~='t' % sample in tpehg
        [tpehg.X(count,:),tpehg.y(count,:),names{count} ] = read_data(k,myDir);
        count_tpehg = count_tpehg+1;
    end   
    
    count = count+1;
    
end



%% figure


figure('defaultAxesFontName','Arial','DefaultAxesFontSize',14)
tiledlayout(4,2,'TileSpacing','loose');
% panel a
nexttile
histogram(all.X(:,1),20)

ylabel('Counts')
xlabel({'Gestational age' 'at recording [days]' '(a)'})


xlim([180 255])
ylim([0 55])

set(gca,'xtick', [190:20:260])
set(gca,'ytick',[0:20:40])

alpha =0.1;
% [H, pValue_gest_age_rec_all, W] = swtest(all.X(:,1), alpha);
sk = skewness(all.X(:,1));

%b
nexttile
histogram(tpehg.X(:,1),20);
xlabel({'Gestational age' 'at recording [days]' '(b)'})
xlim([180 255])
ylim([0 55])
set(gca,'xtick', [190:20:260])
set(gca,'ytick',[0:20:40])
sk = skewness(all.X(:,1));


% c
nexttile
histogram(all.y(:,2),20)

xlabel({'Gestational age' 'at delivery [days]' '(c)'})
ylabel('Counts')

xlim([200 300])
ylim([0 35])

set(gca,'xtick',[210:20:300])
set(gca,'ytick',[0:10:30])

[H, pValue_delivery_age_all, W] = swtest(all.y(:,2), alpha);
% Results=normalitytest(all.y(:,2)')
sk = skewness(all.y(:,2))

%d
nexttile
histogram(tpehg.y(:,2),20)
xlabel({'Gestational age' 'at delivery [days]' '(d)'})
xlim([200 300])
ylim([0 35])
set(gca,'xtick',[210:20:300])
set(gca,'ytick',[0:10:30])

[H, pValue_delivery_age_tpehg, W] = swtest(tpehg.y(:,2), alpha);
% Results=normalitytest(all.y(:,2)')
sk = skewness(tpehg.y(:,2))

%e
nexttile

y_preterm_idx = find(all.y(:,2)<259 );
y_preterm = all.y(y_preterm_idx,2);
histogram(y_preterm,15)
ylabel('Counts')
xlabel({'Gestational age at delivery' 'of preterms [days]' '(e)'})
min(y_preterm)
max(y_preterm)

xlim([205 260])
ylim([0 10])

set(gca,'xtick',[210:20:260])
set(gca,'ytick',[0:5:10])


%f
nexttile

y_preterm_idx = find(tpehg.y(:,2)<259 );
y_preterm = tpehg.y(y_preterm_idx,2);
histogram(y_preterm,15)
xlabel({'Gestational age at delivery' 'of preterms [days]' '(f)'})
min(y_preterm)
max(y_preterm)
xlim([205 260])
ylim([0 10])
set(gca,'xtick',[210:20:260])
set(gca,'ytick',[0:5:10])

%g
nexttile

y_term_idx = find(all.y(:,2)>258 );
y_term = all.y(y_term_idx,2);
histogram(y_term,15)
ylabel('Counts')
xlabel({'Gestational age at delivery' 'of terms [days]' '(g)'})
min(y_term)
max(y_term)
xlim([258 295])
ylim([0 25])
set(gca,'xtick',[260:20:300])
set(gca,'ytick',[0:10:20])

[H, pValue_delivery_age_all_term, W] = swtest(y_term, alpha,true);

%h
nexttile

y_term_idx = find(tpehg.y(:,2)>258 );
y_term = tpehg.y(y_term_idx,2);
histogram(y_term,15)
xlabel({'Gestational age at delivery' 'of terms [days]' '(h)'})
min(y_term)
max(y_term)
xlim([258 295])
ylim([0 25])
set(gca,'xtick',[260:20:300])
set(gca,'ytick',[0:10:20])

[H, pValue_y_term_tpehg, W] = swtest(y_term, alpha);


% annotation('textbox', [0.045, 0.89, 0.05, 0.06], 'String', "a",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.52, 0.89, 0.05, 0.06], 'String', "b",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.045, 0.67, 0.05, 0.06], 'String', "c",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.52, 0.67, 0.05, 0.06], 'String', "d",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.045, 0.44, 0.05, 0.06], 'String', "e",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.52, 0.44, 0.05, 0.06], 'String', "f",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.045, 0.21, 0.05, 0.06], 'String', "g",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
% annotation('textbox', [0.52, 0.21, 0.05, 0.06], 'String', "h",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
annotation('textbox', [0.21, 0.91, 0.5, 0.06], 'String', "All samples",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
annotation('textbox', [0.62, 0.91, 0.5, 0.06], 'String', "TPEHG DB samples",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')



