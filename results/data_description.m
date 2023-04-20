% Computed values in Table 1 of the manuscript, describing the clinical information in the dataset

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


count = 1;
for k = 1:length(headerFiles)
    
    
    [downsample_signal,t_downsample,is_preterm(count),name,gestational_age(count),delivery_age(count)] = preprocessFile_v3(k,myDir);
    [X(count,:),y(count,:),names{count} ] = read_data_2(k,myDir);
    % Exclude recordings done before 26 weeks of gestation
    if  gestational_age(count) < 26*7
        k;
        continue
    end
 
    count = count+1;  
end


find(y(:,1)==1)
names(y(:,1)==1)

%% table info

% number samples
number_samples = size(y,1);

% number of preterms
number_preterms = sum(y(:,1));
fraction_preterms = mean(y(:,1));


% continuous variables: gestational_age, age, weight

% median (interquantile range), AUC with Confidence interval, missing
% entries


% age

is_valid = ~isnan(X(:,2));
valid_ages = X(is_valid,2);

median_age = median(valid_ages);
quantiles_age = quantile(valid_ages,[0.25 0.50 0.75 ])
[A_age,Aci_age] = auc([y(is_valid,1)-0.5,valid_ages],0.05);
missing_entries_age = sum(isnan(X(:,2)));

% gestational_age


is_valid = ~isnan(X(:,1));
valid_gest_ages = X(is_valid,1);
median_gest_age = median(valid_gest_ages);
median_gest_age_weeks = floor(median_gest_age/7);
median_gest_age_days = mod((median_gest_age  - floor(median_gest_age/7)*7),7);


quantiles_gest_age_weeks = floor(quantile(valid_gest_ages,[0.25 0.50 0.75 ])/7);
quantiles_gest_age_days = mod(quantile(valid_gest_ages,[0.25 0.50 0.75 ])  - quantiles_gest_age_weeks*7 ,7)
missing_entries_gest_age = sum(isnan(X(:,1)));

% weight



is_valid = ~isnan(X(:,5));
valid_weight = X(is_valid,5);
median_weight = median(valid_weight);
quantiles_weight = quantile(valid_weight,[0.25 0.50 0.75 ])
[A_weight,Aci_weight] = auc([-(y(is_valid,1)-0.5),valid_weight],0.05);
missing_entries_weight = sum(isnan(X(:,5)));


% discrete variables

% parity
is_valid =  ~isnan(X(:,3));
parity = X(is_valid,3)>0;
number_parous = sum(parity);
fraction_parous = mean(parity)*100;
oddds_matrix = [sum(parity.*y(is_valid,1)), sum(parity.*(~y(is_valid,1)));...
    sum(~parity.*y(is_valid,1)), sum(~parity.*~y(is_valid,1)) ];

odds(oddds_matrix)
missing_parity =  sum(isnan(X(:,3)));

% abortions
is_valid =  ~isnan(X(:,4));
abortions = X(is_valid,4)>0;
number_abortions = sum(abortions);
fraction_abortions = mean(abortions)*100;


oddds_matrix = [sum(abortions.*y(is_valid,1)), sum(abortions.*(~y(is_valid,1)));...
    sum(~abortions.*y(is_valid,1)), sum(~abortions.*~y(is_valid,1)) ];

odds(oddds_matrix)

missing_abortions =  sum(isnan(X(:,4)));

% hypertension 

is_valid =  ~isnan(X(:,6));
hypertension = X(is_valid,6)>0;

number_hypertension = sum(hypertension);
fraction_hypertension = mean(hypertension);


oddds_matrix = [sum(hypertension.*y(is_valid,1)), sum(hypertension.*(~y(is_valid,1)));...
    sum(~hypertension.*y(is_valid,1)), sum(~hypertension.*~y(is_valid,1)) ];

% odds(oddds_matrix); skip due to a lack of diagnoses of hypertension

missing_hypertension =  sum(isnan(X(:,6)));

% diabetes 
is_valid =  ~isnan(X(:,7));
diabetes = X(is_valid,7)>0;

number_diabetes = sum(diabetes);
fraction_diabetes = mean(diabetes)*100;


oddds_matrix = [sum(diabetes.*y(is_valid,1)), sum(diabetes.*(~y(is_valid,1)));...
    sum(~diabetes.*y(is_valid,1)), sum(~diabetes.*~y(is_valid,1)) ];


% odds(oddds_matrix); skip due to a lack of diagnoses of diabetes

missing_diabetes =  sum(isnan(X(:,7)));

% placental_position ( anterior placenta = front)

is_valid =  ~isnan(X(:,8));
placenta = X(is_valid,8)>0;

number_placenta = sum(placenta);
fraction_placenta = mean(placenta)*100;

oddds_matrix = [sum(placenta.*y(is_valid,1)), sum(placenta.*(~y(is_valid,1)));...
    sum(~placenta.*y(is_valid,1)), sum(~placenta.*~y(is_valid,1)) ];


odds(oddds_matrix)

missing_placenta =  sum(isnan(X(:,8)));


% bleeding1

is_valid =  ~isnan(X(:,9));
bleeding1 = X(is_valid,9)>0;

number_bleeding1 = sum(bleeding1);
fraction_bleeding1 = mean(bleeding1)*100;
oddds_matrix = [sum(bleeding1.*y(is_valid,1)), sum(bleeding1.*(~y(is_valid,1)));...
    sum(~bleeding1.*y(is_valid,1)), sum(~bleeding1.*~y(is_valid,1)) ];

odds(oddds_matrix)
missing_bleeding1 =  sum(isnan(X(:,9)));

% bleeding2

is_valid =  ~isnan(X(:,10));
bleeding2 = X(is_valid,10)>0;

number_bleeding2 = sum(bleeding2);
fraction_bleeding2 = mean(bleeding2)*100;
oddds_matrix = [sum(bleeding2.*y(is_valid,1)), sum(bleeding2.*(~y(is_valid,1)));...
    sum(~bleeding2.*y(is_valid,1)), sum(~bleeding2.*~y(is_valid,1)) ];

odds(oddds_matrix)

missing_bleeding2 =  sum(isnan(X(:,10)));

% funneling

is_valid =  ~isnan(X(:,11));
funneling = X(is_valid,11)>0;
number_funneling = sum(funneling);
fraction_funneling = mean(funneling)*100;
oddds_matrix = [sum(funneling.*y(is_valid,1)), sum(funneling.*(~y(is_valid,1)));...
    sum(~funneling.*y(is_valid,1)), sum(~funneling.*~y(is_valid,1)) ];

odds(oddds_matrix)

missing_funneling =  sum(isnan(X(:,11)));


% smoker

is_valid =  ~isnan(X(:,12));
smoker = X(is_valid,12)>0;
number_smoker = sum(smoker);
fraction_smoker = mean(smoker)*100;
oddds_matrix = [sum(smoker.*y(is_valid,1)), sum(smoker.*(~y(is_valid,1)));...
    sum(~smoker.*y(is_valid,1)), sum(~smoker.*~y(is_valid,1)) ];

odds(oddds_matrix)

missing_smoker =  sum(isnan(X(:,12)));

%% figure timings - S1 Fig.

preterms = (y(:,1) ==1);
terms = (y(:,1) ==0);

figure('defaultAxesFontName','Arial','DefaultAxesFontSize',14)
tiledlayout(1,2,'TileSpacing','loose');
nexttile

plot(X(preterms,1),y(preterms,2),'o','MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize', 8);hold on

plot(X(terms,1),y(terms,2),'d','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize', 8);hold on

yline(258.5,'--k','Linewidth',2)
xlabel({'Gestational age at recording [days]','(a)'})
ylabel('Gestational age at birth [days]')

yl=get(gca,'ylabel');
pyl = get(yl,'position');
pyl(1) = 0.9*pyl(1);
set(yl,'position',pyl)
    
xlim([180 255])
ylim([205 300])

set(gca,'xtick', [190:20:260])
set(gca,'xticklabel', [190:20:260])


set(gca,'ytick',[210:20:300])
set(gca,'yticklabel',[210:20:300])

annotation('textbox', [0.38, 0.75 0.05, 0.06], 'String', "Term",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
annotation('textbox', [0.37, 0.4 0.05, 0.06], 'String', "Preterm",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')

nexttile

plot(X(preterms,1),y(preterms,2)-X(preterms,1),'o','MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319','MarkerSize', 8);hold on
plot(X(terms,1),y(terms,2)-X(terms,1),'d','MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD','MarkerSize', 8);hold on

x_points = [180 255];
y_points = 258.5- x_points;
plot(x_points,y_points,'--k','Linewidth',2)

xlabel({'Gestational age at recording [days]','(b)'})
ylabel({'Time from recording' 'to birth [days]'})

xlim([180 255])
ylim([-3 84])

set(gca,'xtick', [190:20:260])
set(gca,'xticklabel', [190:20:260])


set(gca,'ytick',[0:20:80])


annotation('textbox', [0.83, 0.75 0.05, 0.06], 'String', "Term",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')
annotation('textbox', [0.58, 0.4 0.05, 0.06], 'String', "Preterm",'FontName','Arial','FontSize',14,'LineStyle','none','fontweight', 'bold')

set(gcf,'color','w');
