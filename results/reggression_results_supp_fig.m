clear all; close all; clc
rng('default')
addpath('results/roc_aucs')
addpath(genpath('utils'))

%metadata
results_reg_meta = load('regression_metadata.mat');

%EHG
results_reg_EHG = load('EHG_regression.mat');

%combined
results_reg_combined = load('regression_combined.mat');


%% Open figure
linewidth = 2;
legend_font_size = 11;

figure('defaultAxesFontName','Arial','DefaultAxesFontSize',13.5)
tiledlayout(3,2,'TileSpacing','loose');

%% Clinical information

y_hat_cat = [];
y_test_cat = [];

for i = 1:5
    y_hat_cat = [y_hat_cat;results_reg_meta.y_hat_cell{i}];
    y_test_cat = [y_test_cat;results_reg_meta.y_test_cell{i}];
    
end

nexttile
plot(y_test_cat,y_hat_cat,'o');hold on;

min_val = min([y_test_cat;y_hat_cat]);
max_val = max([y_test_cat;y_hat_cat]);

min_val = 200;
max_val = 300;
xlim([200 300])
ylim([200 300])

xlabel('(a)')
ylabel({'Predicted delivery' 'age [days]'})

[p] = polyfit(y_test_cat,y_hat_cat,1);
fit_line_1 = p(1)*min_val+p(2);
fit_line_2 = p(1)*max_val+p(2);

plot([min_val,max_val ],[fit_line_1,fit_line_2 ],'LineWidth',linewidth,'Color','k'); hold on
plot([min_val,max_val ],[min_val,max_val ],'--','LineWidth',linewidth,'Color','k'); hold on

set(gca,'xtick', [200:20:300])
set(gca,'ytick',[200:40:300])


count = 1;
for fold = 1:5:96
    
    y_hat_cat =[];
    y_test_cat =[];
    
    for i = 0:4
        y_hat_cat = [y_hat_cat;results_reg_meta.y_hat_cell{fold+i}];
        y_test_cat = [y_test_cat;results_reg_meta.y_test_cell{fold+i}];
        
    end
    
    mdl = fitlm(y_test_cat,y_hat_cat);
    R2(count) = mdl.Rsquared.Ordinary;
    RMSE(count) = mdl.RMSE;
    slope(count) = mdl.Coefficients.Estimate(2);
    
    count= count+1;
end

mean(R2);
ci_r2 = CI95(R2);

mean(RMSE);
ci_RMSE = CI95(RMSE);

mean(slope);
ci_slope = CI95(slope);

str = {['RMSE = ', num2str(round(mean(RMSE),2)), ' (', num2str(round(ci_RMSE(1),2)), ', ',num2str(round(ci_RMSE(2),2)), ') days'],...
    ['R^2 = ', num2str(round(mean(R2),2)), ' (', sprintf('%.2f',round(ci_r2(1),2)), ', ',num2str(round(ci_r2(2),2)), ')' ],...
    ['Slope = ' num2str(round(mean(slope),2)), ' (', num2str(round(ci_slope(1),2)), ', ',num2str(round(ci_slope(2),2)), ')' ]};

annotation('textbox', [0.28, 0.72, 0.5, 0.06], 'String', str,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')

% Bland Altman plot
nexttile

y_hat_cat =[];
y_test_cat =[];

for i = 1:5
    y_hat_cat = [y_hat_cat;results_reg_meta.y_hat_cell{i}];
    y_test_cat = [y_test_cat;results_reg_meta.y_test_cell{i}];
    
end


means = mean([y_test_cat,y_hat_cat],2);
diffs = y_test_cat -y_hat_cat;


plot(means,diffs,'o');hold on;

min_val_x = min(means);
max_val_x = max(means);

min_val_y = min(diffs);
max_val_y = max(diffs);

mean_diffs = mean(diffs);
two_Stds = 1.96*std(diffs);

plot([200 300],[mean_diffs, mean_diffs],'LineWidth',linewidth,'Color','k');hold on

plot([200 300],[mean_diffs+two_Stds, mean_diffs+two_Stds],'--','LineWidth',linewidth,'Color','k');hold on
plot([200 300],[mean_diffs-two_Stds, mean_diffs-two_Stds],'--','LineWidth',linewidth,'Color','k');hold on


xlim([200 290])
ylim([-70 55])

set(gca,'xtick', [200:20:300])
set(gca,'ytick',[-60:40:41])
xlabel('(b)')

ylabel({'Delivery age - predicted' 'delivery age [days]'})


count = 1;
for fold = 1:5:96
    
    y_hat_cat =[];
    y_test_cat =[];
    
    for i = 0:4
        y_hat_cat = [y_hat_cat;results_reg_meta.y_hat_cell{fold+i}];
        y_test_cat = [y_test_cat;results_reg_meta.y_test_cell{fold+i}];
        
    end
    
    diffs = y_test_cat -y_hat_cat;
    
    mean_diffs(count) = mean(diffs);
    mean_plus_two_Stds(count) = mean(diffs) + 1.96*std(diffs);
    mean_minus_two_Stds(count) = mean(diffs) - 1.96*std(diffs);
    
    count= count+1;
end

mean(mean_diffs);
ci_mean = CI95(mean_diffs);

mean(mean_plus_two_Stds);
ci_mean_plus = CI95(mean_plus_two_Stds);

mean(mean_minus_two_Stds);
ci_mean_minus = CI95(mean_minus_two_Stds);

str1 = ['+1.96 s.d.: ',num2str(round(mean(mean_plus_two_Stds),2)), ' (',num2str(round(ci_mean_plus(1),2)), ', ',num2str(round(ci_mean_plus(2),2)), ')'];
str2 = ['Mean: ',num2str(round(mean(mean_diffs),2)), ' (', sprintf('%.2f',(round(ci_mean(1),2))), ', ',num2str(round(ci_mean(2),2)), ')'];

str_3_1 = '-1.96 s.d.: ';
str_3_2 = [num2str(round(mean(mean_minus_two_Stds),2)), ' (',sprintf('%.2f',(round(ci_mean_minus(1),2))), ', ',sprintf('%.2f',(round(ci_mean_minus(2),2))), ')'];
str3 = [{str_3_1},{str_3_2}];


annotation('textbox', [0.57, 0.855, 0.5, 0.06], 'String', str1,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')
annotation('textbox', [0.57, 0.795, 0.5, 0.06], 'String', str2,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')
annotation('textbox', [0.57, 0.75, 0.5, 0.06], 'String', str3,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')



%% EHG

y_hat_cat =[];
y_test_cat =[];

for i = 1:5
    y_hat_cat = [y_hat_cat;results_reg_EHG.y_hat{i}];
    y_test_cat = [y_test_cat;results_reg_EHG.y_test_stored{i}];
    
end

nexttile
plot(y_test_cat,y_hat_cat,'o');hold on;

min_val = min([y_test_cat;y_hat_cat]);
max_val = max([y_test_cat;y_hat_cat]);


min_val = 200;
max_val = 300;

xlim([200 300])
ylim([200 300])

xlabel('(c)')
ylabel({'Predicted delivery' 'age [days]'})
[p] = polyfit(y_test_cat,y_hat_cat,1);
fit_line_1 = p(1)*min_val+p(2);
fit_line_2 = p(1)*max_val+p(2);

plot([min_val,max_val ],[fit_line_1,fit_line_2 ],'LineWidth',linewidth,'Color','k'); hold on
plot([min_val,max_val ],[min_val,max_val ],'--','LineWidth',linewidth,'Color','k'); hold on

set(gca,'xtick', [200:20:300])
set(gca,'ytick',[200:40:300])


count = 1;
for fold = 1:5:96
    
    y_hat_cat =[];
    y_test_cat =[];
    
    for i = 0:4
        y_hat_cat = [y_hat_cat;results_reg_EHG.y_hat{fold+i}];
        y_test_cat = [y_test_cat;results_reg_EHG.y_test_stored{fold+i}];
        
    end
    
    mdl = fitlm(y_test_cat,y_hat_cat);
    R2(count) = mdl.Rsquared.Ordinary;
    RMSE(count) = mdl.RMSE;
    slope(count) = mdl.Coefficients.Estimate(2);
    
    count= count+1;
end

mean(R2);
ci_r2 = CI95(R2);

mean(RMSE);
ci_RMSE = CI95(RMSE);

mean(slope);
ci_slope = CI95(slope);


str = {['RMSE = ', num2str(round(mean(RMSE),2)), ' (', num2str(round(ci_RMSE(1),2)), ', ',num2str(round(ci_RMSE(2),2)), ') days'],...
    ['R^2 = ', num2str(round(mean(R2),2)), ' (', num2str(round(ci_r2(1),2)), ', ',num2str(round(ci_r2(2),2)), ')' ],...
    ['Slope = ' num2str(round(mean(slope),2)), ' (', num2str(round(ci_slope(1),2)), ', ',num2str(round(ci_slope(2),2)), ')' ]};

annotation('textbox', [0.28, 0.42, 0.5, 0.06], 'String', str,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')

% Bland altman plot

nexttile

y_hat_cat =[];
y_test_cat =[];

for i = 1:5
    y_hat_cat = [y_hat_cat;results_reg_EHG.y_hat{i}];
    y_test_cat = [y_test_cat;results_reg_EHG.y_test_stored{i}];
    
end

means = mean([y_test_cat,y_hat_cat],2);
diffs = y_test_cat -y_hat_cat;


plot(means,diffs,'o');hold on;

min_val_x = min(means);
max_val_x = max(means);

min_val_y = min(diffs);
max_val_y = max(diffs);

mean_diffs = mean(diffs);
two_Stds = 1.96*std(diffs);

plot([200 300],[mean_diffs, mean_diffs],'LineWidth',linewidth,'Color','k');hold on

plot([200 300],[mean_diffs+two_Stds, mean_diffs+two_Stds],'--','LineWidth',linewidth,'Color','k');hold on
plot([200 300],[mean_diffs-two_Stds, mean_diffs-two_Stds],'--','LineWidth',linewidth,'Color','k');hold on


xlim([200 290])

ylim([-70 55])
xlabel('(d)')

ylabel({'Delivery age - predicted' 'delivery age [days]'})


count = 1;
for fold = 1:5:96
    
    y_hat_cat =[];
    y_test_cat =[];
    
    for i = 0:4
        y_hat_cat = [y_hat_cat;results_reg_EHG.y_hat{fold+i}];
        y_test_cat = [y_test_cat;results_reg_EHG.y_test_stored{fold+i}];
        
    end
    
    diffs = y_test_cat -y_hat_cat;
    
    mean_diffs(count) = mean(diffs);
    mean_plus_two_Stds(count) = mean(diffs) + 1.96*std(diffs);
    mean_minus_two_Stds(count) = mean(diffs) - 1.96*std(diffs);
    
    count= count+1;
end

mean(mean_diffs);
ci_mean = CI95(mean_diffs);

mean(mean_plus_two_Stds);
ci_mean_plus = CI95(mean_plus_two_Stds);


mean(mean_minus_two_Stds);
ci_mean_minus = CI95(mean_minus_two_Stds);


str1 = ['+1.96 s.d.: ',num2str(round(mean(mean_plus_two_Stds),2)), ' (',num2str(round(ci_mean_plus(1),2)), ', ',num2str(round(ci_mean_plus(2),2)), ')'];
str2 = ['Mean: ',num2str(round(mean(mean_diffs),2)), ' (',num2str(round(ci_mean(1),2)), ', ',num2str(round(ci_mean(2),2)), ')'];

str3_1 = '-1.96 s.d.:';
str3_2 = [num2str(round(mean(mean_minus_two_Stds),2)), ' (',num2str(round(ci_mean_minus(1),2)), ', ',num2str(round(ci_mean_minus(2),2)), ')'];
str3 = [{str3_1},{str3_2}];

annotation('textbox', [0.57, 0.555, 0.5, 0.06], 'String', str1,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')
annotation('textbox', [0.57, 0.5, 0.5, 0.06], 'String', str2,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')
annotation('textbox', [0.57, 0.46, 0.5, 0.06], 'String', str3,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')

%% Combined


y_hat_cat =[];
y_test_cat =[];

for i = 1:5
    y_hat_cat = [y_hat_cat;results_reg_combined.y_hat_stored{i}];
    y_test_cat = [y_test_cat;results_reg_combined.y_test_stored{i}];
    
end

nexttile
plot(y_test_cat,y_hat_cat,'o');hold on;

min_val = min([y_test_cat;y_hat_cat]);
max_val = max([y_test_cat;y_hat_cat]);

min_val = 200;
max_val = 300;


xlim([200 300])
ylim([200 300])

xlabel({'Delivery age [days]','(e)'})
ylabel({'Predicted delivery' 'age [days]'})

[p] = polyfit(y_test_cat,y_hat_cat,1);

fit_line_1 = p(1)*min_val+p(2);
fit_line_2 = p(1)*max_val+p(2);

plot([min_val,max_val ],[fit_line_1,fit_line_2 ],'LineWidth',linewidth,'Color','k'); hold on
plot([min_val,max_val ],[min_val,max_val ],'--','LineWidth',linewidth,'Color','k'); hold on

set(gca,'xtick', [200:20:300])
set(gca,'ytick',[200:40:300])


count = 1;
for fold = 1:5:96
    
    y_hat_cat =[];
    y_test_cat =[];
    
    for i = 0:4
        y_hat_cat = [y_hat_cat;results_reg_combined.y_hat_stored{fold+i}];
        y_test_cat = [y_test_cat;results_reg_combined.y_test_stored{fold+i}];
        
    end
    
    mdl = fitlm(y_test_cat,y_hat_cat);
    R2(count) = mdl.Rsquared.Ordinary;
    RMSE(count) = mdl.RMSE;
    slope(count) = mdl.Coefficients.Estimate(2);
    
    count= count+1;
end

mean(R2);
ci_r2 = CI95(R2);

mean(RMSE);
ci_RMSE = CI95(RMSE);

mean(slope);
ci_slope = CI95(slope);


str = {['RMSE = ', num2str(round(mean(RMSE),2)), ' (', num2str(round(ci_RMSE(1),2)), ', ',num2str(round(ci_RMSE(2),2)), ') days'],...
    ['R^2 = ', num2str(round(mean(R2),2)), ' (', num2str(round(ci_r2(1),2)), ', ',num2str(round(ci_r2(2),2)), ')' ],...
    ['Slope = ' num2str(round(mean(slope),2)), ' (', num2str(round(ci_slope(1),2)), ', ',num2str(round(ci_slope(2),2)), ')' ]};


annotation('textbox', [0.28, 0.12, 0.5, 0.06], 'String', str,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')

% Bland Altman plot

nexttile

y_hat_cat =[];
y_test_cat =[];

for i = 1:5
    y_hat_cat = [y_hat_cat;results_reg_combined.y_hat_stored{i}];
    y_test_cat = [y_test_cat;results_reg_combined.y_test_stored{i}];
    
end

means = mean([y_test_cat,y_hat_cat],2);
diffs = y_test_cat -y_hat_cat;


plot(means,diffs,'o');hold on;

min_val_x = min(means);
max_val_x = max(means);

min_val_y = min(diffs);
max_val_y = max(diffs);

mean_diffs = mean(diffs);
two_Stds = 1.96*std(diffs);

plot([200 300],[mean_diffs, mean_diffs],'LineWidth',linewidth,'Color','k');hold on

plot([200 300],[mean_diffs+two_Stds, mean_diffs+two_Stds],'--','LineWidth',linewidth,'Color','k');hold on
plot([200 300],[mean_diffs-two_Stds, mean_diffs-two_Stds],'--','LineWidth',linewidth,'Color','k');hold on



xlim([200 290])

ylim([-70 55])

xlabel({'Mean delivery age and' 'predicted delivery age [days]','(f)'})
ylabel({'Delivery age - predicted' 'delivery age [days]'})


count = 1;
for fold = 1:5:96
    
    y_hat_cat =[];
    y_test_cat =[];
    
    for i = 0:4
        y_hat_cat = [y_hat_cat;results_reg_combined.y_hat_stored{fold+i}];
        y_test_cat = [y_test_cat;results_reg_combined.y_test_stored{fold+i}];
        
    end
    
    diffs = y_test_cat -y_hat_cat;
    
    mean_diffs(count) = mean(diffs);
    mean_plus_two_Stds(count) = mean(diffs) + 1.96*std(diffs);
    mean_minus_two_Stds(count) = mean(diffs) - 1.96*std(diffs);
    
    count= count+1;
end

mean(mean_diffs);
ci_mean = CI95(mean_diffs);

mean(mean_plus_two_Stds);
ci_mean_plus = CI95(mean_plus_two_Stds);


mean(mean_minus_two_Stds);
ci_mean_minus = CI95(mean_minus_two_Stds);



str1 = ['+1.96 s.d.: ',num2str(round(mean(mean_plus_two_Stds),2)), ' (',num2str(round(ci_mean_plus(1),2)), ', ',num2str(round(ci_mean_plus(2),2)), ')'];
str2 = ['Mean: ',num2str(round(mean(mean_diffs),2)), ' (',num2str(round(ci_mean(1),2)), ', ',num2str(round(ci_mean(2),2)), ')'];

str3_1 = '-1.96 s.d.:';
str_3_2 = [ num2str(round(mean(mean_minus_two_Stds),2)), ' (',num2str(round(ci_mean_minus(1),2)), ', ',num2str(round(ci_mean_minus(2),2)), ')'];
str3 = [{str3_1},{str_3_2}];

annotation('textbox', [0.57, 0.26, 0.5, 0.06], 'String', str1,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')
annotation('textbox', [0.57, 0.21, 0.5, 0.06], 'String', str2,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')
annotation('textbox', [0.57, 0.17, 0.5, 0.06], 'String', str3,'FontName','Arial','FontSize',legend_font_size,'LineStyle','none')

annotation('textbox', [0.43, 0.92 0.2, 0.06], 'String', "Metadata alone",'FontName','Arial','FontSize',16,'LineStyle','none','fontweight', 'bold')
annotation('textbox', [0.45, 0.61 0.2, 0.06], 'String', "EHG alone",'FontName','Arial','FontSize',16,'LineStyle','none','fontweight', 'bold')
annotation('textbox', [0.43, 0.32 0.2, 0.06], 'String', "Metadata and EHG",'FontName','Arial','FontSize',16,'LineStyle','none','fontweight', 'bold')

font_size =17;