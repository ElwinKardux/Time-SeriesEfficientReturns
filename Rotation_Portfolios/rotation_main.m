%% Main Rotation portfolio's 
rot = load('rotation_data4.mat');
rot2 = load('rotation_data.mat');
rot_MA = load('rotation_MA.mat');
PCA = load('rot_PCA.mat');
eff_PC = PCA.eff_PC;
rot_PC_eff = rot2.rot_PC_eff;
rot_PC_MA = PCA.rot_PC_MA;
rot_PC_PM = PCA.rot_PC_PM;

%important
rot_old = rot2.rot_old;
rot_ar1_eff = rot2.rot_ar1_eff;
rot_ar12_eff = rot2.rot_ar12_eff;
rot_cross_eff = rot2.rot_cross_eff;

%rot_old_MA = rot.rot_MA_old;
%rot_ar1_MA = rot.rot_MA_eff;
rot_old_MA = rot_MA.rot_old_MA;
rot_ar1_MA = rot_MA.rot_eff_MA;
rot_ar12_MA = rot_MA.rot_ar12_MA;
rot_cross_MA = rot_MA.rot_cross_MA;
rot_old_PM = rot.rot_PM_old;
rot_ar1_PM = rot.rot_PM_ar1;
rot_ar12_PM = rot.rot_ar12_PM;
rot_cross_PM = rot.rot_cross_PM;


n_mom = length(rot_ar12_MA);
n_eff = length(rot_ar1_eff);

%% Efficient 
eff = load('efficient_data.mat');
original = eff.a_org;
eff_inter = eff.a_ar1;
eff_prior = eff.a_ar12;
eff_cross = eff.a_cross;

dataFF5 = csvread("FF5.CSV",4,1);
ff5oud = dataFF5(1:end,1:5);
ff5oudX = [ones(size(ff5oud,1),1) ff5oud];
Rf = dataFF5(1:end,6);
SR_oldFF5 = zeros(1,5);
SR_newFF5 = zeros(1,5);
ff5eff = zeros(size(ff5oud));
XiFF5 = zeros(size(ff5oud));
for i=1:5
    [ff5eff(:,i),XiFF5(:,i), SR_newFF5(i), SR_oldFF5(i)] = efficient(ff5oud(:,i));
end
%ff5effX = [ones(size(ff5eff,1),1) ff5eff];
ff5effX = ff5eff(end-n_eff +1:end,:);

%% OLS regression
% ALL stands for the necessary values for the tables on Alpha Analysis.
ALL = zeros(24,7);
% Prevailing Mean
model_old_PM= fitlm(ff5effX, rot_old_PM);
ALL(1:2,1:6) = table2array(model_old_PM.Coefficients(:,1:2))';
ALL(1,7) = model_old_PM.Rsquared.Ordinary;

model_ar1_PM= fitlm(ff5effX, rot_ar1_PM);
ALL(3:4,1:6) = table2array(model_ar1_PM.Coefficients(:,1:2))';
ALL(3,7) = model_ar1_PM.Rsquared.Ordinary;

model_ar12_PM= fitlm(ff5effX, rot_ar12_PM);
ALL(5:6,1:6) = table2array(model_ar12_PM.Coefficients(:,1:2))';
ALL(5,7) = model_ar12_PM.Rsquared.Ordinary;

model_cross_PM= fitlm(ff5effX, rot_cross_PM);
ALL(7:8,1:6) = table2array(model_cross_PM.Coefficients(:,1:2))';
ALL(7,7) = model_cross_PM.Rsquared.Ordinary;
%Adaptive Lasso

model_old= fitlm(ff5effX, rot_old);
ALL(9:10,1:6) = table2array(model_old.Coefficients(:,1:2))';
ALL(9,7) = model_old.Rsquared.Ordinary;
% 
model_ar1_eff= fitlm(ff5effX, rot_ar1_eff);
ALL(11:12,1:6) = table2array(model_ar1_eff.Coefficients(:,1:2))';
ALL(11,7) = model_ar1_eff.Rsquared.Ordinary;

model_ar12_eff= fitlm(ff5effX, rot_ar12_eff);
ALL(13:14,1:6) = table2array(model_ar12_eff.Coefficients(:,1:2))';
ALL(13,7) = model_ar12_eff.Rsquared.Ordinary;

model_cross_eff= fitlm(ff5effX, rot_cross_eff);
ALL(15:16,1:6) = table2array(model_cross_eff.Coefficients(:,1:2))';
ALL(15,7) = model_cross_eff.Rsquared.Ordinary;


% Cross-sectional
model_old_MA= fitlm(ff5effX, rot_old_MA);
ALL(17:18,1:6) = table2array(model_old_MA.Coefficients(:,1:2))';
ALL(17,7) = model_old_PM.Rsquared.Ordinary;

model_ar1_MA= fitlm(ff5effX, rot_ar1_MA);
ALL(19:20,1:6) = table2array(model_ar1_MA.Coefficients(:,1:2))';
ALL(19,7) = model_ar1_MA.Rsquared.Ordinary;

model_ar12_MA= fitlm(ff5effX, rot_ar12_MA);
ALL(21:22,1:6) = table2array(model_ar12_MA.Coefficients(:,1:2))';
ALL(21,7) = model_ar12_PM.Rsquared.Ordinary;

model_cross_MA= fitlm(ff5effX, rot_cross_MA);
ALL(23:24,1:6) = table2array(model_cross_MA.Coefficients(:,1:2))';
ALL(23,7) = model_cross_MA.Rsquared.Ordinary;

% Annualizing the alpha's
ALL(:,1) = ALL(:,1)*12;

%% Recession: 
% Computes the Recession Analysis results
VIX = csvread("USREC.csv",1,1);
%VI = quantile(VIX,0.80);
%V = VIX >= VI;
V = VIX;
NAI = csvread("NAI.csv",1,1);
NA = quantile(NAI,0.20);
N = NAI <= NA;
B = zeros(8,6);
%% Prevailing Mean
model_old_PM= fitlm(V, rot_old_PM);
B(1:2,1) = table2array(model_old_PM.Coefficients(2,1:2))';
model_old_PM= fitlm(N, rot_old_PM);
B(1:2,2) = table2array(model_old_PM.Coefficients(2,1:2))';

%% 
model_ar1_PM= fitlm(V, rot_ar1_PM);
B(3:4,1) = table2array(model_ar1_PM.Coefficients(2,1:2))';
model_ar1_PM= fitlm(N, rot_ar1_PM);
B(3:4,2) = table2array(model_ar1_PM.Coefficients(2,1:2))';
%%
model_ar12_PM= fitlm(V, rot_ar12_PM);
B(5:6,1) = table2array(model_ar12_PM.Coefficients(2,1:2))';
model_ar12_PM= fitlm(N, rot_ar12_PM);
B(5:6,2) = table2array(model_ar12_PM.Coefficients(2,1:2))';
%%
model_cross_PM= fitlm(V, rot_cross_PM);
B(7:8,1) = table2array(model_cross_PM.Coefficients(2,1:2))';
model_cross_PM= fitlm(N, rot_cross_PM);
B(7:8,2) = table2array(model_cross_PM.Coefficients(2,1:2))';
%%
% Adaptive Lasso
model_old= fitlm(V, rot_old+0.25);
B(1:2,3) = table2array(model_old.Coefficients(2,1:2))';
model_old= fitlm(N, rot_old+0.25);
B(1:2,4) = table2array(model_old.Coefficients(2,1:2))';
%%
model_ar1_eff= fitlm(V, rot_ar1_eff);
B(3:4,3) = table2array(model_ar1_eff.Coefficients(2,1:2))';
model_ar1_eff= fitlm(N, rot_ar1_eff);
B(3:4,4) = table2array(model_ar1_eff.Coefficients(2,1:2))';
%%
model_ar12_eff= fitlm(V, rot_ar12_eff);
B(5:6,3) = table2array(model_ar12_eff.Coefficients(2,1:2))';
model_ar12_eff= fitlm(N, rot_ar12_eff);
B(5:6,4) = table2array(model_ar12_eff.Coefficients(2,1:2))';
%%
model_cross_eff= fitlm(V, rot_cross_eff);
B(7:8,3) = table2array(model_cross_eff.Coefficients(2,1:2))';
model_cross_eff= fitlm(N, rot_cross_eff);
B(7:8,4) = table2array(model_cross_eff.Coefficients(2,1:2))';

%%

% Cross-sectional
model_old_MA= fitlm(V, rot_old_MA);
B(1:2,5) = table2array(model_old_MA.Coefficients(2,1:2))';
model_old_MA= fitlm(N, rot_old_MA);
B(1:2,6) = table2array(model_old_MA.Coefficients(2,1:2))';
%%
model_ar1_MA= fitlm(V, rot_ar1_MA);
B(3:4,5) = table2array(model_ar1_MA.Coefficients(2,1:2))';
model_ar1_MA= fitlm(N, rot_ar1_MA);
B(3:4,6) = table2array(model_ar1_MA.Coefficients(2,1:2))';
%%
model_ar12_MA= fitlm(V, rot_ar12_MA);
B(5:6,5) = table2array(model_ar12_MA.Coefficients(2,1:2))';
model_ar12_MA= fitlm(N, rot_ar12_MA);
B(5:6,6) = table2array(model_ar12_MA.Coefficients(2,1:2))';
%%
model_cross_MA= fitlm(V, rot_cross_MA);
B(7:8,5) = table2array(model_cross_MA.Coefficients(2,1:2))';
model_cross_MA= fitlm(N, rot_cross_MA);
B(7:8,6) = table2array(model_cross_MA.Coefficients(2,1:2))';

% Principal components

ALL2 = zeros(6,7);
% Prevailing Mean
model_PC_PM= fitlm(ff5effX, rot_PC_PM);
ALL2(1:2,1:6) = table2array(model_PC_PM.Coefficients(:,1:2))';
ALL2(1,7) = model_PC_PM.Rsquared.Ordinary;

model_PC_eff= fitlm(ff5effX, rot_PC_eff);
ALL2(3:4,1:6) = table2array(model_PC_eff.Coefficients(:,1:2))';
ALL2(3,7) = model_PC_eff.Rsquared.Ordinary;

model_PC_MA= fitlm(ff5effX, rot_PC_MA);
ALL2(5:6,1:6) = table2array(model_PC_MA.Coefficients(:,1:2))';
ALL2(5,7) = model_PC_MA.Rsquared.Ordinary;

ALL2(:,1) = ALL2(:,1)*12;

PCrec = zeros(2,6);
model_PC_PM= fitlm(V, rot_PC_PM);
PCrec(1:2,1) = table2array(model_PC_PM.Coefficients(2,1:2))';
model_PC_PM= fitlm(N, rot_PC_PM);
PCrec(1:2,2) = table2array(model_PC_PM.Coefficients(2,1:2))';
%
model_PC_eff = fitlm(V, rot_PC_eff);
PCrec(1:2,3) = table2array(model_PC_eff.Coefficients(2,1:2))';
model_PC_eff = fitlm(N, rot_PC_eff);
PCrec(1:2,4) = table2array(model_PC_eff.Coefficients(2,1:2))';
%
model_PC_MA= fitlm(V, rot_PC_MA);
PCrec(1:2,5) = table2array(model_PC_MA.Coefficients(2,1:2))';
model_PC_PM= fitlm(N, rot_PC_MA);
PCrec(1:2,6) = table2array(model_PC_PM.Coefficients(2,1:2))';

%% 
%cumplot4(rot_old+0.25, rot_ar1_eff+0.25, rot_ar12_eff+0.25, rot_cross_eff +0.25, rot_PC_eff + 0.25, '(B) Adaptive Lasso Rotation Portfolio')
%cumplot4(rot_old_MA, rot_ar12_MA, rot_ar1_MA, rot_cross_MA , rot_PC_MA,'(C) Cross-Sectional Momentum Rotation Portfolio')
%cumplot4(rot_old_PM, rot_ar1_PM, rot_ar12_PM, rot_cross_PM,rot_PC_PM, '(A) Prevailing Mean')

function cumplot2(returns, returns2, returns3, titles)
returns11 = returns/std(returns);
returns22 = returns2/std(returns2);
returns33 = returns3/std(returns3);

a = log(cumprod(returns11/100 + 1.0));
a2 = log(cumprod(returns22/100 + 1));
a3 = log(cumprod(returns33/100 + 1));

b = csvread("12_indust_month_value.csv",1 ,0);
b2 = b(end-423+1:end,1);
time = datenum(num2str(b2),'yyyymm');
annRet = num2str(round(mean(returns+ 0)*12,2));
annRet2 = num2str(round(mean(returns2+ 0)*12,2));
annRet3 = num2str(round(mean(returns3+ 0)*12,2));
display(strcat("Annual. Return Rotation eff = ", num2str(annRet)));
display(strcat("Annual. Return Moving Average = ", num2str(annRet2)));
display(strcat("Annual. Return Market factor = ", num2str(annRet3)));

display(strcat("Std Rotation eff = ", num2str(sqrt(12)*std(returns))));
display(strcat("Std Moving Average = ", num2str(sqrt(12)*std(returns2))));
display(strcat("Std Market Factor = ", num2str(sqrt(12)*std(returns3))));

display(strcat("Max DD Rotation eff = ", num2str(round(maxdrawdown(cumprod(returns/100 + 1.0))*100,2))));
display(strcat("Max DD Moving Average = ", num2str(round(maxdrawdown(cumprod(returns2/100 + 1.0))*100,2))));
display(strcat("Max DD Market Factor = ", num2str(round(maxdrawdown(cumprod(returns3/100 + 1.0))*100,2))));
figure(gcf);
plot(time,a,'b', 'LineWidth', 2 );hold on;
plot(time,a2,'r', 'LineWidth', 2);
plot(time,a3,'k', 'LineWidth', 2);
plot([min(xlim()),max(xlim())],[0,0], 'k--');
ax = gca; 
ylim([-0.15, 1.15]);
axis 'auto x'
datetick('x', 'yyyy', 'keepticks')

title(titles);
ax.FontSize = 12;
recessionplot;
set(gcf,'units','points','position',[10,10,400,250])
legend({'Adaptive Lasso','Moving Average', 'All Ind. Equal Weighted'},'Location','northwest','FontSize',10)
end


function cumplot4(returns, returns2, returns3, returns4, returns5, titles)
returns11 = returns/std(returns);
returns22 = returns2/std(returns2);
returns33 = returns3/std(returns3);
returns44 = returns4/std(returns4);
returns55 = returns5/std(returns5);


a = log(cumprod(returns11/100 + 1.0));
a2 = log(cumprod(returns22/100 + 1));
a3 = log(cumprod(returns33/100 + 1));
a4 = log(cumprod(returns44/100 + 1));
a5 = log(cumprod(returns55/100 + 1));
b = csvread("12_indust_month_value.csv",1 ,0);
b2 = b(end-423+1:end,1);
time = datenum(num2str(b2),'yyyymm');

annRet = num2str(round(mean(returns+ 0)*12,2));
annRet2 = num2str(round(mean(returns2+ 0)*12,2));
annRet3 = num2str(round(mean(returns3+ 0)*12,2));
annRet4 = num2str(round(mean(returns4+ 0)*12,2));
annRet5 = num2str(round(mean(returns5+ 0)*12,2));

display(strcat("Annual. Return Original = ", num2str(annRet)));
display(strcat("Annual. Return Inter-Industry = ", num2str(annRet2)));
display(strcat("Annual. Return Prior-year = ", num2str(annRet3)));
display(strcat("Annual. Return Cross-Industry = ", num2str(annRet4)));
display(strcat("Annual. Return Principal Component = ", num2str(annRet5)));

display(strcat("Std Rotation Original = ", num2str(sqrt(12)*std(returns))));
display(strcat("Std Inter = ", num2str(sqrt(12)*std(returns2))));
display(strcat("Std prior year = ", num2str(sqrt(12)*std(returns3))));
display(strcat("Std cross = ", num2str(sqrt(12)*std(returns4))));
display(strcat("Std Principal Component = ", num2str(sqrt(12)*std(returns5))));


display(strcat("Max DD Original = ", num2str(round(maxdrawdown(cumprod(returns/100 + 1.0))*100,2))));
display(strcat("Max DD Inter indust = ", num2str(round(maxdrawdown(cumprod(returns2/100 + 1.0))*100,2))));
display(strcat("Max DD Prior year = ", num2str(round(maxdrawdown(cumprod(returns3/100 + 1.0))*100,2))));
display(strcat("Max DD Cross = ", num2str(round(maxdrawdown(cumprod(returns4/100 + 1.0))*100,2))));
display(strcat("Max DD Principal Component = ", num2str(round(maxdrawdown(cumprod(returns5/100 + 1.0))*100,2))));

figure;
plot(time,a,'k','LineWidth', 1.5 );hold on;
plot(time,a2,'r', 'LineWidth', 1.5);
plot(time,a3,'g', 'LineWidth', 1.5);
plot(time,a4,'b', 'LineWidth', 1.5);
plot(time,a5,'m', 'LineWidth', 1.5);

plot([min(xlim()),max(xlim())],[0,0], 'k--');

ax = gca; 
ylim([-0.31, 0.51]);
axis 'auto x'
datetick('x', 'yyyy', 'keepticks')
ax.FontSize = 12;
recessionplot;
%set(gca,'units','points','position',[10,10,400,250])
set(gcf,'units','points','position',[30,30,400,250])
%legend({'Original returns','Intra-Industry', 'Prior-year', 'Cross-Industry', 'Principal Component'},'Location','eastoutside','FontSize',14)
end


