%Computes the equally weighted portfolio analyses results for the Alpha
%analysis and the Recession analysis.

Alleff = load('all_efficient.mat');
eff_inter = mean((Alleff.Alleff_ar1)');
eff_prior = mean((Alleff.Alleff_ar12)');
eff_cross = mean((Alleff.Alleff_cross)');
original = mean((Alleff.excInd)');
PCA = load('rot_PCA.mat');
eff_PC = mean((PCA.eff_PC)');

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
ff5effX = ff5eff;
%% OLS regression
% ALL contains all the information necessary for the equally weighted
% portfolio Alpha analyses
ALL = zeros(10,7);

model_old_PM= fitlm(ff5effX, original);
ALL(1:2,1:6) = table2array(model_old_PM.Coefficients(:,1:2))';
ALL(1,7) = model_old_PM.Rsquared.Ordinary;

model_ar1_PM= fitlm(ff5effX, eff_inter);
ALL(3:4,1:6) = table2array(model_ar1_PM.Coefficients(:,1:2))';
ALL(3,7) = model_ar1_PM.Rsquared.Ordinary;

model_ar12_PM= fitlm(ff5effX, eff_prior);
ALL(5:6,1:6) = table2array(model_ar12_PM.Coefficients(:,1:2))';
ALL(5,7) = model_ar12_PM.Rsquared.Ordinary;

model_cross_PM= fitlm(ff5effX, eff_cross);
ALL(7:8,1:6) = table2array(model_cross_PM.Coefficients(:,1:2))';
ALL(7,7) = model_cross_PM.Rsquared.Ordinary;

model_cross_PM= fitlm(ff5effX, eff_PC);
ALL(9:10,1:6) = table2array(model_cross_PM.Coefficients(:,1:2))';
ALL(9,7) = model_cross_PM.Rsquared.Ordinary;

ALL(:,1) = ALL(:,1)*12;


%% Recession: 
VIX = csvread("USREC.csv",1,1);
V = VIX;
NAI = csvread("NAI.csv",1,1);
NA = quantile(NAI,0.20);
N = NAI <= NA;
B = zeros(10,2);
%B contains all the necessary information for the Recession analyeses for
%the equally weighted portfolio.
%%
model_old_PM= fitlm(V, original(end-422:end));
B(9:10,1) = table2array(model_old_PM.Coefficients(2,1:2))';
model_old_PM= fitlm(N, original(end-422:end));
B(9:10,2) = table2array(model_old_PM.Coefficients(2,1:2))';

%%
model_old_PM= fitlm(V, eff_inter(end-422:end));
B(1:2,1) = table2array(model_old_PM.Coefficients(2,1:2))';
model_old_PM= fitlm(N, eff_inter(end-422:end));
B(1:2,2) = table2array(model_old_PM.Coefficients(2,1:2))';

%% 
model_ar1_PM= fitlm(V, eff_prior(end-422:end));
B(3:4,1) = table2array(model_ar1_PM.Coefficients(2,1:2))';
model_ar1_PM= fitlm(N, eff_prior(end-422:end));
B(3:4,2) = table2array(model_ar1_PM.Coefficients(2,1:2))';
%%
model_ar12_PM= fitlm(V, eff_cross(end-422:end));
B(5:6,1) = table2array(model_ar12_PM.Coefficients(2,1:2))';
model_ar12_PM= fitlm(N, eff_cross(end-422:end));
B(5:6,2) = table2array(model_ar12_PM.Coefficients(2,1:2))';
%%
model_cross_PM= fitlm(V, eff_PC(end-422:end));
B(7:8,1) = table2array(model_cross_PM.Coefficients(2,1:2))';
model_cross_PM= fitlm(N, eff_PC(end-422:end));
B(7:8,2) = table2array(model_cross_PM.Coefficients(2,1:2))';