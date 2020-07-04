% Computes efficient industry returns using its inter industry lag
% Computes Rotation portfolio based on these efficient returns
% Also computes the predictions for the sharpe ratio improvements.
%% Fama french 5 on 1963:7 - 2018:12 
dataFF5 = csvread("FF5.CSV",4,1);
ff5oud = dataFF5(1:end-15,1:5);
Rf = dataFF5(1:end-15,6);

SR_oldFF5 = zeros(1,5);
SR_newFF5 = zeros(1,5);
ff5eff = zeros(size(ff5oud));
XiFF5 = zeros(size(ff5oud));
zFF5 = zeros(1,5);
for i=1:5
    [ff5eff(:,i),XiFF5(:,i), SR_newFF5(i), SR_oldFF5(i), zFF5(i)] = efficient(ff5oud(:,i));
end
pFF5 = 2*normcdf(-abs(zFF5));
Rf_all = dataFF5(:,6);

SR_oldFF5pred = zeros(1,5);
SR_newFF5pred = zeros(1,5);
z_FF5pred = zeros(1,5);
for i=1:5
    [SR_newFF5pred(i), SR_oldFF5pred(i), z_FF5pred(i)] = efficient_pred(ff5oud(:,i));
end
p_FF5pred = 2*normcdf(-abs(z_FF5pred));
clear dataFF5;

%% industry 1 lag on 1963:7- 2020:03
dataIndustries = csvread("12_indust_month_value.csv",1,1);
excInd = dataIndustries(445:end,:) - Rf_all;
[num_obs, num_indust] = size(excInd);

SR_oldIND = zeros(1,12);
SR_newIND = zeros(1,12);
eff_ar1 = zeros(size(excInd));
Xind = zeros(size(excInd));
zInd = zeros(12,1);
for i=1:12
    [eff_ar1(:,i), Xind(:,i), SR_newIND(i), SR_oldIND(i), zInd(i)] = efficient(excInd(:,i));
end
pInd = 2*normcdf(-abs(zInd));

%% Rotation efficient portfolio
[rot_old, SR_rot_old] = Rotation_AdaptiveLasso(excInd);
[rot_ar1_eff, SR_rot_ar1_eff] = Rotation_AdaptiveLasso(eff_ar1);

[rot_old_MA, SR_rot_old_MA] = Rotation_PrevailingMean(excInd);
[rot_eff_MA, SR_rot_eff_MA] = Rotation_PrevailingMean(eff_ar1);
[rot_old_PM, SR_rot_old_PM] = Rotation_Momentum(excInd);
[rot_ar1_PM, SR_rot_ar1_PM] = Rotation_Momentum(eff_ar1);

%% Prediction Industries 1 lag
SR_oldInd_pred = zeros(1,12);
SR_Ind_pred = zeros(1,12);
z_Ind_pred = zeros(1,12);
for i=1:12
    [SR_Ind_pred(i), SR_oldInd_pred(i), z_Ind_pred(i)] = efficient_pred(excInd(:,i));
end
p_Ind_pred = 2*normcdf(-abs(z_Ind_pred));