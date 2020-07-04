% Computes the Fama french and industry returns prior year efficient
% returns. Also constructs the rotation portfolios.

%% Fama french 5 on 1963:7 - 2018:12 
dataFF5 = csvread("FF5.CSV",4,1);
ff5oud = dataFF5(1:end-15,1:5);
Rf = dataFF5(1:end-15,6);

SR_12FF5old = zeros(1,5);
SR_12FF5new = zeros(1,5);
ff5eff = zeros(size(ff5oud));
Xi12FF5 = zeros(size(ff5oud));
rho__FF5n = zeros(5,1);
z12FF5 = zeros(1,5);
for i=1:5
    [ff5eff(:,i),Xi12FF5(:,i), SR_12FF5new(i), SR_12FF5old(i), rho__FF5n(i), z12FF5(i)] = efficient2(ff5oud(:,i));
end
p12FF5 = 2*normcdf(-abs(z12FF5));
Rf_all = dataFF5(:,6);
clear dataFF5;

%% industry 1 lag on 1963:7- 2018:12
dataIndustries = csvread("12_indust_month_value.csv",1,1);
% 485th observation is 1963:7
excInd = dataIndustries(445:end,:) - Rf_all;
AnnMeanExcRet = 12 *mean(excInd);
AnnStdExcRet = sqrt(12)* std(excInd);
AnnSRExcInd = AnnMeanExcRet / AnnStdExcRet;
[num_obs, num_indust] = size(excInd);

SR12_Indold = zeros(12,1);
SR12_Indnew = zeros(12,1);
eff_ar12 = zeros(size(excInd));
X12ind = zeros(size(excInd));
rho__n = zeros(12,1);
z12_Ind = zeros(12,1);
for i=1:12
    [eff_ar12(:,i), X12ind(:,i), SR12_Indnew(i), SR12_Indold(i), rho__n(i), z12_Ind(i)] = efficient2(excInd(:,i));
end
p12_Ind = 2*normcdf(-abs(z12_Ind));

%% Rotation efficient portfolio
[rot_ar12_eff, SR_rot_ar12_eff] = Rotation_AdaptiveLasso(eff_ar12);
[rot_ar12_MA, SR_rot_ar12_MA] = Rotation_PrevailingMean(eff_ar12);
[rot_ar12_PM, SR_rot_ar12_PM] = Rotation_Momentum(eff_ar12);
