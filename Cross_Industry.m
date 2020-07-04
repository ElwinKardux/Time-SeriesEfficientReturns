% Computes the Cross industry efficient returns and its Sharpe ratio
% improvements. Also computes the rotation portfolios for the
% cross_industry using the rotation files.
dataIndustries = csvread("12_indust_month_value.csv",1,1);
dataFF5 = csvread("FF5.CSV",4,1);
Rf_all = dataFF5(:,6);
% 485th observation is 1963:7
excInd = dataIndustries(445:end,:) - Rf_all;
[num_obs, num_indust] = size(excInd);
XLAG = lagmatrix(excInd,1);
adaptive= zeros(13,12);
lambdas = zeros(12,1);
intercept = zeros(12,1);
lasBeta = zeros(12,12);
% each column is an industry.
for i = 1:12
    [adaptive(:,i), lambdas(i)] = adapLasso(XLAG(2:end,:), excInd(2:end,i), 0);
    intercept(i) = adaptive(1,i);
    lasBeta(:,i) = adaptive(2:13,i);
end
SR_oud = (sqrt(12).*mean(excInd)./std(excInd))';
%rho = autocorr(excInd);
%rho1 = rho(2);
betasum = sum(lasBeta);
betasumsq = sum(lasBeta.^2);
mu_p = mean(excInd);
SR = mean(excInd)./std(excInd);
zeta = (SR.^2 + betasumsq)./(SR.^2 + 1);
mu_con = zeros(num_obs,num_indust);
sigmasq = var(excInd).* (1-betasumsq);
%initialize
X_cross = ones(num_obs,num_indust);   eff_cross = ones(num_obs,num_indust);
SR_cross_new = zeros(num_indust,1);   z_cross = zeros(num_indust,1); 
%Computes the Cross industry weights X for each observation and industry. 
for j = 1:num_indust
    for i = 2:num_obs   
        mu_con(i,j) = mu_p(j).*(1-lasBeta(j,j)) + XLAG(i,:)*lasBeta(:,j);
        X_cross(i,j) = (mu_p(j) ./ zeta(j)) .* (mu_con(i,j) ./ (mu_con(i,j).^2 + sigmasq(j)));
        if X_cross(i,j)>1
            X_cross(i,j) = 1;
        elseif X_cross(i,j) < 0
            X_cross(i,j) = 0;
        end
    end
    eff_cross(:,j) = X_cross(:,j).* excInd(:,j);
    SR_cross_new(j) = sqrt(12)*mean(eff_cross(:,j))/std(eff_cross(:,j));
    covar = cov(excInd(:,j),eff_cross(:,j));
    theta =  (2*var(excInd(:,j))*var(eff_cross(:,j)) - 2* std(excInd(:,j))*std(eff_cross(:,j))*covar(1,2) + ...
            0.5*mean(eff_cross(:,j))^2*var(excInd(:,j)) + 0.5*mean(excInd(:,j))^2*var(eff_cross(:,j)) - ...
            (mean(eff_cross(:,j))*mean(excInd(:,j))*covar(1,2)^2)/(std(excInd(:,j))*std(eff_cross(:,j)))) /( num_obs);
    z_cross(j) = (std(excInd(:,j))*mean(eff_cross(:,j)) - std(eff_cross(:,j))*mean(excInd(:,j)))/ sqrt(theta);
end
p_cross = 2*normcdf(-abs(z_cross));
%% Rotation efficient portfolio
[rot_cross_eff, SR_rot_cross_eff, ZindF, ZSR_F, Zmaxsum, Zminsum] = Rotation_AdaptiveLasso(eff_cross);
[rot_cross_MA, SR_rot_cross_MA] = Rotation_Momentum(eff_cross);
[rot_cross_PM, SR_rot_cross_PM] = Rotation_PrevailingMean(eff_cross);
