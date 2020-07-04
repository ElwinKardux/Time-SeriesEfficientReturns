% Computes the efficient returns for the Principal Component model.
dataIndustries = csvread("12_indust_month_value.csv",1,1);
dataFF5 = csvread("FF5.CSV",4,1);
Rf_all = dataFF5(:,6);
% 445th observation is 1963:7
excInd = dataIndustries(445:end,:) - Rf_all;
[num_obs, num_indust] = size(excInd);

PC = zeros(size(excInd));
PC2 = zeros(size(excInd));
PC3 = zeros(size(excInd));
beta = zeros(5,12);
zeta = zeros(12,1);
SR = mean(excInd)./std(excInd);
mu_con = zeros(num_obs,num_indust);
correlatie = zeros(4,12);
%initialize
X_PC = ones(num_obs,num_indust);   eff_PC = ones(num_obs,num_indust);
SR_PC_new = zeros(num_indust,1);   z_PC = zeros(num_indust,1); 
sigmasq = zeros(12,1);
mu_p = mean(excInd);
for i = 1:12
    A = excInd;
    A(:,i) = [];
    coeff = pca(A);
    PC(:,i) = A*coeff(:,1);
    PC2(:,i) = A*coeff(:,2);
    PC3(:,i) = A*coeff(:,3);
    X = [ones(num_obs,1) lagmatrix(excInd(:,i),1)  lagmatrix(PC(:,i),1) lagmatrix(PC2(:,i),1) lagmatrix(PC3(:,i),1)];
    beta(:,i) = regress(excInd(2:end,i), X(2:end,:));
    rho1 = beta(2,i);
    PC_rho = beta(3,i);
    PC2_rho = beta(4,i);
    PC3_rho = beta(5,i);
    A = lagmatrix(excInd(:,i),1);
    temp = corrcoef(excInd(2:end,i), A(2:end,:));
    correlatie(1,i) = temp(2);
    A = lagmatrix(PC(:,i),1);
    temp = corrcoef(excInd(2:end,i), A(2:end,:));
    correlatie(2,i) = temp(2);
    A =  lagmatrix(PC2(:,i),1);
    temp = corrcoef(excInd(2:end,i), A(2:end,:));
    correlatie(3,i) = temp(2);
    A =  lagmatrix(PC3(:,i),1);
    temp =  corrcoef(excInd(2:end,i), A(2:end,:));
    correlatie(4,i) = temp(2);
    % Computes the Zeta as in the methodolgy is shown.
    PCLAG = lagmatrix(PC(:,i),1);
    temprho = corrcoef(excInd(2:end,i),PCLAG(2:end));
    sigmasq(i) = var(excInd(:,i))* (1-rho1^2 - PC_rho^2 -PC2_rho^2 -PC3_rho^2);
    zeta(i) = (SR(i)^2 + rho1^2 + (PC_rho^2*var(PC)/var(excInd(:,i))) + PC2_rho^2*var(PC2)/var(excInd(:,i))...
        + PC3_rho^2*var(PC3)/var(excInd(:,i)))./(SR(i)^2 + 1);
    for j = 2:num_obs
        lag_ind = excInd(j-1,i);
        lag_PC = PC(j-1,i);
        lag_PC2 = PC2(j-1,i);
        lag_PC3 = PC3(j-1,i);
        mu_con(j,i) = mu_p(i)*(1- rho1) + rho1*lag_ind +  PC_rho*lag_PC + PC2_rho*lag_PC2 + PC3_rho*lag_PC3;
        X_PC(j,i) = (mu_p(i) / zeta(i)) * (mu_con(j,i)/ (mu_con(j,i).^2 + sigmasq(i)));
        %Weight restrictions:
        if X_PC(j,i)>1
            X_PC(j,i) = 1;
        elseif X_PC(j,i) < 0
            X_PC(j,i) = 0;
        end
    end
    eff_PC(:,i) = X_PC(:,i).* excInd(:,i);
    SR_PC_new(i) = sqrt(12)*mean(eff_PC(:,i))/std(eff_PC(:,i));
    covar = cov(excInd(2:end,i),eff_PC(2:end,i));
    theta =  (2*var(excInd(:,i))*var(eff_PC(:,i)) - 2* std(excInd(:,i))*std(eff_PC(:,i))*covar(1,2) + ...
            0.5*mean(eff_PC(:,i))^2*var(excInd(:,i)) + 0.5*mean(excInd(:,i))^2*var(eff_PC(:,i)) - ...
            (mean(eff_PC(:,i))*mean(excInd(:,i))*covar(1,2)^2)/(std(excInd(:,i))*std(eff_PC(:,i)))) /( num_obs);
    z_PC(i) = (std(excInd(:,i))*mean(eff_PC(:,i)) - std(eff_PC(:,i))*mean(excInd(:,i)))/ sqrt(theta);
end
%% Computes the rotation portfolios for the PCA.
p_PC = 2*normcdf(-abs(z_PC));
[rot_PC_eff, SR_rot_PC_eff, ZindF, ZSR_F, Zmaxsum, Zminsum] = Rotation_AdaptiveLasso(eff_PC);
[rot_PC_MA, SR_rot_PC_MA] = Rotation_Momentum(eff_PC);
[rot_PC_PM, SR_rot_PC_PM] = Rotation_PrevailingMean(eff_PC);