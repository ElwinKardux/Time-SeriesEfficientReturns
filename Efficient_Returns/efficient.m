% Function for computing the Intra industry efficient weights and Sharpe
% ratio improvements
function [eff,X, SR_new, SR_oud, z] = efficient(input)
SR_oud = sqrt(12)*mean(input)/std(input);
rho = autocorr(input);
rho1 = rho(2);
mu_p = mean(input);
SR = mean(input)/std(input);
zeta = (SR^2 + rho1^2)/(SR^2 + 1);
mu_con = zeros(length(input),1);
sigma = (std(input))^2* (1-rho1^2);
X = ones(length(input),1);
% Computes the Weights X for each observation. This is done using the
% formulas in the Methodology.
for i = 2:length(input)   
    return_lag = input(i-1);
    mu_con(i) = mu_p*(1-rho1) + rho1*return_lag;
    X(i) = (mu_p / zeta)* (mu_con(i) / (mu_con(i)^2 + sigma));
    %Weight restrictions
    if X(i)>1
        X(i) = 1;
    elseif X(i) < 0
        X(i) = 0;
    end
end
eff = X .* input;
SR_new = sqrt(12)*mean(eff)/std(eff);
covar = cov(input,eff);
theta =  (2*var(input)*var(eff) - 2* std(input)*std(eff)*covar(1,2) + ...
        0.5*mean(eff)^2*var(input) + 0.5*mean(input)^2*var(eff) - ...
        (mean(eff)*mean(input)*covar(1,2)^2)/(std(input)*std(eff))) /( length(input));
z = (std(input)*mean(eff) - std(eff)*mean(input))/ sqrt(theta);
end