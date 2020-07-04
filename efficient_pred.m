% This file contains a function for computing the fama french Sharpe ratio
% improvement predictions
function [SR_new, SR_oud, z] = efficient_pred(input)
SR_oud = sqrt(12)*mean(input)/std(input);
rho = autocorr(input);
rho1 = rho(2);
mu_p = mean(input);
SR = mean(input)/std(input);
zeta = (SR^2 + rho1^2)/(SR^2 + 1);
sigma_p = sqrt(mu_p^2 * ((1/zeta) - 1)); 
SR_new = (mu_p / sigma_p) * sqrt(12);
covar = (mu_p^2*(1-rho1^2))/(SR^2 + rho1^2);
theta =  (2*var(input)*sigma_p^2 - 2* std(input)*sigma_p*covar + ...
        0.5*mean(input)^2*var(input) + 0.5*mean(input)^2*sigma_p^2 - ...
        (mean(input)*mean(input)*covar^2)/(std(input)*sigma_p)) /( length(input));
z = (std(input)*mean(input) - sigma_p*mean(input))/ sqrt(theta);
end