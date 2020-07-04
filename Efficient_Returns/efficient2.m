function [eff,X, SR_new, SR_oud, rho_n, z] = efficient2(input)
% Computes the Prior year efficient returns and Sharpe ratio improvements
av12lag = zeros(1,length(input));
for i = 2:length(input)
    %if i < 13
    if i > 12
    %    av12lag(i) = mean(input(1:i-1));
    %else
        av12lag(i) = sum(input(i-12:i-1))/12;
    end
end
SR_oud = sqrt(12)*mean(input)/std(input);
rhos = corrcoef(av12lag, input);
rho_n = rhos(1,2);
rho = autocorr(input);
rho1 = rho(2);
%betasumsq = sum(beta.^2);
mu_p = mean(input);
SR = mean(input)/std(input);
zeta = (SR^2 + rho_n^2)/(SR^2 + 1);
%zeta = (SR^2 + betasumsq)/(SR^2 + 1);
mu_con = zeros(1,length(input));
sigma = (std(input))^2* (1-rho1^2);
X = ones(length(input),1);
b_n = regress(input, av12lag');
% Computes the Weights X for each observation. This is done using the
% formulas in the Methodology.
for i = 13:length(input)
    mu_con(i) = mu_p*(1-rho1) + b_n*av12lag(i);
    X(i) = (mu_p ./ zeta)* (mu_con(i) / (mu_con(i)^2 + sigma));
    % Weight restrictions
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