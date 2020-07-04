%Function for computing the Sharpe ratio improvements with the original and
%the efficient returns as inputs.
function [SR,p] = SR_improve(input,eff)
%theta and z are given in equation (7) and (8).
covar = cov(input,eff);
theta =  (2*var(input)*var(eff) - 2* std(input)*std(eff)*covar(1,2) + ...
        0.5*mean(eff)^2*var(input) + 0.5*mean(input)^2*var(eff) - ...
        (mean(eff)*mean(input)*covar(1,2)^2)/(std(input)*std(eff))) /( length(input));
z = (std(input)*mean(eff) - std(eff)*mean(input))/ sqrt(theta);
SR = sqrt(12)*mean(eff)/std(eff);
p = 2*normcdf(-abs(z));
end
