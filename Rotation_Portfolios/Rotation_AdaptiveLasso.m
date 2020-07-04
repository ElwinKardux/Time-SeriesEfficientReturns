% Rotation portfolio using the predictive Adaptive Lasso method.
function [Rotat, SR, indF, SR_F, maxsum, minsum] = Rotation_AdaptiveLasso(returns)
start_out = 259;
% observation 259 = January 1985;
[num_obs, num_indust] = size(returns);
XLAG = lagmatrix(returns,1);

forec = zeros(size(returns));
%These lambdas are chosen from the whole sample cross validation. 
%Using cross validation on the moving window costs too much time.
lambdas = [0.638803536093476;1.65511965219752;0.634787510178556;1.42207031250000;...
    1.66464843750000;2.63320312500000;1.36582031250000;1.17597656250000;...
    0.856982740539195;0.684644096098084;1.26322443286443;0.579732438011590];
for i = 1:12
    for j = start_out: num_obs
        [adaptive, ~] = adapLasso(XLAG(j-start_out +1: j-1 ,:), returns(j-start_out +1: j-1 ,i), lambdas(i));
        intercept = adaptive(1);
        lasBeta = adaptive(2:13);
        forec(j,i) = intercept +  XLAG(j,:)*lasBeta;
    end
end
indF = forec(start_out:end,:);
SR_F = sqrt(12) .* mean(forec(start_out+1:end,:))./std(forec(start_out+1:end,:));
maxsum = zeros(1,num_indust);
minsum = zeros(1,num_indust);
portf = zeros(num_obs, 4);
%This for loop computes the 4 necessary returns for the rotation portfolio.
%portf 1 amd 2 are the 2 largest predicted returns, and 3 and 4 the lowest.
for i = start_out:num_obs
    cur = forec(i,:);
    [~, Imax] = maxk(cur,2);
    [~, Imin] = mink(cur,2);
    maxsum(Imax(1)) = maxsum(Imax(1)) + 1;
    maxsum(Imax(2)) = maxsum(Imax(2)) + 1;
    minsum(Imin(1)) = minsum(Imin(1)) + 1;
    minsum(Imin(2)) = minsum(Imin(2)) + 1;
    portf(i,1) = returns(i, Imax(1));
    portf(i,2) = returns(i, Imax(2));
    portf(i,3) = returns(i, Imin(1));
    portf(i,4) = returns(i, Imin(2));
end
portf = portf(start_out:num_obs,:);
Rotat = 0.5 * portf(:,1) + 0.5 * portf(:,2) - 0.5 * portf(:,3) - 0.5 * portf(:,4);
SR = sqrt(12)*mean(Rotat(2:end))/std(Rotat(2:end));

end
