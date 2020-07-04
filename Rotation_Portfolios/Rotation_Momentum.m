% Rotation portfolio using Cross-Sectional momentum!
function [Rotat, SR] = Rotation_Momentum(returns)
start_out = 259;
% 258 = January 1985;
[num_obs, num_indust] = size(returns);
MA = zeros(num_obs, num_indust);
%This portfolio forecasts using the cumulative product of the past 12
%months. 
for j = 1:num_indust
    for i = start_out:num_obs
        temp = cumprod(  returns(i-13:i-1,j)/100+1  );
        MA(i,j) = temp(end);
    end
end
maxsum = zeros(1,num_indust);
minsum = zeros(1,num_indust);

portf = zeros(num_obs, 4);
%This for loop computes the 4 necessary returns for the rotation portfolio.
%portf 1 amd 2 are the 2 largest predicted returns, and 3 and 4 the lowest.
for i = start_out: num_obs
    cur = MA(i,:);
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

Rotat = 0.5 *   portf(:,1) + 0.5 * portf(:,2) - 0.5 * portf(:,3) - 0.5 * portf(:,4);
Rotat = Rotat(start_out:end);
SR = sqrt(12)*mean(Rotat)/std(Rotat);


end