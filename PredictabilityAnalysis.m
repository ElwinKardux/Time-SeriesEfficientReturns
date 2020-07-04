%% Fama french 5 on 1963:7 - 2018:12 
dataFF5 = csvread("FF5.CSV",4,1);
ff5oud = dataFF5(1:end-15,1:5);
Rf = dataFF5(1:end-15,6);
Rf_all = dataFF5(:,6);

%% industry 1 lag on 1963:7- 2018:12
dataIndustries = csvread("12_indust_month_value.csv",1,1);
% 485th observation is 1963:7
%industry = dataIndustries(445:end,:);
excInd = dataIndustries(445:end,:) - Rf_all;


%% FF5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the Fama french five factor low and high values, with their
% t-values. 
ZlowHigh = zeros(6,5);
ZvarLowHigh = zeros(4,5);
% column 1: Low  
% column 2: High
% column 3: H-L
% column 4: Low t-stat
% column 5: High t-stat
% column 6: H-L t-stat
for i = 1:5
    [ZlowHigh(1,i), ZlowHigh(2,i), ZvarLowHigh(1,i), ZvarLowHigh(2,i), ZlowHigh(4,i),...
        ZlowHigh(5,i), ZlowHigh(6,i), ZvarLowHigh(4,i)] = panel(ff5oud(:,i));
    ZlowHigh(3,i) = ZlowHigh(2,i) - ZlowHigh(1,i);
    ZvarLowHigh(3,i) = (ZvarLowHigh(2,i)/ZvarLowHigh(1,i));
end

%% FF5  12Month: Computes the Fama French five factor low and high values using the past 12 months.
ZlowHigh12 = zeros(3,5);
ZvarLowHigh12 = zeros(4,5);
% column 1: Low  
% column 2: High
% column 3: H-L
% column 4: Low t-stat
% column 5: High t-stat
% column 6: H-L t-stat
for i = 1:5
    [ZlowHigh12(1,i), ZlowHigh12(2,i), ZvarLowHigh12(1,i), ZvarLowHigh12(2,i),...
        ZlowHigh12(4,i), ZlowHigh12(5,i), ZlowHigh12(6,i), ZvarLowHigh12(4,i)] = panel12(ff5oud(:,i));
    ZlowHigh12(3,i) = ZlowHigh12(2,i) - ZlowHigh12(1,i);
    ZvarLowHigh12(3,i) = (ZvarLowHigh12(2,i)/ZvarLowHigh12(1,i));
end

%% Industries!: Industry low and high using previous month.
lowHighInd = zeros(6,12);
varLowHighInd = zeros(4,12);
% column 1: Low  
% column 2: High
% column 3: H-L
% column 4: Low t-stat
% column 5: High t-stat
% column 6: H-L t-stat
for i = 1:12
    [lowHighInd(1,i), lowHighInd(2,i), varLowHighInd(1,i), varLowHighInd(2,i),...
        lowHighInd(4,i), lowHighInd(5,i), lowHighInd(6,i) , varLowHighInd(4,i)] = panel(excInd(:,i));
    lowHighInd(3,i) = lowHighInd(2,i) - lowHighInd(1,i);
    varLowHighInd(3,i) = (varLowHighInd(2,i)/varLowHighInd(1,i));
end
lowHighInd = lowHighInd';
varLowHighInd = varLowHighInd';
%% Industries 12 : Industry low and high using average of prior 12 months.
lowHigh12Ind = zeros(6,12);
varLowHigh12Ind = zeros(4,12);
% column 1: Low  
% column 2: High
% column 3: H-L
% column 4: Low t-stat
% column 5: High t-stat
% column 6: H-L t-stat
for i = 1:12
    [lowHigh12Ind(1,i), lowHigh12Ind(2,i), varLowHigh12Ind(1,i), varLowHigh12Ind(2,i), ...
        lowHigh12Ind(4,i), lowHigh12Ind(5,i), lowHigh12Ind(6,i), varLowHigh12Ind(4,i)] = panel12(excInd(:,i));
    lowHigh12Ind(3,i) = lowHigh12Ind(2,i) - lowHigh12Ind(1,i);
    varLowHigh12Ind(3,i) = (varLowHigh12Ind(2,i)/varLowHigh12Ind(1,i));
end
lowHigh12Ind = lowHigh12Ind';
varLowHigh12Ind = varLowHigh12Ind';
clear i

%% Average FF5: Computes the bottom `Combined' results for the Fama french factors.
ZAVlowHigh = zeros(6,1);
ZAVvarLowHigh = zeros(4,1);
% column 1: Low  
% column 2: High
% column 3: H-L
% column 4: Low t-stat
% column 5: High t-stat
% column 6: H-L t-stat
[ZAVlowHigh(1), ZAVlowHigh(2), ZAVvarLowHigh(1), ZAVvarLowHigh(2), ZAVlowHigh(4),...
    ZAVlowHigh(5), ZAVlowHigh(6), ZAVvarLowHigh(4)] = panelAv(ff5oud);
ZAVlowHigh(3) = ZAVlowHigh(2) - ZAVlowHigh(1);
ZAVvarLowHigh(3) = (ZAVvarLowHigh(2)./ZAVvarLowHigh(1));

ZAV12lowHigh = zeros(6,1);
ZAV12varLowHigh = zeros(4,1);
[ZAV12lowHigh(1), ZAV12lowHigh(2), ZAV12varLowHigh(1), ZAV12varLowHigh(2), ZAV12lowHigh(4),...
    ZAV12lowHigh(5), ZAV12lowHigh(6), ZAV12varLowHigh(4)] = panelAv12(ff5oud);
ZAV12lowHigh(3) = ZAV12lowHigh(2) - ZAV12lowHigh(1);
ZAV12varLowHigh(3) = (ZAV12varLowHigh(2)./ZAV12varLowHigh(1));

%% Average Industries: Computes the `combined' row results for the industry returns.
AVlowHighIND = zeros(6,1);
AVvarLowHighIND = zeros(4,1);
% column 1: Low  
% column 2: High
% column 3: H-L
% column 4: Low t-stat
% column 5: High t-stat
% column 6: H-L t-stat
[AVlowHighIND(1), AVlowHighIND(2), AVvarLowHighIND(1), AVvarLowHighIND(2), ...
    AVlowHighIND(4), AVlowHighIND(5), AVlowHighIND(6), AVvarLowHighIND(4)] = panelAv(excInd);
AVlowHighIND(3) = AVlowHighIND(2) - AVlowHighIND(1);
AVvarLowHighIND(3) = (AVvarLowHighIND(2)./AVvarLowHighIND(1));

AV12lowHighIND = zeros(6,1);
AV12varLowHighIND = zeros(4,1);
[AV12lowHighIND(1), AV12lowHighIND(2), AV12varLowHighIND(1), AV12varLowHighIND(2),...
    AV12lowHighIND(4), AV12lowHighIND(5), AV12lowHighIND(6),AV12varLowHighIND(4)] = panelAv12(excInd);
AV12lowHighIND(3) = AV12lowHighIND(2) - AV12lowHighIND(1);
AV12varLowHighIND(3) = (AV12varLowHighIND(2)./AV12varLowHighIND(1));

%% Cross Industry : Computes the low and high group results using the cross-industry model
[crossLow, crossHigh, crossVarLow, crossVarHigh, crossLowTstat, crossHighTstat,...
    crossHLTstat,crossVarFtest] = panelCross(excInd);

%% Principal component: Computes the low and high group results using the Principal Component model
[PCLow, PCHigh, PCVarLow, PCVarHigh, PCLowTstat, PCHighTstat,...
    PCHLTstat,PCVarFtest] = panelPC(excInd);
PCAVlowHighIND = zeros(6,1);
PCAVvarLowHighIND = zeros(4,1);
[PCAVlowHighIND(1), PCAVlowHighIND(2), PCAVvarLowHighIND(1), PCAVvarLowHighIND(2), ...
    PCAVlowHighIND(4), PCAVlowHighIND(5), PCAVlowHighIND(6), PCAVvarLowHighIND(4)] = panelPCAV(excInd);
%[a1, a2, a3, a4, a5 ,a6, a7, a8] 
PCAVlowHighIND(3) = PCAVlowHighIND(2) - PCAVlowHighIND(1);
PCAVvarLowHighIND(3) = (PCAVvarLowHighIND(2)./PCAVvarLowHighIND(1));

%% Functions:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1 lag Industry" Intra- industry model
function [low, high, var_low, var_high, tval_low, tval_high, tval_HL, var_ftest] = panel(factor)
lows = factor < median(factor);
highs = factor > median(factor);
count = 0;
sums = 0;
sums2= 0;
count2= 0;
vec_low = zeros(1,length(factor));
vec_high = zeros(1,length(factor));

for i = 2:length(factor)
    if lows(i-1) == 1
        sums = sums + factor(i);
        count = count +1;
        vec_low(i) = factor(i);
    end
    if highs(i-1) == 1
        sums2 = sums2 + factor(i);
        count2= count2 + 1;
        vec_high(i) = factor(i);
    end
end
low = sums/count;
high= sums2/count2;
vec_lows = vec_low(vec_low~=0);
vec_highs = vec_high(vec_high~=0);

var_low = var(vec_low)*(length(factor)-1)/(sum(lows)-1);
var_high = var(vec_high)*(length(factor)-1)/(sum(lows)-1);
tval_low = low/sqrt(var_low/count);
tval_high = high/sqrt(var_high/count2);
tval_HL =  (high - low)/sqrt(  (var_low/count) + (var_high/count2) );

[~,var_ftest] = vartest2(vec_highs,vec_lows);
end


%% 12 lagged returns!
% Prior 12 months low and high function!
function [low, high, var_low, var_high, tval_low, tval_high, tval_HL, var_ftest] = panel12(factor)
av12lag = zeros(1,length(factor));
for i = 2:length(factor)
    if i < 13
    %if i > 12
       av12lag(i) = mean(factor(1:i-1));
    else
        av12lag(i) = sum(factor(i-12:i-1))/12;
    end
end
lows = av12lag <= median(av12lag);
highs = av12lag >= median(av12lag);
count =0;
sums = 0;
sums2=0;
count2=0;
vec_low = zeros(1,length(factor));
vec_high = zeros(1,length(factor));
for i = 1:length(factor)
    if lows(i) == 1
        sums = sums + factor(i);
        count = count +1;
        vec_low(i) = factor(i);
    end
    if highs(i) == 1
        sums2 = sums2 + factor(i);
        count2= count2 + 1;
        vec_high(i) = factor(i);
    end
end
low = sums/count;
high= sums2/count2;
var_low = var(vec_low)*(length(factor)-1)/(sum(lows)-1);
var_high = var(vec_high)*(length(factor)-1)/(sum(lows)-1);

tval_low = low/sqrt(var_low/count);
tval_high = high/sqrt(var_high/count2);
tval_HL =  (high - low)/sqrt(  (var_low/count) + (var_high/count2) );
vec_lows = vec_low(vec_low~=0);
vec_highs = vec_high(vec_high~=0);
[~,var_ftest] = vartest2(vec_highs,vec_lows);
end


%% Average
%Average” at the bottom of the table is computed by
%assigning each month for a factor into high or low a group and then computing the average returns
%across all factors. The high or low groups in month t can be empty if no factor gets assigned into
%it.
function [low, high, var_low, var_high, tval_low, tval_high, tval_HL, var_ftest] = panelAv(data)
[num_obs, num_factors] = size(data);
lows = data < median(data);
highs = data > median(data);
sums = 0;
sums2= 0;
noneLow = 0;
noneHigh = 0;
vec_low = zeros(1,num_obs);
vec_high = zeros(1, num_obs);
for i = 2:num_obs
    sums = sums + sum(data(i,:).*lows(i-1,:))./(sum(lows(i-1,:))+(10^(-6)));
    sums2 = sums2 + sum(data(i,:).*highs(i-1,:))./(sum(highs(i-1,:))+(10^(-6)));
    vec_low(i) = sum(data(i,:).*lows(i-1,:))./(sum(lows(i-1,:))+(10^(-6)));
    vec_high(i) = sum(data(i,:).*highs(i-1,:))./(sum(highs(i-1,:))+(10^(-6)));
    % 10^-6 to avoid division by zero
    if sum(lows(i-1,:))==0
        noneLow = noneLow + 1;
    end
    if sum(highs(i-1,:))==0
        noneHigh = noneHigh + 1;
    end
end
low = sums/(num_obs-noneLow);
high= sums2/(num_obs-noneHigh);
var_low = var(vec_low);
var_high = var(vec_high);
tval_low = low./sqrt(var_low./(num_obs-noneLow));
tval_high = high./sqrt(var_high./(num_obs-noneHigh));
tval_HL =  (high - low)./sqrt(  (var_low./(num_obs-noneLow)) + (var_high./(num_obs-noneHigh)) );

vec_lows = vec_low(vec_low~=0);
vec_highs = vec_high(vec_high~=0);
[~,var_ftest] = vartest2(vec_highs,vec_lows);
end


%% Average 12 Lags
%Average” at the bottom of the table is computed by
%assigning each month for a factor into high or low a group and then computing the average returns
%across all factors. The high or low groups in month t can be empty if no factor gets assigned into
%it.
function [low, high, var_low, var_high, tval_low, tval_high, tval_HL, var_ftest] = panelAv12(data)
[num_obs, num_factors] = size(data);
av12lag = zeros(num_obs, num_factors);
for i = 2:num_obs
    if i < 13
    %if i > 12
       av12lag(i,:) = mean(data(1:i-1,:));
    else
        av12lag(i,:) = sum(data(i-12:i-1,:))./12;
    end
end

lows = av12lag < median(av12lag);
highs = av12lag > median(av12lag);
sums = 0;
sums2= 0;
noneLow = 0;
noneHigh = 0;
vec_low = zeros(1,num_obs);
vec_high = zeros(1, num_obs);
for i = 13:num_obs
    sums = sums + sum(data(i,:).*lows(i,:))./(sum(lows(i,:))+(10^(-6)));
    sums2 = sums2 + sum(data(i,:).*highs(i,:))./(sum(highs(i,:))+(10^(-6)));
    vec_low(i) = sum(data(i,:).*lows(i,:))./(sum(lows(i,:))+(10^(-6)));
    vec_high(i) = sum(data(i,:).*highs(i,:))./(sum(highs(i,:))+(10^(-6)));
    if sum(lows(i,:))==0
        noneLow = noneLow + 1;
    end
    if sum(highs(i,:))==0
        noneHigh = noneHigh + 1;
    end
end
low = sums/(num_obs-noneLow);
high= sums2/(num_obs-noneHigh);
var_low = var(vec_low);
var_high = var(vec_high);
tval_low = low./sqrt(var_low./(num_obs-noneLow));
tval_high = high./sqrt(var_high./(num_obs-noneHigh));
tval_HL =  (high - low)./sqrt(  (var_low./(num_obs-noneLow)) + (var_high./(num_obs-noneHigh)) );
vec_lows = vec_low(vec_low~=0);
vec_highs = vec_high(vec_high~=0);
[~,var_ftest] = vartest2(vec_highs,vec_lows);
end



%% Cross-Industry lagged returns!
% Predictability function for the Cross-Industry model.
function [low, high, var_lows, var_highs, tval_low, tval_high, tval_HL, var_ftest] = panelCross(data)
[num_obs, num_indust] = size(data);
XLAG = lagmatrix(data,1);
forec = zeros(num_obs-1,num_indust);
low = zeros(num_indust,1);
high = zeros(num_indust,1);
var_lows = zeros(num_indust,1);
var_highs = zeros(num_indust,1);
tval_low = zeros(num_indust,1);
tval_high = zeros(num_indust,1);
tval_HL = zeros(num_indust,1);
var_ftest = zeros(num_indust,1);
for i = 1:num_indust
    [adaptive, ~] = adapLasso(XLAG(2:end,:), data(2:end,i), 0);
    intercept = adaptive(1);
    lasBeta = adaptive(2:13);
    for j = 2:num_obs
         forec(j,i) = intercept +  XLAG(j,:)*lasBeta;
    end
    lows = forec(:,i) <= median(forec(:,i));
    highs = forec(:,i) >= median(forec(:,i));
    count =0;
    sums = 0;
    sums2=0;
    count2=0;
    vec_low = zeros(num_obs-1, 1);
    vec_high = zeros(num_obs-1, 1);
    for j = 2:num_obs
        if lows(j) == 1
            sums = sums + data(j,i);
            count = count +1;
            vec_low(j) = data(j,i);
        end
        if highs(j) == 1
            sums2 = sums2 + data(j);
            count2= count2 + 1;
            vec_high(j) = data(j);
        end
    end
    low(i) = sums/count;
    high(i) = sums2/count2;
    
    vec_lows = vec_low(vec_low~=0);
    vec_highs = vec_high(vec_high~=0);
    var_lows(i) = var(vec_lows);
    var_highs(i) = var(vec_highs);
    tval_low(i) = low(i)/sqrt(var_lows(i)/count);
    tval_high(i) = high(i)/sqrt(var_highs(i)/count2);
    tval_HL(i) =  (high(i) - low(i))/sqrt(  (var_lows(i)/count) + (var_highs(i)/count2) );
    [~,var_ftest(i)] = vartest2(vec_highs,vec_lows);
end
end

%% Principal components!
% Preditability function for the Principal Component model.
function [low, high, var_lows, var_highs, tval_low, tval_high, tval_HL, var_ftest] = panelPC(data)
[num_obs, num_indust] = size(data);

%Initialisatie
PC = zeros(size(data)); PC2 = zeros(size(data)); PC3 = zeros(size(data));
beta = zeros(5,12);  forec = zeros(num_obs-1,num_indust);
low = zeros(num_indust,1); high = zeros(num_indust,1);
var_lows = zeros(num_indust,1); var_highs = zeros(num_indust,1);
tval_low = zeros(num_indust,1); tval_high = zeros(num_indust,1); 
tval_HL = zeros(num_indust,1); var_ftest = zeros(num_indust,1);
for i = 1:num_indust
    for j = 3:num_obs
        A = data;
        A(:,i) = [];
        coeff = pca(A(1:j,:));
        PC(:,i) = A*coeff(:,1);
        PC2(:,i) = A*coeff(:,2);
        PC3(:,i) = A*coeff(:,3);
        X = [ones(num_obs,1) lagmatrix(data(:,i),1)  lagmatrix(PC(:,i),1) lagmatrix(PC2(:,i),1) lagmatrix(PC3(:,i),1)];
        beta(:,i) = regress(data(2:end,i), X(2:end,:));
        rho1 = beta(2,i);
        PC_rho = beta(3,i);
        PC2_rho = beta(4,i);
        PC3_rho = beta(5,i);
        
        forec(j,i) = beta(1,i) + rho1*data(j-1,i) + PC_rho*PC(j-1,i) +PC2_rho*PC2(j-1,i)+PC3_rho*PC3(j-1,i);
    end
    lows = forec(:,i) <= median(forec(:,i));
    highs = forec(:,i) >= median(forec(:,i));
    count =0;
    sums = 0;
    sums2=0;
    count2=0;
    vec_low = zeros(num_obs-1, 1);
    vec_high = zeros(num_obs-1, 1);
    for j = 2:num_obs
        if lows(j) == 1
            sums = sums + data(j,i);
            count = count +1;
            vec_low(j) = data(j,i);
        end
        if highs(j) == 1
            sums2 = sums2 + data(j);
            count2= count2 + 1;
            vec_high(j) = data(j);
        end
    end
    low(i) = sums/count;
    high(i) = sums2/count2;
    
    vec_lows = vec_low(vec_low~=0);
    vec_highs = vec_high(vec_high~=0);
    var_lows(i) = var(vec_lows);
    var_highs(i) = var(vec_highs);
    tval_low(i) = low(i)/sqrt(var_lows(i)/count);
    tval_high(i) = high(i)/sqrt(var_highs(i)/count2);
    tval_HL(i) =  (high(i) - low(i))/sqrt(  (var_lows(i)/count) + (var_highs(i)/count2) );
    [~,var_ftest(i)] = vartest2(vec_highs,vec_lows);
end
end


%% Principal Component model Average results!
function [low, high, var_lows, var_highs, tval_low, tval_high, tval_HL, var_ftest] = panelPCAV(data)
[num_obs, num_indust] = size(data);
%Initialisation
PC = zeros(size(data)); PC2 = zeros(size(data)); PC3 = zeros(size(data));
beta = zeros(5,12);  forec = zeros(num_obs-1,num_indust);
for i = 1:num_indust
    %mdl = glm_logistic(excInd(:,i),XLAG, "nointercept");
    A = data;
    A(:,i) = [];
    coeff = pca(A);
    PC(:,i) = A*coeff(:,1);
    PC2(:,i) = A*coeff(:,2);
    PC3(:,i) = A*coeff(:,3);
    X = [ones(num_obs,1) lagmatrix(data(:,i),1)  lagmatrix(PC(:,i),1) lagmatrix(PC2(:,i),1) lagmatrix(PC3(:,i),1)];
    beta(:,i) = regress(data(2:end,i), X(2:end,:));
    rho1 = beta(2,i);
    PC_rho = beta(3,i);
    PC2_rho = beta(4,i);
    PC3_rho = beta(5,i);
    for j = 2:num_obs
         forec(j,i) = beta(1,i) + rho1*data(j-1,i) + PC_rho*PC(j-1,i) +PC2_rho*PC2(j-1,i)+PC3_rho*PC3(j-1,i);
    end
end

lows = forec <= median(forec);
highs = forec >= median(forec);
count =0;
sums = 0;
sums2=0;
count2=0;
vec_low = zeros(num_obs, 1);
vec_high = zeros(num_obs, 1);
noneLow = 0;
noneHigh = 0;
    
for i = 2:num_obs
    sums = sums + sum(data(i,:).*lows(i,:))./(sum(lows(i,:))+(10^(-6)));
    sums2 = sums2 + sum(data(i,:).*highs(i,:))./(sum(highs(i,:))+(10^(-6)));
    vec_low(i) = sum(data(i,:).*lows(i,:))./(sum(lows(i,:))+(10^(-6)));
    vec_high(i) = sum(data(i,:).*highs(i,:))./(sum(highs(i,:))+(10^(-6)));
    if sum(lows(i,:))==0
        noneLow = noneLow + 1;
    end
    if sum(highs(i,:))==0
        noneHigh = noneHigh + 1;
    end
end
low = sums/(num_obs-noneLow);
high= sums2/(num_obs-noneHigh);
vec_lows = vec_low(vec_low~=0);
vec_highs = vec_high(vec_high~=0);
var_lows = var(vec_lows);
var_highs = var(vec_highs);
tval_low = low./sqrt(var_lows./(num_obs-noneLow));
tval_high = high./sqrt(var_highs./(num_obs-noneHigh));
tval_HL =  (high - low)./sqrt(  (var_lows./(num_obs-noneLow)) + (var_highs./(num_obs-noneHigh)) );

[~,var_ftest] = vartest2(vec_highs,vec_lows);
end