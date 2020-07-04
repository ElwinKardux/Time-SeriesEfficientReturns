% Necessary function for the Cross-Industry model. This function uses the
% penalized package.
function [beta, lambdas] = adapLasso(X2,y2, minLambda)
ols_model= fitlm(X2,y2);
ols_coeff = ols_model.Coefficients;
beta_ols = table2array(ols_coeff(:,1));
model2 = glm_gaussian(y2,X2);
if minLambda == 0
    cv2 = cv_penalized(model2, @p_lasso, "folds", 5);
    lamInd = cv2.minlambda;
%usually the minimum lambda minus one standard deviation gives a better
%solution, hence I multiply it by 0.90
    minLam = cv2.lambda(lamInd)*0.90;
else
    minLam = minLambda;
end
%for more details on this package, please visit:
%https://www.jstatsoft.org/article/view/v072i06
fit4 = penalized(model2, @p_adaptive, "gamma",  0.01, "lambda", minLam, "adaptivewt", {beta_ols});
beta = fit4.beta;
lambdas = fit4.lambda;
end