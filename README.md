# Time-Series Efficient Industry Returns
Exploiting Inter- and Intradependencies in Industry Portfolios
A short description for each file:

* **PredictabilityAnalysis.m:**
Creates the predictability results for the Industry returns and for the Fama-French five factor returns

* **efficient.m:**
Function for creating the Intra-Industry efficient returns. It also computes its Sharpe Ratio improvements.

* **efficient2.m:**
Function for creating the Prior year Industry efficient returns. It also computes its Sharpe Ratio improvements.

* **Cross\_Industry.m:**
Creates the Cross-Industry efficient returns and computes its Sharpe Ratio improvements.

* **PCA.m:**
Creates the Principal component efficient returns and computes its Sharpe Ratio improvements.

* **SR\_improve.m:**
Function for computing the Sharpe ratio improvement z-value and p-value.

* **Intra\_efficient.m:**
Uses the efficient function of the efficient.m file for creating the Intra-Industry efficient returns for the Fama-french five factors and the twelve industry returns. Also constructs its rotation portfolios and computes the predictions for the Sharpe ratio improvements using the function in efficient\_pred.m. These predictions are reported in this research only for the Fama-french five factors.

* **efficient\_pred.m:**
This file contains a function for computing the efficient Fama French five factor Sharpe ratio improvement predictions.

* **PriorYear\_efficient.m:**
Uses the efficient2 function of the efficient2.m file for creating the Prior year Industry efficient returns for the Fama-french five factors and the twelve industry returns.

* **Rotation\_PrevailingMean.m**
Constructs industry rotation portfolios using the Prevailing mean forecasts.

* **Rotation\_AdaptiveLasso.m**
Constructs industry rotation portfolios using the predictive Adaptive Lasso forecasts.

* **Rotation\_Momentum.m**
Constructs industry rotation portfolios using the Cross-Sectional Momentum forecasts.

* **rotation\_main.m:**
This file computes all the necessary results for the rotation portfolios. The file is divided in three parts: Alpha Analysis, Recession Analysis and plotting the figures.

* **equal\_weighted\_analyses.m:**
This file computes the results for the equally weighted portfolios in the original and the efficient returns. The analyses are the Alpha Analysis and the Recession Analysis.

* **logcumplot.m:**
Creates the figures for the equally weighted portfolios in the original and the efficient industry returns.

* **adapLasso.m**
Necessary function for computing the adaptive lasso parameters for the cross-industry model.

If you want to run the **Cross\_Industry.m** file, you have to download the penalized toolbox package first. The download link and the details on this package can be found in https://www.jstatsoft.org/article/view/v072i06.

McIlhagga, W. (2016). penalized: A MATLAB Toolbox for Fitting Generalized Linear Models with Penalties. Journal of Statistical Software, 72(6), 1 - 21. doi:http://dx.doi.org/10.18637/jss.v072.i06

