function [DM,p_value] = dmtest_modified(e1, e2,dfun, h)
%DMTEST: Retrieves the Diebold-Mariano test statistic (1995) for the 
% equality of forecast accuracy of two forecasts under general assumptions.
%
%   DM = dmtest(e1, e2, ...) calculates the D-M test statistic on the base 
%   of the loss differential defined by the function handle dfun, i.e., d =
%   dfun(e1,e2). For example, dfun = @(x,y) x.^2-y.^2 compares the squared
%   loss differential. 
%
%   In particular, with the DM statistic one can test the null hypothesis: 
%   H0: E(d) = 0. The Diebold-Mariano test assumes that the loss 
%   differential process 'd' is stationary and defines the statistic as:
%   DM = mean(d) / sqrt[ (1/T) * VAR(d) ]  ~ N(0,1),
%   where VAR(d) is an estimate of the unconditional variance of 'd'.
%
%   This function also corrects for the autocorrelation that multi-period 
%   forecast errors usually exhibit. Note that an efficient h-period 
%   forecast will have forecast errors following MA(h-1) processes. 
%   Diebold-Mariano use a Newey-West type estimator for sample variance of
%   the loss differential to account for this concern.
%
%   'e1' is a 'T-by-1' vector of the forecast errors from the first model
%   'e2' is a 'T-by-1' vector of the forecast errors from the second model
%   'dfun' is a function handle specifying the functional form of the loss 
%   differential. Default choice: dfun = @(x,y) x.^2-y.^2
%
%   DM = DMTEST(e1, e2, dfun, 'h') allows you to specify an additional parameter 
%   value 'h' to account for the autocorrelation in the loss differential 
%   for multi-period ahead forecasts.   
%       'h'         the forecast horizon, initially set equal to 1
%
%   DM = DMTEST(...) returns a constant:
%       'DM'      the Diebold-Mariano (1995) test statistic
%
%  Semin Ibisevic (2011) Yifan Li (2024)
%
% -------------------------------------------------------------------------
% References
% Diebold, F.X. and R.S. Mariano (1995), "Comparing predictive accuracy", 
% Journal of Business & Economic Statistics, 13, 253-263.
% -------------------------------------------------------------------------

if nargin <3
   dfun = @(x,y) x.^2-y.^2; 
end
if nargin < 4
   h = 1; 
end
if size(e1,1) ~= size(e2,1) || size(e1,2) ~= size(e2,2)
    error('dmtest:InvalidInput','Vectors should be of equal length');
end


% Initialization
T = size(e1,1);

% Define the loss differential
d = dfun(e1,e2);

% Ralculate the variance of the loss differential, taking into account
% autocorrelation.
dMean = mean(d);
gamma0 = var(d);
if h > 1
    gamma = zeros(h-1,1);
    for i = 1:h-1
        sampleCov = cov(d(1+i:T),d(1:T-i),1);
        gamma(i) = sampleCov(2);
    end
    varD = gamma0 + 2*sum(gamma);
else
    varD = gamma0;
end
%k is calculated to adjuste the statistic as per Harvey, Leybourne, and Newbold (1997) 
k = ((T+1-2*h+(((h)*(h-1))/T))/T)^(1/2);
% Retrieve the diebold mariano statistic DM ~N(0,1)
    DM = (dMean / sqrt ( (varD/T) ))*k;    
%P_VALUE is calculated    
    p_value = 2*tcdf(-abs(DM),T-1);
end
    


