function mu_w = mean_wtd(X,wts)
% mu_w = mean_wtd(X,wts)
% 
% Function to calculate a weighted mean, whereby the mean of the data
% vector, X, is calculated by weighting each datum by the value in the
% corresponding weights vector, wts.

if nargin<2 || isempty(wts) || ~any(wts~=1)
    mu_w = mean(X);
else
    totXw = sum(X.*wts);
    totw = sum(wts);
    mu_w = totXw./totw;
end



end