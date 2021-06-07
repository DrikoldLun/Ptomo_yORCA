function [ out ] = sample_Jeffreys( Nsamp,k_max,k_min )
%[ out ] = sample_Jeffreys( Nsamp,k_max,k_min )
%   Take a random sample from a discrete random variable that follows a
%   Jeffrey's distribution between max/min (inclusive) values of k_max and
%   k_min. The Jeffrey's distribution weights the likelihood of each value
%   k as 1/k. Nsamp is the number of samples to draw.

if nargin < 3 || isempty(k_min)
    k_min = 1;
end

% possible values
k_vals = [k_min:k_max];
% non-normalised probability of each value
k_prob = 1./k_vals;
% normalised probability
k_probn = k_prob./sum(k_prob);
% cumulative sum, including zero
k_cumsum = [0,cumsum(k_probn)];

% random sample somewhere along the distribution
r = rand(Nsamp,1);
% output value is the largest integer in the cumulative probability
% distribution that r is greater than.
out = max(double(r > k_cumsum).*(ones(Nsamp,1)*[k_vals,nan]),[],2);


end

