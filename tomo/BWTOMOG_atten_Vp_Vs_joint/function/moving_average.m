function [ B ] = moving_average( A,n,dim )
%[ B ] = moving_average( A,n[=5],dim[=1] )
% tapered, n-point moving average in direction dim
% n must be integer, and will be rounded up if even
% 
%   Z. Eilon,   March 2017, edited Sept 2019

if nargin < 2 || isempty(n)
    n = 5; 
end
if nargin < 3 || isempty(dim)
    dim = 1; 
end


if ~isodd(n)
    n = round_level(n,2)+1;
end

if isrow(A), A = A(:); end

dy = 0.5*(n-1);
L = size(A,dim);

% mid full moving average section
C = toeplitz([1,zeros(1,L-n)],[ones(1,n),zeros(1,L-n)]);

% taper section % wish I could think of a way to avoid the loop!
D = zeros(dy,L);
f = 2*[1:dy]'-1;
for ii = 1:length(f)
    D(ii,1:f(ii)) = 1;
end

%complete moving average filter (not normalised)
G = [D;C;rot90(D,2)];

% do in right direction!
if dim == 1
    B = G*A./sum(G,2); % multiply and normalise
elseif dim == 2
    B = [G*A'./sum(G,2)]'; % multiply and normalise
end


end

