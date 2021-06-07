function [ XX ] = voronoi( X,val,xrng )
% [ XX ] = voronoi( X,val,xrng )
% 
%  Function to make a velocity profile from a set of voronoi elements. 
%  Only in 1-D for now.

if nargin<3
    xrng = [0,max(X)];
end
xrng = sort(xrng);
val = val(:);
X = X(:);
% 
% X = sort(50*rand(5,1));
% val = [3.2, 3.9, 4.1, 4.5, 4.6]';
% xrng = [0,100];

N = length(val);
Xtop = [xrng(1);midpts(X)];
Xbot = [midpts(X);xrng(2)];
XX = struct('top',Xtop,...
            'bot',Xbot,...
            'val',val,...
            'Xplot',reshape([Xtop';Xbot'],[2*N,1]),...
            'Vplot',reshape([val';val'],[2*N,1]));



end

