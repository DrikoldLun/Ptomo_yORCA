function [maxA,x,y] = maxgrid(A)
% [maxA,x,y] = maxgrid(A)
%MAXGRID This script finds the (x,y) location and value of the maximum of a
%MxN matrix, e.g. an error grid where:
% x is the column number of the maximum - the x coordinate
% y is the row number of the maximum - the y coordinate
[~,x]=max(max(A));
[maxA,y]=(max(A(:,x)));
maxA;

end