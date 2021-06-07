function [minA,x,y] = mingridc(A, confint)
% [minA,x,y] = mingrid(A)
%MINGRID This script finds the (x,y) location and value of the minimum of a
%MxN matrix, e.g. an error grid where:
% x is the column number of the minimum - the x cooridinate
% y is the row number of the minimum - the y coordinate
% in addition, this script returns the bounds of the confidence interval
% about the minimum defined by being less than some multiple of the minimum
% e.g. if the minimum value were 10 and the confint were 1.1, the function
% would also output the indices for the points along horizontal and
% vertical lines that were less than or equal to 11.
% the outputs are then x=[lbx minx ubx] and y=[lby miny uby] (upper and lower
% bounds)
    [~,minx]=min(min(A));
    [minA,miny]=(min(A(:,minx)));
if nargin==1
    x=minx;
    y=miny;
elseif nargin==2
horzrange=find(A(miny,:) <= minA*confint);
vertrange=find(A(:,minx) <= minA*confint);
if isempty(horzrange)==1; horzrange=minx; end
if isempty(vertrange)==1; vertrange=miny; end
x=[min(horzrange) minx max(horzrange)];
y=[min(vertrange) miny max(vertrange)];
end

