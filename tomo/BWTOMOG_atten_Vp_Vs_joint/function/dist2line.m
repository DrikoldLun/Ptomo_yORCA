function [ distance ] = dist2line( Q1,Q2,P )
%[ distance ] = dist2line( Q1,Q2,P )
%   calculates the distance from a line between points Q1 and Q2 of a
%   series of points, P in a NxD vector (where D is 2 for a plane, and 3
%   for 3-space
% Q1 and Q2 are 1xD vectors defining the line ends
% 
% Written by Zach Eilon, 2013


N = size(P,1);

if iscolumn(Q1), Q1=Q1'; end
if iscolumn(Q2), Q2=Q2'; end

if numel(Q1)==2; Q1= [Q1,0]; end
if numel(Q2)==2; Q2= [Q2,0]; end
if size(P,2)<3; P = [P,zeros(N,1)]; end

A = ones(N,1)*(Q2-Q1);
B = P - ones(N,1)*Q1;
C = cross(A,B);

distance = -C(:,3)/norm(Q2-Q1);

end

