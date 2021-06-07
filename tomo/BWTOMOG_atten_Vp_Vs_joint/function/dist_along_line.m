function [ distance ] = dist_along_line( Q1,Q2,P )
%[ distance ] = dist_along_line( Q1,Q2,P )
%   projects a series of points P onto a line between enpoints Q1 and Q2
%   and calculates the distance each of the projected points from Q1
%   P in a NxD vector (where D is 2 for a plane, and 3 for 3-space.
%   Q1 and Q2 are 1xD vectors defining the line ends
% 
% Written by Zach Eilon, 2015


N = size(P,1);

Q1=Q1(:)';
Q2=Q2(:)';

if numel(Q1)==2; Q1= [Q1,0]; end
if numel(Q2)==2; Q2= [Q2,0]; end
if size(P,2)<3; P = [P,zeros(N,1)]; end




distance = dot(ones(N,1)*(Q2-Q1),P-(ones(N,1)*Q1),2)/norm(Q2-Q1);


end

