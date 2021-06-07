function [X,Y,H,R] = rosediag(thetas,binw,scale)
% USEAGE:   [ X , Y [,H] [,R] ] = rosediag( thetas [,binwidth] )
% 
% INPUTS
%   thetas - 1xN cell array containing Mx1 vectors of angles, in degrees. 
%          - Each column will be a different colour/histogram, all the
%           columns will be layered in triangles on an area-proportional
%           rose histogram (i.e. simplest if theta is just column vector)
%   binw   - width of histogram bins, in degrees. Must divide 360 evenly
%           into L bins (if not, will be rounded to nearest integer L)
%          - optional argument - default is 30deg
%   scale  - histogram scale length - i.e. bin heights (really bin radii).
% 
% OUTPUTS
%   X      - 1xN cell array, each cell contains Lx5 matrix of x coordinates
%           for rose diagram patches for that one of the N data series
%   X      - 1xN cell array, each cell contains Lx5 matrix of y coordinates
%           for rose diagram patches for that one of the N data series
%   H      - LxN matrix of number of entries in each of L bins for each of
%           N data series 
%   R      - LxN matrix of rose diag radii in each of L bins for each of N
%           data series 
% 
%  Written by Z. Eilon    14 July 2015

if nargin < 2
    binw = 30;
end
if nargin < 3;
    scale = 1;
end
if ~iscell(thetas)
    junk{1} = thetas;
    thetas=junk;
end

%% some prelims
N = size(thetas,2);
L = round(360/binw);
hw = 180/L;
bincs = linspace(hw,360-hw,L)';
schw = sind(hw)*cosd(hw);
phi1 = bincs - hw;   phi2 = bincs + hw;
cp1 = cosd(phi1); sp1 = sind(phi1);    
cp2 = cosd(phi2); sp2 = sind(phi2);

%% results structures
H = zeros(L,N);
R = zeros(L,N);


%% histogram
for kk = 1:N
H(:,kk) = hist(thetas{kk},bincs);
end

%% rose diag Radii for each bin
R(:,1) = sqrt(scale*H(:,1)/schw );
for kk = 2:N
    R(:,kk) = quadratic_solve(1, 2.*sum(R(:,1:kk-1),2), -scale.*H(:,kk)./schw);
end
sumR = [zeros(L,1), cumsum(R,2)];


for kk = 1:N
    Xtemp = [sumR(:,kk).*sp1 sumR(:,kk+1).*sp1 sumR(:,kk+1).*sp2 sumR(:,kk).*sp2 sumR(:,kk).*sp1];
    Ytemp = [sumR(:,kk).*cp1 sumR(:,kk+1).*cp1 sumR(:,kk+1).*cp2 sumR(:,kk).*cp2 sumR(:,kk).*cp1];
    
    X{kk} = Xtemp;
    Y{kk} = Ytemp;
end

%% to plot, uncomment this:
% figure(1),clf, hold on, axis equal
% for kk = 1:N
% for ii = 1:L
% fill(X{kk}(ii,:),Y{kk}(ii,:),colour_get(kk,1,N)')
% end
% end

if N==1;
    X = X{1};
    Y = Y{1};
end




% end