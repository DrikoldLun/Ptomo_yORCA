function [ cc ] = cmap_makecustom( c1, c2, wtfrac )
%[ cc ] = cmap_makecustom( c1, c2, wtfrac )
%   Makes custom colourmap scaling linearly between the two endmember
%   colours, c1 and c2, through white, with an optional fraction of
%   whitespace in between
% 
% INPUT:
%   c1      Upper endmember colour as 3-component vector [r g b]
%   c2      Lower endmember colour as 3-component vector [r g b]
%   wtfrac  fraction of white space of total [default is 0]
%             e.g. for -5<x<5, wtfrac=0.2 means region -1<x<1 will be white
% 
% OUTPUT:
%   cc      colourmap: a (200 x 3) matrix, can be called with colormap(cc)
% 
% Written by Zach Eilon, 2012


if nargin < 3
    wtfrac=0;
end

N = 200;
W = N*wtfrac/2;
M = N/2 - W;

rampu = linspace(0,1,M+1)'; rampu = rampu(1:end-1);
rampd = flipud(rampu);
flatm = ones(M,1);
whitw = ones(2*W,1);

d1 = 1-c1;
d2 = 1-c2;

rr = [c1(1)*flatm + d1(1)*rampu; whitw; c2(1)*flatm + d2(1)*rampd];
gg = [c1(2)*flatm + d1(2)*rampu; whitw; c2(2)*flatm + d2(2)*rampd];
bb = [c1(3)*flatm + d1(3)*rampu; whitw; c2(3)*flatm + d2(3)*rampd];

cc = [rr gg bb];
end

