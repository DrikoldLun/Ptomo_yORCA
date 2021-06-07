function [ZZ,ZX,ZY] = gauss2D(X,sX,Y,sY,xx,yy)
% function [ZZ,ZX,ZY] = gauss2D(X,sX,Y,sY,xx,yy)
% 
% function to create a surface with a peak at values X and Y with x-width
% described by sX and y-width described by sY and resolved onto a grid with
% dimensions xx and yy.

[ZX,ZY] = meshgrid(xx,yy);

ZZ = exp(-( ((ZX - X).^2)/(2*sX^2) + ((ZY - Y).^2)/(2*sY^2) ) );