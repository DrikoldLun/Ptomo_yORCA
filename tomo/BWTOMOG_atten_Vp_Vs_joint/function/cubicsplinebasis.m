function [f,f_norm] = cubicsplinebasis(x,x0,xw)
% [f,f_norm] = cubicsplinebasis(x,x0,xw)
% 
%  function to compute cubic splines that comprise basis set spanning the
%  space of input vector x. Each of the splines is a two-sided function of
%  the form 3x^2 - 2x^3. The edges have one-sided functions. The sum of the
%  splines is one at every point. The basis is non-orthogonal. 
% 
% INPUTS:
%  x   - (npts,1) vector of independent variable space;
%  xo  - (nspl,1) vector of spline centres (strongly recommend this
%        includes min(x) and max(x) and has spacing wx/2
%  wx  - scalar or (nspl,1) vector of full spline widths.
% 
%       recommended x etc:
%        x = [minx:dx:maxx]'; % dx is the independent variable spacing
%        xw = 2*(max(x)-min(x))/(nspl-1)
%        x0 = [min(x):xw/2:max(x)];
% 
% OUTPUTS:
%  f      - (npts,nspl) matrix of spline values, zero except within each
%           column's bin
%  f_norm - (npts,nspl) matrix of normalised spline values, zero except
%           within each column's bin, and summing to one along each column\
% 
%       Resolving a function to this basis: if you have a function, y, the
%       same size as x, you can integrate to get the basis coefficients:
%           coeffs = f_norm'*y;
%       and then simply multiply these by the basis functions and sum:
%           yint = sum((ones(npts,1)*coeffs').*f,2); 
%       where the outer product is the handy MATLAB trick to matricise two
%       vectors
% 
% 
% Z. Eilon 08/2016




if length(x0)>1;
    [x0,xx] = meshgrid(x0,x);
    if length(xw)>1;
        xw = xw(:)';
        [xw,xx] = meshgrid(xw,x);
    end      
x = xx;
end

f = 1 - 3*(abs(x-x0)./xw/0.5).^2 + 2*(abs(x-x0)./xw/0.5).^3;
f(x > x0+xw/2) = 0;
f(x < x0-xw/2) = 0;  

f_norm = f./(ones(length(x),1)*sum(f)); % normalise by integral of each basis fn


end

