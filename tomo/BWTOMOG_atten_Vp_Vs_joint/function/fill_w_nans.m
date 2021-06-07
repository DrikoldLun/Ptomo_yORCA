function fill_w_nans(x,y,varargin)
% fill_w_nans(x,y,c,varargin)
%
% quick n dirty wrapper function to take vectors x,y, that represent broken
% up areas, separated by nans, and plot them with fill
indx = find(isnan(x));

ih = ishold;
if ~ih, hold on, end

for ii = 1:length(indx)-1
    ind1 = indx(ii)+1;
    ind2 = indx(ii+1)-1;
    fill(x(ind1:ind2),y(ind1:ind2),varargin{:})
end

if ~ih, hold off, end
