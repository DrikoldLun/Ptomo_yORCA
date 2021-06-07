function h = add_state_boundaries(ax,latlims,lonlims,col)
% h = add_state_boundaries(ax,latlims,lonlims,colour)
%   function to plot light grey state boundaries on maps

if nargin ==0 || isempty(ax)
    ax = gca;
end

if nargin < 2 || isempty(latlims)
    axl = axis(ax);
    latlims = axl(3:4);
end

if nargin < 3 || isempty(lonlims)
    axl = axis(ax);
    lonlims = axl(1:2);
end
if nargin<4 || isempty(col)
    col = [0.6 0.6 0.6];
end

    
latlims = latlims(:); % make column
lonlims = lonlims(:); % make column

states = shaperead('usastatehi','UseGeoCoords', true, 'BoundingBox', [lonlims,latlims]);

hrat = get(gca,'DataAspectRatio');
h = geoshow(ax, states, 'Edgecolor', col,'linestyle','-','facecolor','none','linewidth',2);
set(gca,'DataAspectRatio',hrat);



end

