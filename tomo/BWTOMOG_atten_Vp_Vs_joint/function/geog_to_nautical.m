function [ nautical_lat_or_lon ] = geog_to_nautical( geo_lat_or_lon )
%[ nautical_lat_or_lon ] = geog_to_nautical( geo_lat_or_lon )
%   Quick function to take a coordinate in decimal degrees and convert it
%   to a string of degrees, minutes and decimal minutes

deg = fix(geo_lat_or_lon);
decdeg = geo_lat_or_lon - deg;
mindec = 60*abs(decdeg);

nautical_lat_or_lon = [num2str(deg),' ',num2str(mindec)];


end

