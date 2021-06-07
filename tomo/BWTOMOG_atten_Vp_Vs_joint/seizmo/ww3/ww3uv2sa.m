function [s]=ww3uv2sa(s,varargin)
%WW3UV2SA    Converts wind in u/v to speed & azimuth
%
%    Usage:    s=ww3uv2sa(s)
%
%    Description:
%     S=WW3UV2SA(S) converts the wind data in WW3 struct S from meters per
%     second in the U/V directions to speed in meters per second and
%     azimuth in degrees from 0-360.
%
%    Notes:
%
%    Examples:
%     % Read, Convert, Map:
%     ww3map(ww3uv2sa(ww3struct));
%
%    See also: WW3STRUCT, WW3REC, WW3CAT, WW3MAP, WW3MAPMOV, WW3MOV,
%              PLOTWW3, PLOTWW3TS, WW3BAZ2AZ

%     Version History:
%        Feb.  5, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  5, 2014 at 00:40 GMT

% todo:

% check ww3 input
if(nargin==0) % gui selection of grib file
    s=ww3struct();
elseif(isstruct(s))
    valid={'path' 'name' 'description' 'units' 'data' ...
        'lat' 'lon' 'time' 'latstep' 'lonstep' 'timestep'};
    if(any(~ismember(valid,fieldnames(s))))
        error('seizmo:ww3uv2sa:badWW3',...
            'S must be a struct as generated by WW3STRUCT!');
    end
elseif(ischar(s)) % filename given
    s=ww3struct(s,varargin{:});
else
    error('seizmo:ww3uv2sa:badWW3',...
        'FILE must be a string!');
end

% loop over struct elements and convert wind data when found
for i=1:numel(s)
    % find wind data from description field values
    if(isequal(s(i).description,{'u wind'  'v wind'}) ...
            || isequal(s(i).description,...
            {'U-component of wind'  'V-component of wind'}))
        s(i).description={'wind speed'  'wind azimuth'};
        s(i).units={'m/s' 'degrees'};
        s(i).data={sqrt(s(i).data{1}.^2+s(i).data{2}.^2) ...
            180/pi*angle(s(i).data{2}+1i.*s(i).data{1})};
        s(i).data{2}(s(i).data{2}<0)=s(i).data{2}(s(i).data{2}<0)+360;
    % components switched?
    elseif(isequal(s(i).description,{'v wind'  'u wind'}) ...
            || isequal(s(i).description,...
            {'V-component of wind'  'U-component of wind'}))
        s(i).description={'wind speed'  'wind azimuth'};
        s(i).units={'m/s' 'degrees'};
        s(i).data={sqrt(s(i).data{1}.^2+s(i).data{2}.^2) ...
            180/pi*angle(s(i).data{1}+1i.*s(i).data{2})};
        s(i).data{2}(s(i).data{2}<0)=s(i).data{2}(s(i).data{2}<0)+360;
    end
end

end
