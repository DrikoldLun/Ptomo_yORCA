function [ h ] = title_custom( title_string,title_y,title_x, varargin )
%[ ax ] = title_custom( title_text,title_y,title_x,varargin{pairs} )
%   function to plot a title in a custom position within the plot.
%  usage:
%  title_custom('titstring',Name,Value...)
%  title_custom('titstring',title_y,Name,Value...)
%  title_custom('titstring',title_y,title_x,Name,Value...)

% defaults
ypos = 0.95;
xpos = 0.5;
fontsize = 20;
interpreter = 'latex';
fontweight = 'bold';
horal = 'center';
veral = 'middle';


if nargin<2, 
    title_x = xpos;
    title_y = ypos;
end
if nargin<3
    title_x = xpos;
end
if nargin>=3
    if ischar(title_y)
        varargin = {title_y,title_x,varargin{:}};
        title_y = ypos;
        title_x = xpos;
    elseif ischar(title_x)
        varargin = {title_x,varargin{:}};
        title_x = xpos;
    end
    if rem(length(varargin),2)~=0, error('need even # of varargin'); end
    for iv = 1:length(varargin)/2;
        switch varargin{2*iv-1}
            case 'fontsize'
                fontsize = varargin{2*iv};
            case 'fontweight'
                fontweight = varargin{2*iv};
            case 'interpreter'
                interpreter = varargin{2*iv};
            case 'horizontalalignment'
                horal = varargin{2*iv};
            case 'varticalalignment'
                veral = varargin{2*iv};
        end
    end
        
    
end
ypos = title_y;
xpos = title_x;

if strcmp(interpreter,'latex'), title_string = regexprep(title_string,'_','\\_'); end

ax99 = axes('pos',[0 0 1 1]);
set(ax99,'visible','off')
if strcmp(interpreter,'latex') && strcmp(fontweight,'bold')
    h = text(ax99,xpos,ypos,['\textbf{',title_string,'}'],...
        'fontsize',fontsize,'interpreter',interpreter,...
        'horizontalalignment',horal,'verticalalignment',veral);
else
   h = text(ax99,xpos,ypos,title_string,...
    'fontsize',fontsize,'interpreter',interpreter,'fontweight',fontweight,...
    'horizontalalignment',horal,'verticalalignment',veral); 
end




end

