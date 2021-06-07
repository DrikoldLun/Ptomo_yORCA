function [MM,DD] = calday(YYYY,JJJ,mode)
%[MM,DD] = calday(YYYY,JJJ,mode)
% converts day of year into calendar date (month & day)
% if mode is 'm', will output only month (as a double)
% if mode is 'd', will output only day (as a double)
% if mode is 's', will output two 2-character strings
% if mode is 'i', will output two numbers
%
% Written by Zach Eilon, 2011
% edited to include only month/only day output option in 2018


if nargin==1
    JJJ=YYYY;
    YYYY=2010;
end

% vectorise if only one year or only one day
if length(YYYY)==1 && length(JJJ)>1
    YYYY = YYYY*ones(size(JJJ));
elseif length(YYYY)>1 && length(JJJ)==1
    JJJ = JJJ*ones(size(YYYY));
end

% normal year matrix
A=ones(31,12);
dpm = [31 28 31 30 31 30 31 31 30 31 30 31]';
for k=1:length(dpm)
    A(dpm(k)+1:end,k)=0;
end
A(A>0)=[1:sum(dpm)]';
Anorm = A;

% leap year matrix
A=ones(31,12);
dpm = [31 29 31 30 31 30 31 31 30 31 30 31]';
for k=1:length(dpm)
    A(dpm(k)+1:end,k)=0;
end
A(A>0)=[1:sum(dpm)]';
Aleap = A;

DD = nan(size(JJJ));
MM = nan(size(JJJ));

for ii = 1:numel(JJJ);
    if isleap(YYYY(ii))
        [DD(ii),MM(ii)]=find(Aleap==JJJ(ii));
    else
        [DD(ii),MM(ii)]=find(Anorm==JJJ(ii));
    end
end

% % leap?
% ifleap = isleap(YYYY);
% 
% [~,ian] = intersect(Anorm,JJJ(~ifleap),'stable');
% [~,ial] = intersect(Aleap,JJJ( ifleap),'stable');
% 
% 
% DD(ifleap) = rem(ial,31);
% MM(ifleap) = ceil(ial/31);
% 
% DD(~ifleap) = rem(ian,31);
% MM(~ifleap) = ceil(ian/31);

% [DD,MM]=find(A==JJJ);

if nargin==3
    if strcmp(mode,'s')
        if MM >= 10
            MM=sprintf('%s',num2str(MM));
        elseif MM < 10 
            MM=sprintf('0%s',num2str(MM));
        end
        if DD >= 10
            DD=sprintf('%s',num2str(DD));
        elseif DD < 10 
            DD=sprintf('0%s',num2str(DD));
        end
    elseif strcmp(mode,'m')
        DD = [];
    elseif strcmp(mode,'d')
        MM = DD;
        DD = [];
    else
        error('\tMode "%s" not recognised',mode);
        return
    end
end


if nargout<2
    MM = [MM,DD];
end

end