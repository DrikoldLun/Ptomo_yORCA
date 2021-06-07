function [filepath] = save2eps(fignumber,filename,odir)
% function [filepath] = save2eps(fignumber,filename,odir)
% 
% function to save a figure to a good-quality jpg
%
% INPUTS:
% fignumber - (required) number of figure to save
% filename - (optional, default is figN) name of saved file
% odir - (optional, default is current dir) dir to save file in

if nargin < 3
    odir = pwd;
end
if nargin < 2
    filename = sprintf('fig%u',fignumber);
end

% thisdir = pwd;
% odir = cd(odir);
% cd(thisdir)

if strcmp(odir(end),'/')~=1
    odir = strcat(odir,'/');
end

h = figure(fignumber);
print(h,'-painters','-depsc2','-r200',strcat(odir,filename))

filepath = strcat(odir,filename,'.eps');
end