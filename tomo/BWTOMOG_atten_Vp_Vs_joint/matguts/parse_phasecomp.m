function [datfiles,stafiles,phases] = parse_phasecomp(tomodatdir,phasecomp,dattype)
%[datfile,stafile,phases] = parse_phasecomp(tomodatdir,phasecomp,dattype)
%
% function to parse phases and collect all files to be read

if ~strcmp(tomodatdir(end),'/'), tomodatdir = [tomodatdir,'/']; end

a = textscan(phasecomp,'%s','delimiter','_'); a = a{:};

if length(a) == 1 %single component
    datfiles{1} = sprintf('%sdata_%s_%s.dat',tomodatdir,dattype,phasecomp);
    stafiles{1}  = sprintf('%sstations_%s.dat',tomodatdir,phasecomp);  
    phases{1} = a{1}(1);
else
    for ic = 1:length(a)-1
        datfiles{ic,1} = sprintf('%sdata_%s_%s.dat',tomodatdir,dattype,a{ic+1});
        stafiles{ic,1} = sprintf('%sstations_%s.dat',tomodatdir,a{ic+1});
        phases{ic,1} = a{ic+1}(1);
    end
end


end

