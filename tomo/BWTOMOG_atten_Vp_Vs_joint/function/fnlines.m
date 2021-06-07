function [ nlines ] = fnlines( fid )
% FNLINES Get number of lines in file
%     NLINES = FNLINES(FID) returns the number of lines in the specified file.
%     
%     FID is an integer file identifier obtained from FOPEN.

nlines=0;
while fgets(fid)~=-1
    nlines=nlines+1;
end

% end

