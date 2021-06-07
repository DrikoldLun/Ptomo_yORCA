function a = minmax(t)
%  a = minmax(t)
% 
%  Function to spit out 2-element vector a containing [min;max] of a time
%  series, tt. 
% 
% If tt is a matrix, a is a 2xNcols matrix with column-wise minima of in
% first row, and maxima in second row.

    a = [min(t);max(t)];

end