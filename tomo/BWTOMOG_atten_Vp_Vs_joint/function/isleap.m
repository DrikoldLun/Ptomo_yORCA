function ifleap = isleap(year)
% ifleap = isleap(year)
% 
% Function to determine if the input year is a leap year according to the
% Gregorian calendar. Can accommodate vectors or matrices of "year". 

ifleap = ((rem(year,4)==0 & rem(year,100)~=0) | rem(year,400)==0);
