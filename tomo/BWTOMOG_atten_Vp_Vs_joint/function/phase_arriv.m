function [ arrtime ] = phase_arriv( phase,gcarc,edep,arrn )
%function [ arrtime ] = phase_arriv( phase,gcarc,edep,arrn )
%   
% Function to quickly output the predicted arrival time of a seismic phase
% at a station, given the following:
% INPUTS:
% phase = name of phase (string)
% gcarc = great circle distance between earthquake and station
% edep = depth of event
% arrn (optional) = if multiple arrivals with name "phase" then which one

if nargin<4
    arrn=1;
end

tt=taupTime([],edep,phase,'deg',gcarc);
arrtime=tt(arrn).time;

end

