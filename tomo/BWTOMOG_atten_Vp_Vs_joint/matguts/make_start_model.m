function [ model, par ] = make_start_model( par,data)
%[ model, par ] = make_start_model( par,nevts,nstas  )
%  make the starting model for the inversion


nevts = data.evt.nevts;
nstas = data.stn.nstas;

% default dimensions
moddim = 1; % multiple of nodes to model parms
sdim = 1;   % multiple of stas to station terms

if par.PS == 3
    if par.Rdvpdvs==0
        moddim = 2;
    end
    sdim = 2;
end

model = struct('mval',zeros(par.nmodel*moddim,1));
           
% positive e_stat or s_stat will decrease predicted time.
% therefore these are additional delays - positive means slow structure
% beneath just that station.
           
% static terms
model.estatic = zeros(nevts,1); 
model.sstatic = zeros(nstas*sdim,1); 

par.start_model = model;
end

