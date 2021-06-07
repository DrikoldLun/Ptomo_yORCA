function [ synth_model ] = make_synth_model( par,data )

fprintf('>  Creating synthetic model\n');

if par.PS == 3
    moddim = 2;
else
    moddim = 1;
end

if strcmp(par.sym.opt,'custom')

    mval = zeros(par.nmodel,1);

    na =  length(par.sym.adval);
    for ia = 1:na
        ind1 = par.mx <= par.sym.acx(ia) + 0.5*par.sym.awx(ia);
        ind2 = par.mx >= par.sym.acx(ia) - 0.5*par.sym.awx(ia);
        ind3 = par.my <= par.sym.acy(ia) + 0.5*par.sym.awy(ia);
        ind4 = par.my >= par.sym.acy(ia) - 0.5*par.sym.awy(ia);
        ind5 = par.mz <= par.sym.acz(ia) + 0.5*par.sym.awz(ia);
        ind6 = par.mz >= par.sym.acz(ia) - 0.5*par.sym.awz(ia);
        ind = ind1 & ind2 & ind3 & ind4 & ind5 & ind6;

        mval(ind) = mval(ind) + 0.01*par.sym.adval(ia);
    end
    
elseif strcmp(par.sym.opt,'checker')

    checker_model = make_checker_model(par,par.sym.dval);
    mval = checker_model.mval;

elseif strcmp(par.sym.opt,'halves')

    mval = zeros(par.nmodel,1);
    % left side
    mval(par.mln <= -5) = mval(par.mln <= -5) -0.10;
    % right side 
    mval(par.mln >= 5) = mval(par.mln >= 5) + 0.10;

end

%% account for P/S
if par.PS == 3 && par.Rdvpdvs == 0 % P and S and they're not just scaled 
    % assume that model defined above is the S model
    mvalS = mval;
    % scale to P
    Rdvpdvs = 0.55; % this is the value of dlnVp/dlnVs
    mvalP = (Rdvpdvs*mvalS)./(1 + (1-Rdvpdvs)*mvalS); %precise scaling for slownesses (since denominator close to 1, will be ignored in inversion scaling) 
    % concatenate, P model first
    mval = [mvalP;mvalS];
end
    
%add a little random noise
mval = mval + random('norm',0,0.00001,size(mval));
    
%% insert into model
synth_model.mval = mval;

% random static terms
synth_model.estatic = random('norm',0,par.sym.estatic_sd,data.evt.nevts,1);
synth_model.sstatic = random('norm',0,par.sym.sstatic_sd,data.stn.nstas*moddim,1);


end

