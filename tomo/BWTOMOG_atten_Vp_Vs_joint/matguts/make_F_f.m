function [ F,f,data ] = make_F_f( G,d,data,par )
% [ F,f,data ] = make_F_f( G,d,data,par )
% 
% create regularised, weighted matrices to invert, using the principle
% 
% [Wd * G]     [Wd*d]       Wd is the data uncertainty
% [      ] m = [    ]
% [Wm * H]     [Wm*h]       Wm is the weighting for the regularisation
% 
% or,    F m = f


if par.PS == 3, sdim = 2; else, sdim = 1; end % account for extra station terms

nstas = data.stn.nstas;   
nrays = data.ray.nrays;
nevts = data.evt.nevts;


%% H_damp
[H_damp,h_damp]  = make_dampmat(data,par);

%% H_smooth
global dbname
if par.t_ts == 1
    Hsmthfile = [dbname,'/H_smth_',par.phasecomp,'_v']; 
elseif par.t_ts == 2
    Hsmthfile = [dbname,'/H_smth_',par.phasecomp,'_q'];
end

if par.build_smooth == 0
    load(Hsmthfile);
elseif par.build_smooth == 1
   [ H_smth ] = make_smoothmat( par ); 
   save (Hsmthfile,'H_smth','-v7.3');
end  
% Account for P and S if needed
if par.PS == 3 && par.Rdvpdvs==0
    % smooth both P and S
    Z_smth = sparse(size(H_smth,1),size(H_smth,2));
    wt = par.VpVs_smooth_same;
    H_smth = [H_smth  Z_smth   % P smoothing
              Z_smth  H_smth   % S smoothing
          wt*[H_smth -H_smth]]; % P and S second deriv the same
end
H_smth = [H_smth sparse(size(H_smth,1),nevts + nstas*sdim)]; % add station+evt columns



%% Data kernel
WG = spdiags(data.ray.wt,0,nrays,nrays)*G;
Wd  = data.ray.wt.*d;

% scale regularisation
if par.scalereg
    LWdGdm = svds(WG,1);
    LH_damp = svds(H_damp,1);
    LH_smth  = normest(H_smth,1e-3);
% LH_smth = svds(H_smth,1);
    
    A = LWdGdm/LH_damp;
    B = LWdGdm/LH_smth;
else
    A = 1;
	B = 1;
end

% fprintf('norm dGdm = %.3f\n',LWdGdm);
% fprintf('norm LH_damp = %.3f\n',LH_damp);
% fprintf('norm LH_smth = %.3f\n',LH_smth);
% fprintf('scale_damp = %.3f\n',par.damp*A);
% fprintf('scale_smth = %.3f\n',par.smooth*B);
% fprintf('norm damp = %.3f\n',par.damp*A*LH_damp);
% fprintf('norm smth = %.3f\n',par.smooth*B*LH_smth);



%% Make F, f
F = [WG; A*H_damp; B*par.smooth*H_smth];
f = [Wd; A*h_damp; zeros(size(H_smth,1),1)];

% pt_k = [45,-125,80]
% [ r_k ] = Rmatrix( F, par, pt_k )
% pause
end

