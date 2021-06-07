function plot_hitq( par,opt )
% plot_hitq( par, [opt=1])
%   plot the hit quality (opt=1) or count (opt=2)


if nargin < 2 || isempty(opt)
    opt=1;
end
if par.PS == 1 
    Nph = 1; compstr = {'P'};
elseif par.PS == 2
    Nph = 1; compstr = {'S'};
elseif par.PS == 3
    Nph = 2; compstr= {'P','S'};   
end

for ip = 1:Nph

    shape = [par.ny,par.nx,par.nz];

    lats= reshape(par.mlt,shape);
    lons= reshape(par.mln,shape);
    hitq = reshape(par.hitq(:,ip),shape);
    hitc = reshape(sum(par.hitc(:,:,ip),2),shape);

    if opt == 1
        hit = 100*hitq;
        str = 'Hit Quality ';
    elseif opt == 2
        hit = hitc;
        str = 'Hit Count ';
    end
    datstr= {'$\delta T$  ','$\Delta t^*$ '};    


    figure(35+ip), clf
        set(gcf,'Position',[0,0,800,800])
    [nxsub,nysub] = nsubplots(par.nz);
    for iz = 2:par.nz-1
        subplot(nxsub,nysub,iz-1)
        hold on
        contourf(lons(:,:,iz),lats(:,:,iz),hit(:,:,iz),[0:10:100],'linestyle','none')
    %     xlabel('Longitude'); 
    %     ylabel('Latitude');
        title(sprintf('Depth %.0f km',par.zz(iz)),'FontSize',14)
        xlim(par.plot_lonlims); ylim(par.plot_latlims)
        shading flat 
        caxis([0 100])


        % label the type of data
        if iz==par.nz-1
            xlabel([str,datstr{par.t_ts},compstr{ip}],...
                'interpreter','latex','FontSize',20)
        end


    end % loop on depths
    
    subplot(nxsub,nysub,nxsub*nysub)
    colorbar

end % loop on phases

end

