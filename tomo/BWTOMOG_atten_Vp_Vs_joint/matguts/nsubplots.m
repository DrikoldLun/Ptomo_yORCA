function [nxsub,nysub] = nsubplots(nz)
% [nxsub,nysub] = nsubplots(nz)
    if nz <= 3 
        nxsub = nz;
        nysub = 1;
    elseif nz <= 6
        nxsub = ceil(nz/2);
        nysub = 2;
    elseif nz <= 12
        nxsub = ceil(nz/3);
        nysub = 3;        
    elseif nz <= 15
        nxsub = ceil(nz/3);
        nysub = 3;        
    elseif nz == 16
        nxsub = 4;
        nysub = 4;        
    elseif nz <= 20
        nxsub = ceil(nz/5);
        nysub = 4;        
    end
end 