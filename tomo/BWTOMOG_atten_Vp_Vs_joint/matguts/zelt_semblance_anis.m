function [semblance,npts] = zelt_semblance_anis(par)
radius = 100;
semblance = zeros(par.nmodel,1);
npts = zeros(par.nmodel,1);
for i = 1:par.nmodel
    dist = sqrt((par.mx-par.mx(i)).^2+(par.my-par.my(i)).^2+(par.mz-par.mz(i)).^2);
    ind = find(dist<=radius);
    npts(i) = length(ind);
    upper = sum((par.si_model.mval(ind)+par.so_model.mval(ind)).^2); 
    lower = sum(par.si_model.mval(ind).^2+par.so_model.mval(ind).^2);
    semblance(i) = 0.5*upper/lower;
end