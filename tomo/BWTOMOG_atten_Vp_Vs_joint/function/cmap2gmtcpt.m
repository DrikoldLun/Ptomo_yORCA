function ofile = cmap2gmtcpt(vmax,vmin,colourmap,ofile)
% ofile = cmap2gmtcpt(vmax,vmin,colourmap,ofile)
% 
% function to take a matlab colour map, with some bounds, and turn it into
% a GMT-readabel .cpt file.

if nargin<3
    colourmap = parula;
end

if nargin<4
    ofile = 'matlab.cpt';
end

vals = vmin:vmax;
cols = colour_get(vals,vmax,vmin,colourmap);

fid = fopen(ofile,'w');
for ival = 1:length(vals)-1
    fprintf(fid,'%f %.0f/%.0f/%.0f %f %.0f/%.0f/%.0f\n',...
        vals(ival),cols(ival,:)*255,vals(ival+1),cols(ival+1,:)*255);
end

end