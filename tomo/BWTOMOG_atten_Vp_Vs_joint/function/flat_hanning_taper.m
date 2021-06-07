function win = flat_hanning_taper(vec_t,tapertime)
% win = flat_hanning_taper(vec_t,tapertime)
% function to create a box-car window with a hanning taper end.
% INPUTS:
%     vec_t     = Nx1 vector of input times, equally spaced by dt
%     tapertime = amount of time ON EACH END to taper over
% OUTPUTS:
%     win       = window ramping up from 0 to 1 over tapertime, flat
%                 in the middle, ramping down over tapertime at the end

win=ones(size(vec_t));

%ndata=length(vec_t);
dt=vec_t(2)-vec_t(1);


taperlength = floor(tapertime/dt);

taper = hanning(taperlength*2);
win(1:taperlength) = taper(1:taperlength);
win(end-taperlength+1:end) = taper(taperlength+1:end);

end
