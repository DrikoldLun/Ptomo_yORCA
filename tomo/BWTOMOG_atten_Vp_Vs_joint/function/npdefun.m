function [ npde ] = npdefun(Kest,theta,hscale,thrange)
% function [ npde ] = npdefun(Kest,theta,thrange)
% This function makes a non-parametric pdf and cdf using the algorithm in
% Fisher's Statistical analysis of circular data, section 2.2
% Inputs:
% Kest  ? the estimated concentration parameter
% theta ? the set of angles being paramaterised
% hscale ? a scaling factor for h, in the range 0.25-1.5 (smooths the dibn)
% thrange ? (1)==> -? to +? (2)==> 0 to +2?
% 
% the output is as follows:
% C1 are the angles, in radians
% C2 are the probability density values for radians: MUST CONVERT using d2r if
% you convert angles from r2d for the plot - to preserve area under curve
% C3 are densities plotted around a circle of unit radius (r-value)
% C4 are the cumulative density values

n=length(theta);
zeta = 1/sqrt(Kest);
h0 = sqrt(7)*zeta/(n^0.2);
h = hscale*h0;
npde=zeros(201,4); % c1 is 0<th<2pi, c2 is f(th), c3 is f*(th) - f scaled to circle of radius r
for j=0:200
    th=j*2*pi/200;
    i=0;
    SUM=0;
    while i < n
        i=i+1;
        d=abs(th-theta(i));
        e=min(d,2*pi - d);
        if e < h
        SUM=SUM + (1 - e^2/h^2)^2;
        end
    end
    npde(j+1,1)=th;
    npde(j+1,2)=0.9375*SUM/(n*h);
end
r=1;
npde(:,3) = r*(1 + pi*npde(:,2)).^2 - r;
% shift range to ±? if necessary 
if thrange==1
    npde(:,1)=mp2pp(npde(:,1));
    npde=sortrows(npde,1);
end
npde(:,4) = cumtrapz(r2d(npde(:,1)),d2r(npde(:,2)));
end

