function [q] = VMquantile(P,K,mu)
%function [q] = VMquantile(P,K,mu)
%Quantile function for the von Mises distribution - for a cumulative
%probability, P, the quantile angle is output
% K=5;
I0K=besseli(0,K);
f = 0.5;
t = 0;
c = log(I0K);
epsilon = 1e-4;
d=10;
while d > epsilon
    g = f - P;
    d = exp(log(abs(g)) + c - K*cos(t));
    t = t - sign(g)*d;
    f = VMcumfreq(t,K,mu);
%     v = abs(f-P);
%     disp(v)
end
q=t;

