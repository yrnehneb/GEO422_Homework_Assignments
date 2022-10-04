function GEO422_HW1_Q2(N)
% GEO422_HW1_Q2
%
% INPUT
%
% N     The number of successive convolutions
%
% Last modified by bshenry@princeton.edu, 10/04/2022

%Intial Input
a = [0 0 0 1000 2 40 7 0 170 97 32];
b = [0 1 0 2 0 3 0 4 0 5 0 6 0 7];

%Computation
[x,c] = again(a,b);

%Now do it again (successive convolutions)
for index=1:N-1
    [x,c] = again(a,c);
end

%Compare with the Gaussian of same mean and variance
mc = trapz(x,x.*c);
vc = trapz(x, (x-mc).^2.*c);
%Now make a pdf for comparison
cp = normpdf(x, mc, sqrt(vc));

%Graphics
plot(x,c);
hold on
plot(x, cp)
hold off

%Cosmetics
title(sprintf('Result after %i convolutions', N))
grid on

%Sub-function
function [x,c] = again(a,b)
%Computation
c = conv(a,b);

%Turn this into a proper distribution 
x = linspace(-4,4,length(c));

%Normalize (pdf integrates to one)
n = trapz(x,c); c = c/n;


