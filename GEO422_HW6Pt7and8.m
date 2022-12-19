% GEO422:Homework #6 Part 7 & 8
%
% This script addresses Problem 7 & 8 in HW Assignment #6. This script
% targets a set of specific frequencies and fit sines and cosines of those 
% target frequencies in a chosen signal using the tools of inverse theory,
% and assesses misfits
%
% NOTE: I'm essentially using the code discussed in class by Frederik,
% except now I'm using my own synthetic signal. (dec6.m)
%
% Last modified by bshenry@princeton.edu, 12/17/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A time column
t=linspace(0,100,100)';

% The true period
T=12;

freq = [1 3 5 7 9]; %create an array of frequencies for synthetic signal

%Choose/Create the Signal (Using My Synthetic Signal From Part 1)
% A deterministic "harmonic" periodic signal - a column
s = sin(2*pi*freq(1)*t/T) + cos(2*pi*freq(2)*t/T) + ...
    sin(2*pi*freq(3)*t/T) + cos(2*pi*freq(4)*t/T) + ...
    sin(2*pi*freq(5)*t/T);

% A noise-to-signal fraction
f=0.5;
% Add some random noise
s=s+f*randn(size(s))*std(s);

% Pretend you did NOT know the period. I used a set of specific periods
% (just frequencies if you invert it)
Tt=1:1:100;

for i = 1:length(Tt)

    % Now try to find the coefficients at the known period T
    gmat=[sin(2*pi*freq(1)*t/Tt(i)) cos(2*pi*freq(2)*t/Tt(i)) ...
        sin(2*pi*freq(3)*t/Tt(i)) cos(2*pi*freq(4)*t/Tt(i)) ...
        sin(2*pi*freq(5)*t/Tt(i))];

    % Now find the coefficients mhat that estimate m=[a b c]';
    mhat=pinv(gmat)*s;
    % Plot the signal prediction
    shat=gmat*mhat;
    % Asses the misfit - in percent
    phi(i)=100*var(shat-s)/var(s);
    % Asses the power
    Shat(i)=sqrt(mhat(1)^2+mhat(2)^2+mhat(3)^2+mhat(4)^2+mhat(5)^2);
    S(i)=sqrt(5*1^2);
end

%Plotting
figure(1)
subplot(3,1,1)
plot(t,s); grid on
title('Original Data with Noise')
ylabel('Signal')
xlabel('Time (indexed)')

subplot(3,1,2)
plot(1./Tt,Shat)
title('Power vs. Period')
ylabel('Power')
xlabel('Frequency')
ylim([0 max(Shat)*1.2])

subplot(3,1,3)
plot(1./Tt,phi)
title('Assessing Misfit')
ylabel('Phi')
xlabel('Frequency')
ylim([0 max(phi)*1.2])


%Now try Periodogram
figure(3)
[PSD,X] = periodogram(s,[],1./Tt,max(256,2^nextpow2(length(s))),"onesided","power");
plot(X, PSD)
title('Power Spectral Density of Synthetic Signal')
ylabel('PSD')
xlabel('Frequency (Hz)')



