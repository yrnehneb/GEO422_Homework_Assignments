% GEO422:Homework #6 Part 1
%
% This script addresses Problem 1 in HW Assignment #6. This script
% generates a chosen yet interesting signal, and samples it using various
% frequency intervals (either below, at or above the Nyquist Frequency) to
% explore the effects of aliasing. Finally, this script plots thee results
% in the time domain and also plots the power spectral density using
% MATLAB's periodogram function.
%
% 
% Last modified by bshenry@princeton.edu, 12/17/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Preliminaries
freq = [1 3 5 7 9]; %create an array of frequencies for synthetic signal

%Create Independent Variable
syms freq1(t) freq2(t) freq3(t) freq4(t) freq5(t)

%Creating the Super Interesting Synthetic Signal!
freq1(t) = sin(2*pi*freq(1)*t);
freq2(t) = cos(2*pi*freq(2)*t);
freq3(t) = sin(2*pi*freq(3)*t);
freq4(t) = cos(2*pi*freq(4)*t);
freq5(t) = sin(2*pi*freq(5)*t);
TS = freq1 + freq2 + freq3 + freq4 + freq5; %Final Product! 

%Set Sampling Rates 
NQ_rate = 2*max(freq); %The Nyquist Rate
NQ_above = NQ_rate+10;  %Higher than Nyquist Rate
NQ_below = NQ_rate-10;  %Lower than Nyquist Rate

%Create Interval of Sampling (Same interval length, different number of
%samples based on different increments)
t_NQ = -1:1/NQ_rate:1;
t_Above = -1:1/NQ_above:1;
t_Below = -1:1/NQ_below:1;


%Plotting the Time Domain Signal
figure(1)
subplot(2,2,1)
%Plot the True Signal
t = -1:0.001:1; %Give it something to plot on
TS = double(TS(t)); %Plotting the True Signal "TS"
plot(t,TS,'Color','k')
title('True Signal in Time Domain')
ylabel('Signal')
xlabel('Time')
ylim([-4.75 4.75])

%Plot the Nyquist Sampling Rate
subplot(2,2,2)
S_NQ = double(freq1(t_NQ) + freq2(t_NQ) + freq3(t_NQ) + freq4(t_NQ) + freq5(t_NQ));
plot(t,TS,'Color','k')
hold on
plot(t_NQ,S_NQ,'.-','Color','r')
title('True Signal Sampled at Nyquist Rate in Time Domain')
ylabel('Signal')
xlabel('Time')
ylim([-4.75 4.75])
legend('True Signal', 'Sampled @ Nyquist Rate','Location','northeast')

%Plot the Above Nyquist Sampling Rate
subplot(2,2,3)
S_AboveNQ = double(freq1(t_Above) + freq2(t_Above) + freq3(t_Above) + ...
    freq4(t_Above) + freq5(t_Above));
plot(t,TS,'Color','k')
hold on
plot(t_Above,S_AboveNQ,'.-','Color','r')
title('True Signal Sampled at Above Nyquist Rate in Time Domain')
ylabel('Signal')
xlabel('Time')
ylim([-4.75 4.75])
legend('True Signal', 'Sampled @ Above Nyquist Rate','Location','northeast')

%Plot the Below Nyquist Sampling Rate
subplot(2,2,4)
S_BelowNQ = double(freq1(t_Below) + freq2(t_Below) + freq3(t_Below) + ...
    freq4(t_Below) + freq5(t_Below));
plot(t,TS,'Color','k')
hold on
plot(t_Below,S_BelowNQ,'.-','Color','r')
title('True Signal Sampled at Below Nyquist Rate in Time Domain')
ylabel('Signal')
xlabel('Time')
ylim([-4.75 4.75])
legend('True Signal', 'Sampled @ Below Nyquist Rate','Location','northeast')

%Find the Power Spectral Density of the 3 Signal Estimates

%Below Nyquist Frequency
figure(2)
subplot(3,1,1)
[PSD_Below,X_Below] = periodogram(S_BelowNQ,[],max(256,2^nextpow2(length(S_BelowNQ))),NQ_below,"onesided","psd");
plot(X_Below, PSD_Below)
title('Below Nyquist Frequency')
ylabel('Power Spectral Density (PSD)')
xlabel('frequency (Hz)')

%At Nyquist Frequency
subplot(3,1,2)
[PSD_NQ,X_NQ] = periodogram(S_NQ,[],max(256,2^nextpow2(length(S_NQ))),NQ_rate,"onesided","psd");
plot(X_NQ,PSD_NQ)
title('At Nyquist Frequency')
ylabel('Power Spectral Density (PSD)')
xlabel('frequency (Hz)')
xlim([min(X_Below) max(X_Below)])
ylim([0 1.5])


%Above Nyquist Frequency
subplot(3,1,3)
[PSD_Above, X_Above] = periodogram(S_AboveNQ,[],max(256,2^nextpow2(length(S_AboveNQ))),NQ_above,"onesided","psd");
plot(X_Above, PSD_Above)
title('Above Nyquist Frequency')
ylabel('Power Spectral Density (PSD)')
xlabel('frequency (Hz)')
xlim([min(X_Below) max(X_Below)])
ylim([0 1.5])
