% GEO422:Homework #6 Part 2
%
% This script addresses Problem 2 in HW Assignment #6. This script
% studies the power spectral density of the mysterious signal of yearly 
% observations contained in the data vector YEARLY.PLT using periodogram.
%
% 
% Last modified by bshenry@princeton.edu, 12/17/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read in the YEARLY.PLT Data 
data = importdata('yearly.txt');

%Create Variables 
years = data(:,1); %years
meas = data(:,2);  %measurements

%Get a sense of the data 
figure(1)
plot(years, meas)
title('Raw Data')
ylabel('Measurement')
xlabel('Time (years)')

%Normalize the Data 
meas_fixed = (meas - mean(meas))/std(meas);
dt = 31536000; %time interval
f = 1/dt; %sampling frequency/rate

%Computations
[PSD,X] = periodogram(meas_fixed,[],max(256,2^nextpow2(length(meas_fixed))),f,"onesided","psd");
Period = log10(1./(X*31536000)); %Converts Frequency (Hz) into Years

%Plotting
figure(2)
subplot(2,1,1)
plot(X,PSD)
title('PSD vs Frequency')
ylabel('PSD')
xlabel('Frequency (Hz)')
subplot(2,1,2)
plot(Period,PSD)
title('PSD vs Period')
ylabel('PSD')
xlabel('Period (Years - Log10)')
