% GEO422:Homework #6 Part 4
%
% This script addresses Problem 4 in HW Assignment #6. This script
% explores different window functions and their periodograms, and how they 
% provide spectral estimates of the input YEARLY.PLT.
%
% 
% Last modified by bshenry@princeton.edu, 12/17/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read in the YEARLY.PLT Data 
data = importdata('yearly.txt');

%Create Variables 
years = data(:,1); %years
meas = data(:,2);  %measurements

%Normalize the Data 
meas_fixed = (meas - mean(meas))/std(meas);
dt = 31536000; %time interval (s)
f = 1/dt; %sampling frequency
n = length(meas_fixed);

%Plot Normalized Data
figure(1)
plot(years, meas_fixed)
title('Normalized Data')
ylabel('Measurements Relative to the Mean Value')
xlabel('Times (Years)')

%Create The Window Functions 
BL = bartlett(n);
CH = chebwin(n);
HA = hann(n);

%Plot The Window Functions
wvtool(BL)
wvtool(CH)
wvtool(HA)

%Create the Modified Periodofunctions
%Frequency Axis
figure(2)
[PSD_BL,X_BL] = periodogram(meas_fixed,BL,max(256,2^nextpow2(length(meas_fixed))),f,"onesided","psd");
plot(X_BL,PSD_BL)
hold on
[PSD_CH,X_CH] = periodogram(meas_fixed,CH,max(256,2^nextpow2(length(meas_fixed))),f,"onesided","psd");
plot(X_CH,PSD_CH)
hold on
[PSD_HA,X_HA] = periodogram(meas_fixed,HA,max(256,2^nextpow2(length(meas_fixed))),f,"onesided","psd");
plot(X_HA,PSD_HA)
title('Power Spectral Density Estimate w/ Different Window Functions')
ylabel('PSD')
xlabel('Frequency')
legend('Bartlett','Chebwin','Hann')

