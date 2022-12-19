% GEO422:Homework #6 Part 6
%
% This script addresses Problem 6 in HW Assignment #6. This script
% studies the power spectral density of the mysterious signal of yearly 
% observations contained in the data vector YEARLY.PLT from scratch.
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
f = 1/dt; %sampling frequency

%Computations
n = length(meas_fixed); %length of the signal
a = hann(n); %window function (chosen to be boxcar by default)
nfft = 2^nextpow2(n); %Number of points for FFT
FT = fft(meas_fixed.*a,nfft); %Fourier Transform of Signal
FT = FT(1:((nfft/2)+1)); %We only need 1 half of the FFT (chose 1st half)
PSD = (1/(n*f))*abs(FT).^2; %Power Spectral Density
int = f/(nfft)*(0:(nfft/2)); %Frequency Axis
period = log10(1./(int*31536000)); %Time (period) axis and converted back to years

%Plotting
%Power Spectral Density vs. Frequency
figure(2)
subplot(2,1,1)
plot(int, PSD)
title('Power vs. Frequency')
ylabel('PSD (FFT 1-Side)')
xlabel('Frequency (Hz)')

%Power Spectral Density
subplot(2,1,2)
plot(period, PSD)
title('Power vs. Period')
ylabel('PSD (FFT 1-Side)')
xlabel('Period (Years - Log Scale)')


