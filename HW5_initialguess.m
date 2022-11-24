function [M0] = initialguess(tS,S)
% [M0]=INITIALGUESS(tS,S)
%
% Function Used to Create an Initial Guess to be used in HW5_Geiger_Method
%
% INPUT:
%
% tS   Arrival times at the different stations (vector)
% S    [xS yS zS] Station coordinates (matrix)
%      These are the "data" and are "given"
% 
%
% OUTPUT:
%
% M0   [x0 y0 z0 t0] Earthquake coordinates and origin times
% This is the initial "guess"
%
% Last modified by bshenry@princeton.edu, 11/21/22

% The Good Stuff 

% Create one cohesive dataset then sort it by min to max arrival times
Data = [S(:,1) S(:,2), S(:,3) tS];
Data = sortrows(Data, 4); 

%Plotting: Getting a Lay of the Land
figure(1)
subplot(2,2,1)
plot(Data(:,4),Data(:,1))
title('X-Coordinates vs. Time')
xlabel('time(s)')
ylabel('Distance in X-Direction [L]')

subplot(2,2,2)
plot(Data(:,4),Data(:,2))
title('Y-Coordinates vs. Time')
xlabel('time(s)')
ylabel('Distance in Y-Direction [L]')

subplot(2,2,3)
plot(Data(:,4),Data(:,3))
title('Z-Coordinates vs. Time')
xlabel('time(s)')
ylabel('Distance in Z-Direction [L]')

subplot(2,2,4)
plot(Data(:,4))
title('Timeseries of Arrival Times')
xlabel('Rank = Min to Max Measurements')
ylabel('Time (s)')

%Coming Up With Average Slopes
xslope = (Data(end,1)-Data(1,1))/(Data(end,4)-Data(1,4));
yslope = (Data(end,2)-Data(1,2))/(Data(end,4)-Data(1,4));
zslope = (Data(end,3)-Data(1,3))/(Data(end,4)-Data(1,4));
disp(fprintf('xslope = %3.3i ' , xslope));
disp(fprintf('yslope = %3.3i' , yslope));
disp(fprintf('zslope = %3.3i' , zslope));

%Scaling Back from the 1st station to find M0
M0 = [Data(1,1)-xslope*Data(1,4) Data(1,2)-yslope*Data(1,4) ...
    Data(1,3)-zslope*Data(1,4) 0];
%disp(fprintf('M0 = %3.3i', M0));





