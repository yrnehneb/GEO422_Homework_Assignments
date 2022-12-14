%GEO422:Homework #4
%
%
%
% 
% Last modified by bshenry@princeton.edu, 11/14/2022


%Preliminaries
frame = 25;  %fps

%Picking images in which ball is in motion only
i = 33:1:54; 
for index = 1:22
    url = strcat('http://geoweb.princeton.edu/people/simons/GOLFBALL/000000', num2str(i(index)), '.jpg');
    b = imread(url);
    image(b)
    [x(index,:),y(index,:)] = ginput(1); 
end

%Pixel to Meter Conversion
d = 0.04; %ball diameter in m
%Find diameter of ball in pixels
for index = 1:2
    b = imread(strcat('http://geoweb.princeton.edu/people/simons/GOLFBALL/000000', num2str(i(1)), '.jpg'));
    image(b)
    [garbage(index,:),edge(index,:)] = ginput(1); 
end
conv = d/(edge(1) - edge(2)); %pixels to meter

%convert height array
h = y*conv;
h = -1*(h - max(h));
for index = 1:length(i)
    if index == 1
        h(index) = 0;
    else
        h(index) = (h(index) - h(index-1));
        %add change in time to positively increment
        h(index) = h(index-1) + h(index);
     end
end

%convert time array
t = 1:length(i);
t = t/frame;

%Model Assumption/Constraints
% h(t) = m1 + m2*t + -1/2*m3*t^2

%G Matrix
G = zeros(length(i),3);
G(:,1) = 1;
G(:,2) = t;
G(:,3) = (-1/2)*t.^2; 

% Estimated Model Parameter Vector
G_g = inv(G.'* G)*G.';
% Find Model 1 Parameters
m1(:,1) = G_g*h;
% Estimate for gravity
a1 = m1(3); %Finding Acceleration under this model to be around 9.2 m/s^2

%Assign Unncertainties to the Measurements: Assuming the uncertainties
%arise from the difference in the data (h) from the modeled data (h_hat)

%Create Covarianace Matrix: Cov = (d - E(d))*(d - E(d))' 
C_y = (h-G*m1(:,1))*(h-G*m1(:,1))';
%Create Model Covariance Matrix 
C_m = G_g*C_y*G_g'./length(t);

%New Model Using Covariance Weight
G_gnew = inv(G.'*inv(C_y)*G)*G.'*inv(C_y);
mnew(:,1) = G_gnew*h;
%New Estimate For Gravity
a2 = mnew(3); 

%Not sure why my a1 and a2 values are the same, but will move forward
%anyway

%Chi-Squared Test 
dof = length(i)-3; %We have length(i) measurements and 3 parameters
h_hat(:,1) = G*mnew(:,1);   %predicted data based on model
data_error = h - h_hat(:,1); %residuals
sig2 = var(data_error);
x2 = zeros(length(i));
for index = 1:length(i)
    x2(index) = data_error(index).^2./sig2;
end
Chi2 = sum(x2); %Chi-Squared Value
incr = 0:0.01:100;
%Plot the Chi2 results
figure(2)
hold on
plot(incr, chi2pdf(incr, dof))
title('Chi^2 Plot (dof = 19)')

%Confidance Interval
c = chi2cdf(Chi2(1), dof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOTTING%

%Raw Converted Data Points
figure(3)
plot(t,h,'.','LineWidth',2.0)
title('Height vs. Time Plot')
ylabel('Height (meters)')
xlabel('Time (seconds)')
ylim([0 0.85])

%Plottting the Data Based on Model Parameters
figure(4)
plot(t,h,'x','Color','k')
hold on
plot(t,G*m1(:,1),'-.','Color','r','LineWidth',1.2)
plot(t,h_hat,'Color','b')
title('Height vs. Time')
ylabel('Height (meters)')
xlabel('Time (seconds)')
legend('Raw Data Points', 'old model (m1)', 'new model (mnew)')
ylim([0 0.85])

%Plotting h and h_hat
figure(5)
plot(t,h,'x','Color','k')
hold on 
plot(t, h_hat,'x','Color','r')
title('Height vs. Time')
ylabel('Height (meters)')
xlabel('Time (seconds)')
legend('Measured (y)', 'Predicted (y-hat)')
ylim([0 0.85])









