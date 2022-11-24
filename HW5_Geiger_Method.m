function [M,X,PX]=geigerlive(tS,S,v,M0)
% [M,X,PX]=GEIGER(tS,S,v,M0)
% SS=GEIGER(NaN,S,v,M0)
%
% Linearized inversion for earthquake location OR
% Forward model computation for testing purposes
%
% INPUT:
%
% tS   Arrival times at the different stations (vector)
% S    [xS yS zS] Station coordinates (matrix)
%      These are the "data" and are "given"
% v    Velocity in medium assumed to be homogenous (km/s)
% M0   [x0 y0 z0 t0] Earthquake coordinates and origin times
%      This is the initial "guess"
%
% OUTPUT:
%
% M    [xE yE zE tE] Best-fit solution of earthquake parameters OR
% SS   A cell array, output only if input tS is NaN, with:
%      SS{1}=[xS yS zS tS]; % The station coordinates/times
%      SS{2}=[xE yE zE tE]; % The earthquake coordinates/time
%      SS{3}=v; % The velocity used
% X    Results for M of previous iterations than the last and best one 
% PX   Misfit criterion for the previous iterations
%
% Last modified by bshenry@princeton.edu, 11/24/22

%Personal Note: Always set up vectors as columns 

%Let's Get an Initial Idea of the Lay of the Land

%PlOTTING: Station Locations

%XY-Plane of Stations
figure(1)
plot(S(:,1), S(:,2), 'v')
title('Station Locations in the XY-Plane')
xlabel('x')
ylabel('y')

figure(2)
%XZ-Plane
subplot(2,1,1)
plot(S(:,1), S(:,3), 'v')
title('Station Locations in the XZ-Plane')
xlabel('x')
ylabel('z')

%YZ-Plane
subplot(2,1,2)
plot(S(:,2), S(:,3), 'v')
title('Station Locations in the YZ-Plane')
xlabel('y')
ylabel('z')

%Plot Arrival Times
figure(3)
plot(tS, '-o')
title('Arrival Times')
ylabel('Time (s)')
xlabel('Station #')

%The Good Stuff

% Misfit criterion, in seconds (Made Small to Ensure the Model will
% continue to improve)
Ptol = 1e-16;

% FORWARD MODEL to calculate station arrival times
t0 = forward(S,M0,size(S,1),v);

% Initial Guess
iter = 0;
M = M0;
t = t0;
P = misfit(t,tS);
dM = 1e9;

% Iterate until the misfit falls below the tolerance or until the
% improvement becomes marginal; location in km, time in s; same order then
while P > Ptol & sum(abs(dM))/4 > Ptol & iter < 6
    % Save location and the misfit 
    X(iter+1,:) = M;
    PX(iter+1) = P;
    disp(fprintf('Iteration %3.3i ; misfit %8.3e', iter, P));
    % Update Iteration
    iter = iter + 1;
    % Calculate Sensitivity Matrix
    G = sensitivity(S,M,size(S,1),v,t);
    % Calculate new model parameters 
    dM= inverse(G,t,tS)';
    M = M + dM;
    % Calculate New Travel Times 
    t = forward(S,M,size(S,1),v);
    % Calculate New Misft 
    P = misfit(t,tS);
end
PX (iter+1) = P;
disp(fprintf('Iteration %3.3i ; misfit %8.3e', iter, P));
disp(fprintf('Best Estimate = %i', M));

%Error/Uncertainty Analysis
resid = t - tS;
meanerror = repmat(mean(resid), size(S,1), 1);
%varerror = [repmat(var(resid), size(S,1), 1) -repmat(var(resid), size(S,1), 1)];
stderror = [repmat(std(resid), size(S,1), 1) -repmat(std(resid), size(S,1), 1)];
disp(fprintf('Mean Error = %i',meanerror));
disp(fprintf('Uncertainty = %i',stderror));

% Plot the Observed vs. Predicted Arrival Times
figure(4)
hold on
plot(tS, '-o')
plot(t, '-x')
title('Observed vs. Predicted Arrival Times')
ylabel('Time (s)')
xlabel('Station Number (Unranked)')
legend('predicted', 'observed')
hold off


%Plot the Residuals
figure(5)
hold on
plot(resid, '-o')
plot(meanerror, '--')
plot(stderror, ':', 'Color','magenta')
title('Residuals Plot')
xlabel('Station Number (Unranked)')
ylabel('Difference between Observed and Predicted Arrival Times (s)')
legend('residuals', 'mean', 'standard deviation')
hold off


%Plot the Estimated Earthquake Location

%XY-Plane
figure(6)
hold on
plot(S(:,1), S(:,2), 'v')
plot(M(1,1), M(1,2), 'o','Color','r')
title('Station Locations in the XY-Plane')
xlabel('x')
ylabel('y')
legend('Station Locations', 'Earthquake Location Estimate')
hold off

figure(7)
%XZ-Plane
subplot(2,1,1)
plot(S(:,1), S(:,3), 'v')
hold on 
plot(M(1,1), M(1,3), 'o','Color','r')
title('Station Locations in the XZ-Plane')
xlabel('x')
ylabel('z')
legend('Station Locations', 'Earthquake Location Estimate')
hold off

%YZ-Plane
subplot(2,1,2)
plot(S(:,2), S(:,3), 'v')
hold on
plot(M(1,2), M(1,3), 'o','Color','r')
title('Station Locations in the YZ-Plane')
xlabel('y')
ylabel('z')
legend('Station Locations', 'Earthquake Location Estimate')
hold off
end

% The Functions: Where the Magic Happens
% FORWARD MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=forward(S,M,N,v)
% S is the station location matrix [x y z]
% M is the event location matrix [x y z t]
% N is the number of stations
T = M(4)+ sqrt(sum((S-repmat(M(1:3),N,1)).^2,2))/v;
end

% MISFIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=misfit(t,tS)
% t is the station arrival time (observed)
% tS is the station arrival time (predicted)
% W is the weight matrix (default value: ones)

%W = eye(length(t));

F = (tS-t)'*(tS-t);
F = F/length(tS);

end
 
% SENSITIVTY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K=sensitivity(S,M,N,v,t)
% S is the station location matrix [x y z]
% M is the event location matrix [x y z t]
% N is the number of stations
% v is the velocity
% t is the station arrival time
K = [-(S-repmat(M(1:3),N,1))/v^2./repmat(t-M(4),1,3) ones(N,1)];
end

% INVERSE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I=inverse(G,t,tS)
I = pinv(G)*(tS-t);
end
