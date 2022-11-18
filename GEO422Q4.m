%function varargout=GEO422_HW2_Q3(mu,sig,M,N)
%Function creates X^2 values given GEO422_HW2_Q3(mu,sig,M,N)
%
%INPUT:
%
%mu         Population Expectation
%sig        Population Standard Deviation
%M          Number of Identical Experiments 
%N          Number of Samples in Every Experiment
%
%OUTPUT:
%
%X^2        M observed sums of squares of N experimentally generated 
%           variables which theory says should be distributed as 
%           chi-squared with N degrees of freedom
%
%Last modified by bshenry@princeton.edu, 10/15/2022

%Preliminaries
distrib = 'unif';

%INPUTS

%Expectation
mu=10;
%Standard Deviation
sig=2;

%Number of the experiment
M = 300;
%Size of the Sample
N = 200;

%Create M rows of experiments with N numbers each
zp = random(distrib,sig,mu, [M N]);

%Standardize the variables zp->z
z = (zp-mu)./sig;

%Choose the number of Bins using Sturges' Rule
Y = 3/2*ceil(log2(N));
%K bins = rounded Y
k = round(Y);

%Create an array for areas of histogram, theoretical, and chi-squared
exp_areas = zeros(1,k);
the_areas = zeros(1,k);
x2= zeros(1,k);

%Array to store all the chi-squared values from M experiments
x2_tot = zeros(1, M);
F_i = zeros(1,k);



%Compute chi-squared for each experiment to end up with M values
for index=1:M
    %create histogram for each experiment
    [a,b] = histcounts(z(index,:),k);
    %Calculate the frequencies of the experiment
    f_i = a/sum(a);
    %Create the real(parent) pdf distribution
    GaussPdf = normpdf(b,0,1);
    %Calculate the frequencies of the True parent distribution
    Gauss = normcdf(b,0,1);
    for x=1:k
        F_i(x) = Gauss(x+1)-Gauss(x);
    end
    %Compute chi-squared values for experiment
    for j=1:k
    x2(j) = (F_i(j)-f_i(j)).^2/F_i(j);
    if j==k
         %Record chi-squared for each experiment
         x2_tot(index) = sum(x2);
    end
    end
   
end

%Calculate the p-values with calculation x^2 values
p_values = zeros(1,M);
dof = k-3;
for g =1:M
    p_values(g) = chi2cdf(x2_tot(g), dof);
end

%PLOTTING
% Three of the random variables - PICKED AT RANDOM
for index=1:fmax-1
  % The output is the figure "handle"
  ah(index)=subplot(2,2,index);
  % Make the histogram of the standardized variables
  % Pick three random rows
  rr=randi(M,fmax-1);
  % Think about the bin size, using a Sturges like rule
  [a,b]=hist(z(rr(index,1),:),Y);
  % Plot the histogram
  bh(index)=bar(b,a,1);
  t(index)=title(sprintf('m = %i out of M = %i, N = %i',...
  rr(index),M,N));
  % Fake it a bit
  hold on; f=plot(-3:1:3, normpdf(-3:1:3),'w.'); hold off
  leg(index)=legend(f,...
  sprintf('x2 = %5.1f',x2_tot(rr(index))),'Location', ...
  'NorthEast');
  set(leg(index),'EdgeColor','w')
  % Prepare for limit calculation
  axis tight
  yls(index,:)=ylim;
  xls(index,:)=xlim;
  xlabel('standard normal variates')
  ylabel('frequency')
end

%Significance Test
signif = 1 - mean(p_values);
if signif > 0.05
    fprintf('Greater than level of significance, accept hypothesis')
elseif signif <= 0.05
    fprintf('Less than level of significance, reject hypothesis')
end


%Create the true chi-squared distribution
figure(2)
hold on
dof = k-3;
x_chi2 = 0:.01:40;
chi2 = chi2pdf(x_chi2,dof);
%PLOTTING
plot(x_chi2, chi2)
%Make it Pretty
title('True Chi^2 Distribution')
ylabel('probability density')
xlabel('Chi^2 values')
hold off

figure(3)
hist(x2_tot, Y)
title('Distribution of x^2 values')
ylabel('probabilty')
xlabel('x^2 values')
%NOTE: I'm aware that there is some graphical error, and I can't figure out
%how to make the correct plot. I hope my methods in calculating p are more
%correct than my last attempt. Thank you Srijan for all of your help.
    
