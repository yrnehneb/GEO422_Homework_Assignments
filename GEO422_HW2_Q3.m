function varargout=GEO422_HW2_Q3(mu,sig,M,N)
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
distrib = 'norm';

%INPUTS

%Expectation
%mu=10;
%Standard Deviation
%sig=2;

%Number of the experiment
%M = 300;
%Size of the Sample
%N = 200;

%Create M rows of experiments with N numbers each
zp = random(distrib,mu,sig, [M N]);

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

%Compute chi-squared for each experiment to end up with M values
for index=1:M
    %create histogram for each experiment
    [a,b] = histcounts(z(index,:),k);
    a = a/sum(a);
    for ind=1:k
        %Histogram Dimensions
        bin_width = b(ind+1)-b(ind);
        bin_height = a(ind);
        %Find area of each histogram bin
        exp_areas(ind) = bin_width*bin_height;
    end
    for i=1:k
    %Find areas of theoretical distribution
    the_areas(i) = cdf('Normal',b(i+1), 0,1)-cdf('Normal',b(i), 0,1);
    end
    %Compute chi-squared values for experiment
    for j=1:k
    x2(j) = ((exp_areas(j)-the_areas(j)).^2)./the_areas(j);
    end
    %Record chi-squared for each experiment
    x2_tot(index) = sum(x2);
end
%Level of Significance Test
LvlofSig = prctile(x2, 95); 
if x2(1)>LvlofSig
    disp('From Differing Parent Population');
else
    disp('From Same Parent Population (Prediction is True)');
end
    
%Create the true chi-squared distribution
hold on
dof = k-3;
x_chi2 = 0:.01:1;
chi2 = chi2pdf(x_chi2,dof);
div = x_chi2.*chi2;
chi2 = chi2./(div*M);
%PLOTTING
plot(x_chi2, chi2)
hist(x2)
%Make it Pretty
legend('Theoretical', 'Experiment/Measured')
title('Differences between the True Chi^2 Values & Measured x^2 values')
ylabel('# of occurances')
xlabel('x^2 values')
hold off







    




