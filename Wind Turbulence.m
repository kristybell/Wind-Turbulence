% Data (u) are presented first - Matlab script may be found starting from line 7510
clear 
clear variables
close all
clc

%% Loading the data, finding and plotting the along-wind turbulence 
load('wind_data_part1.mat')

phi_average=mean(phi);      % mean wind direction [deg], determined from Matlab data

U=0.44704*V.*cosd(phi-phi_average); % Conversion from mph to m/s

Umean=mean(U);              % along-wind mean speed [m/s], determined from Matlab data
u=U-Umean;                  % Wind turbulence data u [m/s] - determined from Matlab input data file

Dt=1/25;                     % Data sampling - time step [seconds]
time1=Dt*(0:1:length(u)-1);  % time vector (resampled, just for verification)

figure(1)
Fts=18;
plot(time1,u)
title('\bfAlong-wind turbulence','Fontname','Times New Roman','Fontsize',Fts)
xlabel('Time [seconds]','Fontname','Times New Roman','Fontsize',Fts-2)
ylabel('Along-wind turbulence  \itu\rm [m/s]','Fontname','Times New Roman','Fontsize',Fts-2)
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Fts-2)


%% Calculating and plotting Autocovariance function of wind turbulence

Ru=xcov(u,'unbiased');
Varu=max(Ru);
Ruu=Ru/Varu;

% Time lag vector, tau
tau=[-(length(u)-1)*Dt:Dt:+(length(u)-1)*Dt];

figure(2)
plot(tau,Ruu)
title('\bfAutocovariance function','Fontname','Times New Roman','Fontsize',Fts)
xlabel('Time lag \it\tau \rm [seconds]','Fontname','Times New Roman','Fontsize',Fts-2)
ylabel('Auto-covariance function  \itR_u_u\rm/\it\sigma_u\rm^2','Fontname','Times New Roman','Fontsize',Fts-2)
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Fts-2)
legend('\itu \rm(from Matlab data)')
legend boxoff

%% Estimating integral scale of turbulence
npt=length(u);      % number of data points
idelta=250;         % max index of points within interval 0<tau<10 seconds

disp('Calculating Integral scale of turbulence (m):')
Lu=(sum(Ruu(npt:npt+idelta))+0.5*(Ruu(npt)+Ruu(npt+idelta)))*Dt*Umean

%% Find PSD(One-sided Spectrum) of the along-wind turbulence

nfftS=4096/4;   % number of Fourier Points (resolution)
Fs=25;          % samples/second

% Suu: ONE-SIDED PSD 
% n:   frequency, Hz
[Su,n]=periodogram(u,[],'onesided',nfftS,Fs); 

f=n*10/Umean;   % Monin or similarity coordinate
Sun=n.*Su/Varu; % Normalized spectrum: nSu/var(u)

% Plotting normalized PSD data on log-log axes
figure(3)
loglog(f,Sun)
hold on 

title('\bfAlong-wind turbulence spectrum (log-log axes)','Fontname','Times New Roman','Fontsize',Fts)
xlabel('Normalized Monin coordinate','Fontname','Times New Roman','Fontsize',Fts-2)
ylabel('One-sided PSD,  \itnS_u_u\rm/\it\sigma_u\rm^2','Fontname','Times New Roman','Fontsize',Fts-2)
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Fts-2)

%% Compare experimental PSD against a Kaimal-like spectrum

kaimal_original=250*f./(1+50*f).^(5/3);  % Original Kaimal model

% Finding the new proportionality constant (by least squares), noting that
% the Kaimal model can be written as Y=AX (A is the sought proportionality. constant)
% with X and Y as shown below

% Restrict data set analysis to the frequency interval 0<n<4Hz (due to instrumental resolution)
[ijx]=find(n<=4);

X=f(ijx)./(1+50*f(ijx)).^(5/3);     % Unit-magnitude Kaimal model (X variable)
Y=Sun(ijx);                         % Normalized PSD data (Y variable)

[P_values]=polyfit(X,Y,1);          % Least-squares, linear regression

disp('Finding a more suitable value for the proportionality constant (Kaimal model):')
A_prop=P_values(1)

kaimal_modified=A_prop*f./(1+50*f).^(5/3);   % modified Kaimal model

loglog(f,kaimal_original,'g--','linewidth',2)
loglog(f,kaimal_modified,'r-.','linewidth',2)

legend('\itu \rm(from data)','Kaimal model','Modified Kaimal model','location','best')





return % END of program

