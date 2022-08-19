% This script is used to generate validation tables using the developed
  % propagator
  
clc
clearvars

addpath(genpath('Subfunctions'))
cmdsize = matlab.desktop.commandwindow.size;
cmdline = repmat('=',1,cmdsize(1));
%% Settings
%Simulation Settings
N =2; %Number of orbits
timesteps = 206; %Time steps per orbit (206 has the best results using Optimize_StepSize2VerifySTK_2Body.m)

%Orbit Settings
perigee_altitude = 555; %km
RAANinit = 100.0330; %deg From STK

%Set Simulation Start Date
yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 1;
min_init = 0; sec_init = 0;
init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];

%Propagator of interest
proptype = '2-Body'; %'2-Body' or 'J2+Drag'
%Drag solar settings
F10abs = 78;
F10avg = 78.71;
Apvec = [3 3 3 2 3 3 4];
%% Calculations
%Initialize Orbit
OE0 = SSO_Earth(perigee_altitude);
OE0(4) = (RAANinit)*pi/180; % Sets RAAN and centers groundtrack on Lawrence

%Grab driving orbital elements
hmag = OE0(1); emag = OE0(2); iang = OE0(3); nuang = OE0(6); RAang = OE0(4); wang = OE0(5);

%Initialize satellite model
SCnow = initiate_RX_model;

%Calculate position and velocity vectors    
[~,Rnow,Vnow] = coe2RV(hmag,emag,iang,nuang,RAang,wang,SCnow.mu);
COEstruct0 = RV2coe(Rnow,Vnow,SCnow.mu);

%Define time array
period_sec = COEstruct0.T_Period;
num_states = timesteps*N; % Basing array on number of orbits and points per orbit
time = linspace(0, N*period_sec, num_states); % define time array
dt = time(2) - time(1); % save time step
tstart = 0; tend = N*period_sec;

%% Propagate Orbit with Runga-Kutta
%Grab correct derivative function
switch proptype
  case '2-Body'
    dervfunc = @(t,X)OrbitDerivFunc_2Body(X,SCnow.Re,SCnow.mu,SCnow.J2,SCnow.CD,...
        SCnow.A,SCnow.m,SCnow.rho0,SCnow.r0,SCnow.H,SCnow.thetadot);
  case 'J2+Drag'
    dervfunc = @(t,X)OrbitDerivFunc_J2Drag(X,SCnow.Re,SCnow.mu,SCnow.J2,SCnow.CD,...
      SCnow.A,SCnow.m,SCnow.rho0,SCnow.r0,SCnow.H,SCnow.thetadot,yr_init,...
      day(datetime(init_utcvec),'dayofyear'),t,F10avg,F10abs,Apvec);
end
fprintf('Running %s Propagation with Runga-Kutta\n',proptype)
%Do the Runge-Kutta
Xstart = [Rnow;Vnow];
[tvec,XRK]=RungeKutta(dervfunc,Xstart,dt,tstart,tend+dt);

%Convert XRK from km to m
XRK = XRK.*1000;
%% Print Results
fprintf('\n%s\n%s\nPrinting Results for %s Propagation\n%s\n%s\n',...
  cmdline,cmdline,proptype,cmdline,cmdline)
fprintf('\nt(sec)\tX(m)\tY(m)\tZ(m)\tVx(m/s)\tVy(m/s)\tVz(m/s)\n%s\n\n',cmdline)
for tid = 1:length(tvec)
  %Grab print components
  tnow = tvec(tid); xnow = XRK(1,tid); ynow = XRK(2,tid); znow = XRK(3,tid); vxnow = XRK(4,tid); vynow = XRK(5,tid); vznow = XRK(6,tid);
  fprintf('%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n',tnow,xnow,ynow,znow,vxnow,vynow,vznow);
end