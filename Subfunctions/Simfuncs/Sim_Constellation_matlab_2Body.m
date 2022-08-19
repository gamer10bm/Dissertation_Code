function [const_struct] = Sim_Constellation_matlab_2Body(const_struct,init_utcvec,N,timesteps)

%Parse inputs
num_sats = length(const_struct);
yr_init = init_utcvec(1);
mnth_init = init_utcvec(2);
day_init = init_utcvec(3);
hr_init = init_utcvec(4);
min_init = init_utcvec(5);
sec_init = init_utcvec(6);
period_sec = const_struct(1).COE0.T_Period;

%% Precalculations
%Define time array
num_states = timesteps*N; % Basing array on number of orbits and points per orbit
time = linspace(0, N*period_sec, num_states); % define time array
dt = time(2) - time(1); % save time step
tstart = 0; tend = N*period_sec;

%% Iteration
%Iterate over defined constellation
for satid = 1:num_sats
    %% Satellite Details
    %Grab constellation details
    satname = const_struct(satid).name;
    SCnow = const_struct(satid).SC;
    Rnow = const_struct(satid).Rinit;
    Vnow = const_struct(satid).Vinit;
    COE0 = const_struct(satid).COE0;
    
    %% Propagate Orbit
    %Do Runga-Kutta
    Xstart = [Rnow;Vnow];

    dervfunc = @(t,X)OrbitDerivFunc_2Body(X,SCnow.Re,SCnow.mu,SCnow.J2,SCnow.CD,...
        SCnow.A,SCnow.m,SCnow.rho0,SCnow.r0,SCnow.H,SCnow.thetadot);
    %Do the Runge-Kutta
    [~,XRK]=RungeKutta(dervfunc,Xstart,dt,tstart,tend+dt);

    %Initialize outputs
    sat_states(num_states) = struct('t',0,'OE',COE0,'R',[0;0;0],'V',[0;0;0],'lat',0,'lon',0,'alt',0);
    %Parse results and do calculations
    for i = 1:num_states
        %Declare time
        sat_states(i).t = time(i);
        %Put RK results into state matrix
        Rnow = XRK(1:3,i); Vnow = XRK(4:6,i);
        sat_states(i).R = Rnow;
        sat_states(i).V = Vnow;
        %Calculate orbital elements
        sat_states(i).OE = RV2coe(Rnow,Vnow,SCnow.mu);
        %Calculate lat and lon
        [sat_states(i).lat, sat_states(i).lon, sat_states(i).alt] =ECEF2latlon(Rnow,sat_states(i).t);
    end
    
    %% Send to outputs
    const_struct(satid).states = sat_states;
end
end