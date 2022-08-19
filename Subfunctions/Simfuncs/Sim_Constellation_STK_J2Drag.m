function [const_struct] = Sim_Constellation_STK_J2Drag(const_struct,init_utcvec,N,timesteps)
%Grab settings
Settings_STK

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
stepsize = time(2) - time(1); % save time step
tstart = 0; tend = N*period_sec;

%Calculate days to end simulation
daystoend = days(seconds(tend));

%Set Simulation Date Strings
dt_init = datetime(init_utcvec);
timestartstr = datestr(dt_init,'dd mmmm yyyy HH:MM:SS.FFF');
dt_end = dt_init+days(daystoend);
timeendstr = datestr(dt_end,'dd mmmm yyyy HH:MM:SS.FFF');
%% Initialize STK
init_STK
%% Simulate Constellation
satnames = cell(1,num_sats);
for satid = 1:num_sats
    %Grab satellite details
    satname = const_struct(satid).name;
    SCnow = const_struct(satid).SC;
    Rnow = const_struct(satid).Rinit;
    Vnow = const_struct(satid).Vinit;
    COE0 = const_struct(satid).COE0;
    
    %% Create and Propagate Satellite Object
    createprop_sat
    %% Export propagation data
    exportdata_STK
    %% Load Outputs    
    const_struct(satid).states = sat_states;
end

