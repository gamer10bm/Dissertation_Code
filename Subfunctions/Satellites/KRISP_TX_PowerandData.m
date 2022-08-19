function [Sat_usepow, batt_pow_inst, data_stored, data_rate, torque_build, SubPowStruct, sat_mode] = ...
    KRISP_TX_PowerandData(SC,state, inst_torque, inst_genpow, gs_windows, radar_windows)
% [Sat_usepow, batt_pow_inst, data_stored, data_rate, torque_build, SubPowStruct] =
%                  KRISP_TX_PowerandData(SC,state, inst_torque, inst_genpow)
%
% The purpose of this function is to take the state structure, torque
% vector, and solar power generation vector to calculate power usage across
% all subsystems of KUbeSat1. 
%
% Inputs:     SC - Spacecraft structure generated in KUbeSat1_Simulator
%             state - A structure generated in KUbeSat1_Simulator that must
%                         have at least t, lat, and lon
%             inst_torque - A vector of instanteneous external torques (N-m)
%             inst_genpow - A vector of instantaneous generated power (W)
%             gs_windows - A structure with fields of Lat, Lon, and BW (deg, deg, deg)
%             radar_windows - A structure with fields of Lat, Lon, and Rad (deg, deg, km)
%
% Outputs: Sat_usepow - Vector of total instantaneous use power (W)
%                batt_pow_inst - Vector of current battery capacity (Wsec)
%                data_stored - Vector of current data stored (kb)
%                data_rate - Vector of current data rates (kb/s)
%                torque_build - Momentum buildup (N-m-s)
%                SubPowStruct - A structure of instant use powers with
%                                         fields for each subsystem (W)
%                sat_mode- A cell array of strings of the current
%                               satellite operational mode
%
% Created by Bailey Miller 3/2/2021, Adapted 1/12/2022


%Get number of iterations from the state structure
n = length(state);
%Get dt from state structure
dt = state(2).t - state(1).t;

%% Usage Settings
%Comms Values
Ant_pow = 0; Radio_pow = 7.4; %W
    %Establish Idle power value (Systems continuously on)
Comms_idle = 0; %Nothing runs when not operating

%Radar Values
Radar_pow = 36.3;%W
    %Establish Idle power value (Systems continuously on)
Radar_idle = 0; %Nothing runs when not operating

%Power Values
EPS_pow = 0.2; batt_heat = 0.8;%W
    %Establish Idle power value (Systems continuously on)
Pow_idle = EPS_pow;

%CDH Values
Proc_CDH = 0.2; %W
    %Establish idle power value (Systems continuously on)
CDH_idle= Proc_CDH;

%Data Values
IRM_health = 1; %kb (file compiled over orbit and only added before TX)
Xband_TX_dat = 150e3; %kb/s (X-Band)
UHFband_TX_dat = 15e3; %kb/s (UHF-Band)
TX_channels = struct('DataRate_kbps',Xband_TX_dat,'Band','X','Power_W',6);
TX_channels(end+1) = struct('DataRate_kbps',UHFband_TX_dat,'Band','UHF','Power_W',1.4);

Radar_dat = 44e3; %kb/s %43e3

%Battery Values
batt_fullcharge = 80*3600;%Wsec based on current battery full charge

%ADCS Values
Reaction_Wheels = 0.2; Mag_Torq = 0.75; %W
Star_Trackers = 1; Gyros = .1; Proc_ADCS = 1; %W
    %Establish Idle power value (Systems continuously on)
ADCS_idle = Star_Trackers + Gyros + Proc_ADCS + Reaction_Wheels;
    %Set constants for momentum buildup conditions
dump_torque = 3e-3; %N-m
max_momentum = 10e-3; %N-m-sec
min_momentum = 0; %N-m-sec (restricts the output momentum dump values)

%Low Power Mode Values
low_pow_enter_Wsec = 25*3600; %Whr Minimum value of battery charge to enter low power
low_pow_exit_Wsec = 78*3600; %Whr Max value of batter charge to exit low power

%Data Purge Mode Values
data_purg_enter_GB = 32; %GB Maximum storage value before executing
data_purg_exit_GB = .05; %GB Minimum storage value before exiting

%% Quick Calcs

steps2dump = ceil(max_momentum/(dump_torque*dt)); %Needed for checking telemetry and momentum dump overlap

%% Iterate over state
%Initialize the usepow vectors
Comms_usepow = zeros(1,n);
Radar_usepow = zeros(1,n);
ADCS_usepow = zeros(1,n);
Pow_usepow = zeros(1,n);
CDH_usepow = zeros(1,n);
Sat_usepow = zeros(1,n);
%Initialize torque vector
torque_build = zeros(1,n+1); %N-m
%Initialize battery power vector
batt_pow_inst = zeros(1,n+1);
batt_pow_inst(1) = batt_fullcharge; %Start with full battery
%Initialize data vectors
data_stored = zeros(1,n+1); %Array values are cummulative (kb)
data_rate = zeros(1,n); %Instant values of data collection (kb/s)
transmit_time = zeros(1,n+1); %Cummulative values of transmit time
%Initialize Mode Cell Array
sat_mode = cell(1,n);
%Initialize changing variables
minnextdump = n; maxRadarlength = 0; Radarlengthnow = 0; maxtxlength = 0; txlengthnow = 0;
low_pow_mode = false;  data_purg_mode = false; sat_mode_now = 'Nominal';
%Do any precalculations
    %Determine valid TX and GS combinations
    TX_channels_cell = struct2cell(TX_channels);
    TX_bands = TX_channels_cell(strcmp(fieldnames(TX_channels),'Band'),:);
    TX_pows = cell2mat(TX_channels_cell(strcmp(fieldnames(TX_channels),'Power_W'),:));
    TX_dats = cell2mat(TX_channels_cell(strcmp(fieldnames(TX_channels),'DataRate_kbps'),:));
    gs_windows_cell = struct2cell(gs_windows);
    gs_bands = gs_windows_cell(strcmp(fieldnames(gs_windows),'Bands'),:);
    bandcheck = zeros(length(TX_channels),length(gs_windows));
    for id = 1:length(gs_windows)
        for kd = 1:length(TX_channels)
            bandcheck(kd,id) = any(strcmp(TX_bands{kd},gs_bands{id}));
        end
    end
%Iterate
for i = 1:n
    lon = state(i).lon; lat = state(i).lat;
    alt_now = norm(state(i).R)-SC.Re; %km above Earth
    %% 1) Comms Check
    %Get geocircle for transmit
    telemcheck = zeros(1,length(gs_windows));
    for gs_id = 1:length(gs_windows)
        %Determine Band/s used
        if any(bandcheck(:,gs_id))
            gs_ant_BW = gs_windows(gs_id).BW;
            wind_radius = alt_now*tand(gs_ant_BW/2);
            lat_gs = gs_windows(gs_id).Lat;
            lon_gs = gs_windows(gs_id).Lon;

            %Check if within circle
            telemcheck(gs_id) = abs(geodist(lat_gs,lon_gs,lat,lon))<wind_radius;
        end
    end
    
    if any(telemcheck) && data_stored(i) > 0
        TX_pownow = sum(TX_pows(logical(bandcheck(:,logical(telemcheck)))));
        %Add TX_pow
        Comms_usepow(i) = Comms_idle + TX_pownow; %W
        txlengthnow = txlengthnow +1;
    else
        Comms_usepow(i) = Comms_idle;
        if txlengthnow > maxtxlength
            maxtxlength = txlengthnow;
        end
        txlengthnow = 0;
    end
    %% 2) Payload Check
    %Get geocircle for Radar pulses
    for wind_id = 1:length(radar_windows)
        wind_radius = radar_windows(wind_id).Rad;
        lat_wind = radar_windows(wind_id).Lat;
        lon_wind = radar_windows(wind_id).Lon;

        %Check if sat is over window
        Radarcheck(wind_id) = abs(geodist(lat_wind,lon_wind,lat,lon))<wind_radius;
    end

   if any(Radarcheck) && ~low_pow_mode && ~data_purg_mode
        Radar_usepow(i) = Radar_idle+Radar_pow;
        Radarlengthnow = Radarlengthnow+1;
   else
       Radar_usepow(i) = Radar_idle;
       if Radarlengthnow > maxRadarlength
           maxRadarlength = Radarlengthnow;
       end
       Radarlengthnow = 0;
   end  
    %% 3) ADCS Calcs
    if strcmp(sat_mode_now,'Momentum Dump')
      sat_mode_now = sat_mode_last;
    end
    %Fill in instantaneous buildup
    torque_build(i+1) = torque_build(i)+inst_torque(i)*dt;
    %Determine steps to next momentum dump
    nextdump = 0; torque_buildfut = torque_build(i+1);
    while torque_buildfut < max_momentum && i+nextdump+1 < n
        nextdump = nextdump +1;
        torque_buildfut = torque_buildfut+inst_torque(i+nextdump)*dt;
    end
    if torque_build(i) <= min_momentum && nextdump < minnextdump
        minnextdump = nextdump; %Gives idea of minimum steps to fill momentum buildup
    end
    
    %Determine if momentum dump is needed now
    momentum_dumpnow_chk = false;
    if nextdump <= 2
        %Do momentum dump now otherwise reaction wheels will be overloaded
        momentum_dumpnow_chk = true;
    end
        
    %Allocate usage power based on momentum dumping now
    if momentum_dumpnow_chk
        %Do momentum dump now
        ADCS_usepow(i) = ADCS_idle + Mag_Torq - Reaction_Wheels; %MAY NEED TO CHANGE
        torque_build(i+1) = torque_build(i+1) - dump_torque*dt;
        %Restrict torque_buildup after momentum dump (for now set to zero)
        if torque_build(i+1) < min_momentum
            torque_build(i+1) = min_momentum;
        end
        sat_mode_last = sat_mode_now;
        sat_mode_now = 'Momentum Dump';
    else
        %Nominal operation
        ADCS_usepow(i) = ADCS_idle;
    end
    %% 4) Power Calcs
    %Check if in eclipse
    if inst_genpow(i)<=0
        Pow_usepow(i) = Pow_idle + batt_heat;
    else
        Pow_usepow(i) = Pow_idle;
    end
    %% 5) CDH Calcs
    CDH_usepow(i) = CDH_idle;
    
    %% 6) Data Rate Calcs
    %Determine the current data rate based on what is operating
    if Comms_usepow(i) ~= Comms_idle
        %Currently transmitting
        data_rate(i) = -sum(TX_dats(logical(bandcheck(:,logical(telemcheck)))));
    elseif Radar_usepow(i) ~= Radar_idle
        %Currently operating Radar
        data_rate(i) = Radar_dat;
    end
    
    %Determine data storage at the next time step
    data_stored(i+1) = data_stored(i)+data_rate(i)*dt;
    if data_stored(i+1) < 0
        %All data has been transferred 
        data_stored(i+1) = 0;
    end
    if data_rate(i) < 0
        %Currently transmitting
        transmit_time(i+1) = transmit_time(i) + dt;
    else
        transmit_time(i+1) = transmit_time(i);
    end
    
    %% 7) Battery Calcs
    Sat_usepow(i) = Comms_usepow(i) +Radar_usepow(i) + ADCS_usepow(i)+CDH_usepow(i)+Pow_usepow(i); %W
    inst_netpow = Sat_usepow(i)-inst_genpow(i);
    %Determine battery power at this instant
    batt_pow_inst(i+1) = batt_pow_inst(i) - inst_netpow*dt; %Wsec
    if batt_pow_inst(i+1) > batt_fullcharge
        batt_pow_inst(i+1) = batt_fullcharge;
    end
    
    %% 8) Low Power Check
    if batt_pow_inst(i+1)<low_pow_enter_Wsec && ~low_pow_mode
        %Enter low power mode
        low_pow_mode = true;
        sat_mode_now = 'Low Power';
    elseif batt_pow_inst(i+1)>low_pow_exit_Wsec && low_pow_mode
        %Exit low power mode
        low_pow_mode = false;
        sat_mode_now = 'Nominal';
    end
   
    %% 9) Data Purge Check (overriden by low power mode)
    if data_stored(i+1) >= data_purg_enter_GB*8*1000*1000 && ~data_purg_mode && ~low_pow_mode
        %Enter data purge mode
        data_purg_mode = true;
        sat_mode_now = 'Data Purge';
    elseif data_stored(i+1) <= data_purg_exit_GB*8*1000*1000 && data_purg_mode
        %Exit data purge mode
        data_purg_mode = false;
        sat_mode_now = 'Nominal';
    end
    
    %% Save outputs
    sat_mode{i} = sat_mode_now;

end
%% Send to outputs
%Format subsystem output structure
SubPowStruct = struct('Sat',Sat_usepow,'ADCS',ADCS_usepow,...
    'CDH',CDH_usepow,'Comms',Comms_usepow,'Power',Pow_usepow,...
    'Radar',Radar_usepow);