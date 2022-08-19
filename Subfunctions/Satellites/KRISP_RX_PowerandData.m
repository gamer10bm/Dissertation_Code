function [Sat_usepow, batt_pow_inst, data_stored, data_rate, torque_build, SubPowStruct] = ...
    KRISP_TX_PowerandData(SC,state, inst_torque, inst_genpow)
% [Sat_usepow, batt_pow_inst, data_stored, data_rate, torque_build, SubPowStruct] =
%                  KRISP_TX_PowerandData(peralt,state, inst_torque, inst_genpow)
%
% The purpose of this function is to take the state structure, torque
% vector, and solar power generation vector to calculate power usage across
% all subsystems of KUbeSat1. 
%
% Inputs: SC - Spacecraft structure generated in KUbeSat1_Simulator
%             state - A structure generated in KUbeSat1_Simulator that must
%                         have at least t, lat, and lon
%             inst_torque - A vector of instanteneous external torques
%                                    (N-m)
%             inst_genpow - A vector of instantaneous generated power (W)
%
% Outputs: Sat_usepow - Vector of total instantaneous use power (W)
%                batt_pow_inst - Vector of current battery capacity (Wsec)
%                data_stored - Vector of current data stored (kb)
%                data_rate - Vector of current data rates (kb/s)
%                torque_build - Momentum buildup (N-m-s)
%                SubPowStruct - A structure of instant use powers with
%                                         fields for each subsystem (W)
%
% Created by Bailey Miller 3/2/2021

%Set lat and lon for Ground station and beamwidth
lat_gs = 38.971669; %(deg) [Lawrence, KS]
lon_gs = -95.23525; %(deg) [Lawrence, KS]
gs_ant_BW = 131; %Groundstation Antenna Beamwidth (deg)

% Set lat and lon for McMurdo and beamwidth
lat_mcmurdo = -82.406; lon_mcmurdo = 0; %(deg) McMurdo Station
mcmurdo_ant_BW = 90;

%Get number of iterations from the state structure
n = length(state);
%Get dt from state structure
dt = state(2).t - state(1).t;
%% Comms Power
Comms_usepow = zeros(1,n);
%Set power usage values per component<<<<<<<<<<<<<<<=======================
Ant_pow = 0; Radio_pow = 2; %W
%Establish Idle power value (Systems continuously on)
Comms_idle = 0; %Nothing runs when not operating
%Determine power based on each time step and given conditions
maxtxlength = 0; txlengthnow = 0;
for i = 1:n
    lon = state(i).lon; lat = state(i).lat;
    alt_now = norm(state(i).R)-SC.Re; %km above Earth
    %Get geocircle for transmit
    capture_radius = alt_now*tand(gs_ant_BW/2);
    [telemcirc_lats, telemcirc_lons] = geocircle(lat_gs,lon_gs,capture_radius);
    %Check if within telemetry circle
    lontelemcheck = ~(all((lon-telemcirc_lons)>0) || all((lon-telemcirc_lons)<0));
    lattelemcheck = ~(all((lat-telemcirc_lats)>0) || all((lat-telemcirc_lats)<0));

    if lontelemcheck && lattelemcheck
        %Add TX_pow
        Comms_usepow(i) = Comms_idle + Ant_pow + Radio_pow; %W
        txlengthnow = txlengthnow +1;
    else
        Comms_usepow(i) = Comms_idle;
        if txlengthnow > maxtxlength
            maxtxlength = txlengthnow;
        end
        txlengthnow = 0;
    end
end

%% HiCal-Lite Power
HiCal_usepow = zeros(1,n);
%Set power usage values per component
HiCal_pow = 0.5;%W
%Establish Idle power value (Systems continuously on)
HiCal_idle = 0; %Nothing runs when not operating
%Determine power based on each time step and given conditions
maxHiCallength = 0; HiCallengthnow = 0;
for i = 1:n
    lon = state(i).lon; lat = state(i).lat;
    alt_now = norm(state(i).R)-SC.Re; %km above Earth
    %Get geocircle for HiCalX pulses
    capture_radius = alt_now*tand(mcmurdo_ant_BW/2);
    [HiCalXcirc_lats, HiCalXcirc_lons] = geocircle(lat_mcmurdo,lon_mcmurdo,capture_radius);

    %Check if sat is over antarctica 
    latHiCalXcheck = ~(all((lat-HiCalXcirc_lats)>0) || all((lat-HiCalXcirc_lats)<0));
    lonHiCalXcheck = ~(all((lon-HiCalXcirc_lons)>0) || all((lon-HiCalXcirc_lons)<0));

   if latHiCalXcheck && lonHiCalXcheck
        HiCal_usepow(i) = HiCal_idle+HiCal_pow;
        HiCallengthnow = HiCallengthnow+1;
   else
       HiCal_usepow(i) = HiCal_idle;
       if HiCallengthnow > maxHiCallength
           maxHiCallength = HiCallengthnow;
       end
       HiCallengthnow = 0;
   end           
end
%% ADCS Power and Momentum Buildup
ADCS_usepow = zeros(1,n);
%Set power usage values per component
Reaction_Wheels = 0.2; Mag_Torq = 0.75; %W
Star_Trackers = 1; Gyros = .1; Proc_ADCS = 1; %W
%Establish Idle power value (Systems continuously on)
ADCS_idle = Star_Trackers + Gyros + Proc_ADCS + Reaction_Wheels;
%Set constants for momentum buildup conditions
dump_torque = 3e-3; %N-m
max_momentum = 10e-3; %N-m-sec
min_momentum = 0; %N-m-sec (restricts the output momentum dump values)
steps2dump = ceil(max_momentum/(dump_torque*dt)); %Needed for checking telemetry and momentum dump overlap
%Determine power based on each time step and given conditions
torque_build = zeros(1,n+1); %N-m
minnextdump = n;
for i = 1:n
    %Fill in instantaneous buildup
    torque_build(i+1) = torque_build(i)+inst_torque(i)*dt;
    %Determine number of steps to transmit
    nexttx = find(Comms_usepow(i:end)>Comms_idle,1,'first');
    if isempty(nexttx)
        nexttx = -maxtxlength; %No more transmit operations left
    end
    %Determine number of steps to HiCal operation
    nextHiCal = find(HiCal_usepow(i:end)>HiCal_idle,1,'first');
    if isempty(nextHiCal)
        nextHiCal = -maxHiCallength; %No more HiCal operations left
    end
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
    if nextdump+steps2dump >= nexttx && nextdump <= nexttx + maxtxlength && nexttx <= steps2dump +1 && nexttx >0
        %Do momentum dump in anticipation of upcoming transmit
        momentum_dumpnow_chk = true;
    elseif nextdump+steps2dump >= nextHiCal && nextdump <= nextHiCal + maxHiCallength && nextHiCal <=steps2dump+1 && nextHiCal >0
        %Do momentum dump in anticipation of upcoming HiCal operation
        momentum_dumpnow_chk = true;
    elseif nextdump == 0
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
    else
        %Nominal operation
        ADCS_usepow(i) = ADCS_idle;
    end
end
%% Power System Power
Pow_usepow = zeros(1,n);
%Set power usage values per component
EPS_pow = 0.2; batt_heat = 0.8;%W
%Establish Idle power value (Systems continuously on)
Pow_idle = EPS_pow;
%Determine power based on each time step and given conditions
for i = 1:n
    %Check if in eclipse
    if inst_genpow(i)<=0
        Pow_usepow(i) = Pow_idle + batt_heat;
    else
        Pow_usepow(i) = Pow_idle;
    end
end

%% PCRD Power
PCRD_usepow = zeros(1,n); 
%Set power usage values per component
PCRD_pow = 4; %W
swap_window = 145; %sec time required to switch operations
%Establish Idle power value (Systems continuously on)
PCRD_idle = 0;
%Determine power based on each time step and given conditions
for i = 1:n
    %PCRD operates when:
        %No transmitting
        %No momentum dump
        %No HiCal
    swap_ind_width = ceil(swap_window/dt);
    ind_strt = i-swap_ind_width;
    if ind_strt <=0
        ind_strt = 1;
    end
    ind_end = i+swap_ind_width;
    if ind_end >= n
        ind_end = n;
    end
    swap_inds = ind_strt:ind_end;
    %Check all other operations in swap window
    if any(Comms_usepow(swap_inds) ~=Comms_idle) || any(HiCal_usepow(swap_inds) ~= HiCal_idle) ...
            || any(ADCS_usepow(swap_inds) ~= ADCS_idle)
        %PCRD not operating
        PCRD_usepow(i) = PCRD_idle;
    else
        %PCRD is operating
        PCRD_usepow(i) = PCRD_idle+PCRD_pow;
    end
end

%% CDH Power
CDH_usepow = zeros(1,n);
%Set power usage values per component
Proc_CDH = 0.2; %W
%Establish idle power value (Systems continuously on)
CDH_idle= Proc_CDH;
%Determine power based on each time step and given conditions
CDH_usepow(1:n) = CDH_idle; %Because CDH should always be running at full power

%% Sum all subsytem values for instant usage power
Sat_usepow = Comms_usepow + HiCal_usepow + ADCS_usepow +...
    PCRD_usepow + CDH_usepow + Pow_usepow; %W
%Format subsystem output structure
SubPowStruct = struct('KUbeSat1',Sat_usepow,'ADCS',ADCS_usepow,...
    'CDH',CDH_usepow,'Comms',Comms_usepow,'Power',Pow_usepow,...
    'HiCalK',HiCal_usepow,'PCRD',PCRD_usepow);
%% Use instant power and generated power to determine battery power storage
%Determine the instanteous power with solar power generation
inst_netpow = Sat_usepow-inst_genpow;
%Determine the battery power storage with time
batt_pow_inst = zeros(1,n+1);
batt_fullcharge = 40*3600;%Wsec based on current battery full charge
batt_pow_inst(1) = batt_fullcharge; %Start with full battery
for i = 1:n
    %Determine battery power at this instant
    batt_pow_inst(i+1) = batt_pow_inst(i) - inst_netpow(i)*dt; %Wsec
    if batt_pow_inst(i+1) > batt_fullcharge
        batt_pow_inst(i+1) = batt_fullcharge;
    end
end
%% Data stored based on usage power
%Iterate through timesteps to determine data storage and instant rate
data_stored = zeros(1,n+1); %Array values are cummulative (kb)
data_rate = zeros(1,n); %Instant values of data collection (kb/s)

IRM_health = 1; %kb (file compiled over orbit and only added before TX)

TX_dat = -250; %kb/s 

HiCal_dat = 40e-3; %kb per pulse (per discord pulse 1 Hz)

PCRD_dat = 2; %kb/s (per discord conversation)

transmit_time = zeros(1,n+1); %Cummulative values of transmit time

for i = 1:n
    %Determine the current data rate based on what is operating
    if Comms_usepow(i) ~= Comms_idle
        %Currently transmitting
        data_rate(i) = TX_dat;
    elseif HiCal_usepow(i) ~= HiCal_idle
        %Currently operating HiCal
        data_rate(i) = HiCal_dat;
    elseif PCRD_usepow(i) ~= PCRD_idle
        %Currently operating PCRD
        data_rate(i) = PCRD_dat;
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
end
end