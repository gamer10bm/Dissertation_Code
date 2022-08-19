%% Settings 
%Constellation Settings
const_vec = {'RX','RX','RX','TX','RX','RX','RX','RX'};
if 0 %The final ideal point
    iangs = 2.556e-06*((1:numel(const_vec))-1);
    RAangs = 1.000e-06*((1:numel(const_vec))-1);
    wangs = 0*((1:numel(const_vec))-1);
    nuangs = 2.222e-08*((1:numel(const_vec))-1);
elseif 0 %The Broad Ground track Plot Opt point
    iangs = 3e-06*((1:numel(const_vec))-1);
    RAangs = 1.1e-06*((1:numel(const_vec))-1);
    wangs = 0*((1:numel(const_vec))-1);
    nuangs = 0*((1:numel(const_vec))-1);
elseif 1 %The Broad Ground track Plot Bad point
    iangs = 12e-06*((1:numel(const_vec))-1);
    RAangs = 1.1e-06*((1:numel(const_vec))-1);
    wangs = 0*((1:numel(const_vec))-1);
    nuangs = 0*((1:numel(const_vec))-1);
elseif 0 %Old and busted
    nuangs = 0.000007*((1:numel(const_vec))-1); %Greatly affects cross-track and along-track
    iangs = 0.00000000*((1:numel(const_vec))-1); %Mainly affects elevation with little cross-track effect
    RAangs = -0.0000001*((1:numel(const_vec))-1);%zeros(size(const_vec));
    wangs = -0.0000069*((1:numel(const_vec))-1);
elseif 0 %From Ideal_Array_Optimize.m %old
    iangs = 4.00e-06*((1:numel(const_vec))-1);
    RAangs = 0.0000000*((1:numel(const_vec))-1);
    wangs = .000002*2*((1:numel(const_vec))-1);
    nuangs = 0.0000000*((1:numel(const_vec))-1);
elseif 1 %From Ideal_Array_Optimize.m %After cluster runs 8/17
    iangs = 2e-5*((1:numel(const_vec))-1);
    RAangs = 0*((1:numel(const_vec))-1);
    wangs = 0*((1:numel(const_vec))-1);
    nuangs = 0*((1:numel(const_vec))-1);
elseif 1 %From Ideal_Array_Optimize.m %After cluster runs 7/8
    iangs = 0.0000001*6*((1:numel(const_vec))-1);
    RAangs = 0*((1:numel(const_vec))-1);
    wangs = 0e-6*((1:numel(const_vec))-1);
    nuangs = 0*((1:numel(const_vec))-1);
end
    

%Orbit Settings
perigee_altitude = 555; %km
RAANinit = 100.0330; %deg From STK

%Set Simulation Start Date
yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 1;
min_init = 0; sec_init = 0;
init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];

%Set allowables for failure
fail_min_off = 0.25; %m
fail_prob = 0.95; %random value that fail check must exceed

%% Constellation Initialization (Auto)
num_sats = length(const_vec);
satids = 1:num_sats;
driving_sat_id = satids(strcmp(const_vec,'TX'));
%Initialize Orbit
OE0 = SSO_Earth(perigee_altitude);
OE0(4) = (RAANinit)*pi/180; % Sets RAAN and centers groundtrack on Lawrence
OE0(6) = 1e-6+pi/2; %Sets True anomaly to non-zero
OE0(2) = 1e-6;
%Grab driving orbital elements
hmag = OE0(1); emag = OE0(2); iang = OE0(3); nuang = OE0(6); RAang = OE0(4); wang = OE0(5);%-pi/100;
%Initialize output
sats(num_sats) = struct('SC',struct(),'states',struct(),'name','','Rinit',[],'Vinit',[],'COE0',struct());

%Iterate over defined constellation
RXnum = 0; TXnum = 0;
for satid = 1:num_sats
    %% Satellite Details
    %Grab constellation details
    satnow_str = const_vec{satid};
    %Name satellite
    if strcmp(satnow_str,'RX')
        satname = sprintf('%s%d',satnow_str,RXnum);
        RXnum = RXnum+1;
        %Add RX satellite
    elseif strcmp(satnow_str,'TX')
        satname = sprintf('%s%d',satnow_str,TXnum);
        TXnum = TXnum+1;
        %Add TX satellite
    else
        satname = satnow_str;
        %Add rando sat
    end
    
    %Initialize satellite model
    SC_model_fh = str2func(sprintf('initiate_%s_model',satnow_str));
    SCnow = SC_model_fh();
    
    %% Calculate Orbit Adjustment
    %Take direct changes to orbital elements
    nuoff = nuangs(satid); ioff = iangs(satid); RAoff = RAangs(satid); woff = wangs(satid);
    
    %% Propagate Orbit
    %Calculate position and velocity vectors    
    [~,Rnow,Vnow] = coe2RV(hmag,emag,iang+ioff,nuang+nuoff,RAang+RAoff,wang+woff,SCnow.mu);
    COEstruct0 = RV2coe(Rnow,Vnow,SCnow.mu);
    
    %% Load Outputs
    sats(satid).name = satname;
    sats(satid).SC = SCnow;
    sats(satid).Rinit = Rnow;
    sats(satid).Vinit = Vnow;
    sats(satid).COE0 = COEstruct0;
end

%Print orbit offsets
if 1
    num_sats = length(sats);
    fprintf('\nSatellite Name\t')
    for i = 1:num_sats
        fprintf('%s\t',sats(i).name);
    end
    fprintf('\nTrue Anomaly (deg.)\t')
    for i = 1:num_sats
        nu_dif = sats(1).COE0.nu_TrueAnomaly-sats(i).COE0.nu_TrueAnomaly;
        fprintf('%.1e\t',rad2deg(nu_dif));
    end
    fprintf('\nRight Ascension(deg.)\t')
    for i = 1:num_sats
        RA_dif = sats(1).COE0.RA_RightAscension-sats(i).COE0.RA_RightAscension;
        fprintf('%.1e\t',rad2deg(RA_dif));
    end
    fprintf('\nArg. of Perigee(deg.)\t')
    for i = 1:num_sats
        w_dif = sats(1).COE0.w_ArgumentofPerigee-sats(i).COE0.w_ArgumentofPerigee;
        fprintf('%.1e\t',rad2deg(w_dif));
    end
    fprintf('\nInclination(deg.)\t')
    for i = 1:num_sats
        i_dif = sats(1).COE0.i_Inclination-sats(i).COE0.i_Inclination;
        fprintf('%.1e\t',rad2deg(i_dif));
    end
    fprintf('\n')
end