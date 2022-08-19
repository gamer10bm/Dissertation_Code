%% Settings 
%Constellation Settings
const_vec = {'TX'};
nuangs = 0.000007*((1:numel(const_vec))-1); %Greatly affects cross-track and along-track
iangs = 0.00000000*((1:numel(const_vec))-1); %Mainly affects elevation with little cross-track effect
RAangs = -0.0000001*((1:numel(const_vec))-1);%zeros(size(const_vec));
wangs = -0.0000069*((1:numel(const_vec))-1);

%Orbit Settings
perigee_altitude = 555; %km
RAANinit = 100.0330; %deg From STK

%Set Simulation Start Date
yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 1;
min_init = 0; sec_init = 0;
init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];

%% Constellation Initialization (Auto)
num_sats = length(const_vec);
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