clc
clearvars

%This script is the amalgam of all of my orbital programming experience to
%date
% [~,hostname] = system('hostname');
% hostname(regexp(hostname,'[\n]')) = [];
% switch hostname
%     case 'MillerXPS' %Personal laptop
%         base_dir = 'C:\Users\Bailey\Documents\GitHub\KUbeSat_Simulations';
%     case 'MillerPC1' %Personal Desktop
%         base_dir = 'C:\Users\bmiller\OneDrive - jayhawksupportitservices\Documents\Coding\MATLAB\KUbeSat_Simulations';
%     otherwise %KU Desktop
%         base_dir = 'E:\Documents\GitHub\KUbeSat_Simulations';
% end
% addpath(fullfile(base_dir,'KUbeSat_Subfunctions'))
% addpath(fullfile(base_dir,'3D_Shape'))
addpath(genpath('Subfunctions'))
%% Settings
%Define settings for simulation
perigee_altitude = 555; %km
N =20; %Number of orbits
timesteps = 300; %Time steps per orbit
yawsteps = 100; %Yaw steps per time step
RAANinit = 100.0330; %deg From STK
% RAANinit = 360+lon_gs-97; %deg First pass is over lawrence

%%%%%%%% Select Sat type %%%%%%%%%%%
sattype = 'TX'; %'RX';

%Set lat and lon for Ground stations
gs_windows = struct('Lat',38.971669,'Lon',-95.23525,'BW',131,'Bands',{{'UHF','S','X'}}); %Lawrence, KS
if 1
    % SSC Space Ground Stations for X-Band Telemetry
    %https://sscspace.com/ssc-worldwide/ground-station-network/
    gs_windows(end+1) = struct('Lat',19+1/60,'Lon',-155-40/60,'BW',131,'Bands',{{'S','X'}}); %South Point
    gs_windows(end+1) = struct('Lat',64+48/60,'Lon',-147-39/60,'BW',131,'Bands',{{'S','X'}}); %North Pole
%     gs_windows(end+1) = struct('Lat',68+24/60,'Lon',-133-30/60,'BW',131,'Bands',{{'S','X'}}); %Inuvik
%     gs_windows(end+1) = struct('Lat',26+44/60,'Lon',-81-2/60,'BW',131,'Bands',{{'S','X'}}); %Clewiston
%     gs_windows(end+1) = struct('Lat',-33-8/60,'Lon',-70-40/60,'BW',131,'Bands',{{'S','C','Ka'}}); %Santiago
    gs_windows(end+1) = struct('Lat',-52-56/60,'Lon',-70-51/60,'BW',131,'Bands',{{'S','X'}}); %Punta Arenas
%     gs_windows(end+1) = struct('Lat',67+53/60,'Lon',21+4/60,'BW',131,'Bands',{{'UHF','S','X'}}); %Esrange Space Center
%     gs_windows(end+1) = struct('Lat',13+6/60,'Lon',100+55/60,'BW',131,'Bands',{{'S','X'}}); %Siracha
%     gs_windows(end+1) = struct('Lat',-29-5/60,'Lon',115+35/60,'BW',131,'Bands',{{'S','X','Ku'}}); %WASC
    
elseif 0
    %Other Options    
    gs_windows(end+1) = struct('Lat',-33.8696,'Lon',151.20695,'BW',131); %Sydney, Australia
    gs_windows(end+1) = struct('Lat',48.86,'Lon',2.349,'BW',131); %Paris, France
    gs_windows(end+1) = struct('Lat',64.8366,'Lon',-147.74,'BW',131); %Fairbanks, Alaska
end

% Set lat and lon for radar operations
radar_windows = struct('Lat',-87,'Lon',90,'Rad',2500); %Antartica
radar_windows(end+1) = struct('Lat',71,'Lon',-43,'Rad',1300); %Greenland
%% Orbit Sim from KUbeSat1
if 1
    % initialize TX 6U Satellite Model
%     SC = initiate_TX6U_model();
    SC = initiate_TX_model();

    % initialize simulation epoch
    yr_init = 2022; mnth_init = 1; day_init = 1; hr_init = 0;
    min_init = 0; sec_init = 0;
    init_utcvec = [yr_init, mnth_init, day_init, hr_init, min_init, sec_init];
    epoch = initial_epoch(yr_init, mnth_init, day_init, hr_init, min_init, sec_init);

    % determine SSO orbit for given altitude and spacecraft
    OE0 = SSO_Earth(perigee_altitude,SC);
    OE0(4) = (RAANinit)*pi/180; % Sets RAAN and centers groundtrack on Lawrence
    [~,Rstart,Vstart] = coe2RV(OE0(1),OE0(2),OE0(3),OE0(6),OE0(4),OE0(5),SC.mu);
    COEstruct0 = RV2coe(Rstart,Vstart,SC.mu);

    % compute orbit period and initial mean anomaly
    period_sec = COEstruct0.T_Period;

    % determine time array parameters
    n = timesteps*N; % Basing array on number of orbits and points per orbit
    time = linspace(0, N*period_sec, n); % define time array
    dt = time(2) - time(1); % save time step

    % initialize vehicle state structure
     if exist('state','var')
        clearvars state
    end
    state(n) = initiate_SC_state(SC);

    % preallocate states structure memory
    for i = n-1:-1:1
        state(i) = state(n);
    end

    %% Do Keplerian Motion Simulation
    COEstructnow = COEstruct0;
    for i = 1:n
        % store time from array
        state(i).t = time(i);

        % find and store vehicle translation states    
        [~,COEstructnow] = FutureAnomaly(state(i).t,COEstructnow);
        state(i).OE = COEstructnow;
        [~,R, V] = coe2RV(COEstructnow);
        state(i).R = R; state(i).V = V;

        % current routines do not use a separate estimation state but save
        % it anyway in case we do in the future
        state(i).R_est = R; state(i).V_est = V;
    end

    %% Do Runga-Kutta Simulation with J2 and Drag
    %Initial conditions
    rstart =state(1).R; %km
    vstart = state(1).V; %km/s

    Xstart = [rstart;vstart];

    dervfunc = @(t,X)OrbitDerivFunc(X,t,SC.Re,SC.mu,SC.J2,SC.CD,SC.A,SC.m,...
      SC.rho0,SC.r0,SC.H,SC.thetadot,year(datetime(init_utcvec)),day(datetime(init_utcvec),'dayofyear'));
    %Do the Runge-Kutta
    tstart = 0; tend = N*period_sec;
    [tRK,XRK]=RungeKutta(dervfunc,Xstart,dt,tstart,tend+dt);

    %Put RK results into state matrix
    for i = 1:n
        Rnow = XRK(1:3,i); Vnow = XRK(4:6,i);
        state(i).R = Rnow;
        state(i).V = Vnow;
    end

    for i = 1:n
          %Determine latitude and longitude (no earth motion)
        R_IJK = state(i).R;
        [state(i).lat, state(i).lon] =ECEF2latlon(R_IJK,state(i).t);
    end
end

%% Satellite Simulation from KUbeSat1
if 1
    % identify the yaw angle array (eliminate redundant boundary)
    n2 = yawsteps;
    Yaw = linspace(0, 2*pi, n2+1); Yaw = Yaw(1:n2);
    %% Do solar panel generation 
    %Iterate over generated states to determine power generation from solar
    %panels
    inst_genpow = zeros(1,n);
    % We need to set the reference value for the quaternion conversion
    % routine. In the initial step, set this to identity [0,0,0,1].
    % Otherwise, set this to the value that produced the optimum last time.
    Q_near = [0; 0; 0; 1];
    for i = 1:n
        % find sun direction at this time
        [RS_au,RS_km] = sun_vector_JD(init_utcvec,state(i).t); %AU

        % initialize maximum solar power
        SP_max = 0;

        % cycle through possible yaw angles
        for j = 1:n2
            % determine applicable quaternion attitude based on yaw angle from
            % array and set pitch and roll to zero. This is based on the local
            % RTN coordinate system, so we need position and velocity vectors
            %
            % This routine takes in a reference quaternion to pick the branch
            % since our possible answer are non-unique. Using the most recent
            % value for quaternion attitude keeps this continuous.
            Q = quaternionC(0, 0, Yaw(j), state(i).R, state(i).V, Q_near);

            % predict solar collection at this yaw angle
%             SP = solar_panel(SC, state(i).R, Q, RS_au);
            SP = solar_illum(state(i).R, RS_km, SC.Re, Q, SC);

            % check if maximum and save
            if (SP > SP_max)
                SP_max = SP ;
                Qmax = Q;
            end
        end

        % save maximum power quaternion to state array attitude command
        if (SP_max > eps) 
            state(i).Q_com = Qmax; 
        else
            state(i).Q_com = Q_near;
        end

        % set maximum value to reference quaternion if no eclipse
        if (SP_max > eps) 
            Q_near = Qmax;
        end

        % save maximum power generated
        inst_genpow(i) = SP_max; %W
    end

    %% Momentum buildup
    inst_torque = zeros(1,n);
    for i = 1:n
        R = state(i).R; V = state(i).V;
        eclip_chk = false; %Default in sunlight
        if inst_genpow(i) <= 0
            eclip_chk = true; %Change to in eclipse
        end
%         inst_torque(i) = momentumbuildup(SC,R,V,eclip_chk); %N-m
        inst_torque(i) = momentum_Daniel(SC,R,V,epoch, state(i).t); %N-m
    end

    %% Power usage and data stored
    switch sattype
        case 'TX'
            [inst_usepow, batt_pow, data_stored, data_rate, torque_build, SubPowStruct, sat_mode] = ...
                KRISP_TX_PowerandData(SC,state,inst_torque,inst_genpow,gs_windows,radar_windows);
        case 'RX'
            [inst_usepow, batt_pow, data_stored, data_rate, torque_build, SubPowStruct, sat_mode] = ...
                KRISP_RX_PowerandData(SC,state,inst_torque,inst_genpow,gs_windows,radar_windows);
    end
    
    %% Print outputs
    %Get and print duty cycles from SubPowstruct
    fnames = fieldnames(SubPowStruct);
    fduty = zeros(1,length(fnames)-1);
    for fid = 2:length(fnames)
        %Get the maximum value for the field
        fmax = max(SubPowStruct.(fnames{fid}));
        %Get number of values at max
        fmax_occurs = length(find(SubPowStruct.(fnames{fid})==fmax));
        %Get duty cycle based on total SubPow length
        fduty(fid) = fmax_occurs./length(SubPowStruct.(fnames{fid}));
        %Print it
        fprintf('\n%s duty cycle = %.2f%s\n',fnames{fid},fduty(fid)*100,'%')
    end
    %Get total collected data
    collect_rate = data_rate;
    collect_logvec = or(collect_rate<0,data_stored(2:end)<=0);
    collect_rate(collect_logvec) = 0;
    collect_data = abs(cumsum(collect_rate.*[0 diff(time)]));
    fprintf('\n%.2f GB of Data Collected in %.0f Orbits\n',collect_data(end)/8/1000/1000,state(end).t/period_sec)
    %Get total telemetered data
    telem_rate = data_rate;
    telem_logvec = or(telem_rate>0,data_stored(2:end)<=0);
    telem_rate(telem_logvec) = 0;
    telem_data = abs(cumsum(telem_rate.*[0 diff(time)]));
    %Correct telemetered data
    i = 1;
    while i<length(data_stored)
        i = i+1;
        if data_stored(i) <= 0
            %Adjust telemetered data to match collected
            data_diff = telem_data(i-1) - collect_data(i-1);
            if data_diff ~= 0
                telem_data(i:end) = telem_data(i:end)-data_diff;
            end
        end
    end
    fprintf('\n%.2f GB of Data Telemetered in %.0f Orbits\n',telem_data(end)/8/1000/1000,state(end).t/period_sec)
    
    %% Plot Section
    save_figs = struct('fig',[],'savepath',[]);
    savedir = fullfile('Output_Figures',sattype);
    if 1
        %Convert time to fraction of orbits
        plot_time = time./period_sec;
        shorttimelim = [0 2];
        longtimelim = [0 state(end).t./period_sec];
        %%%%%%%%%%% Instant Power vs. Time %%%%%%%%%%%%%
        plotname = 'InstantPower_Time';
        figure(3)
        gp = plot(plot_time,inst_genpow);
        hold on
        plot([plot_time(1) plot_time(end)],[mean(inst_genpow) mean(inst_genpow)],':','Color',gp.Color)
        ip = plot(plot_time,inst_usepow);
        plot([plot_time(1) plot_time(end)],[mean(inst_usepow) mean(inst_usepow)],':','Color',ip.Color)
        hold off
        grid on
        legend({'Generated','Avg. Generated','Used','Avg. Used'})
        title('Instantaneous Power')
        xlabel('Orbits (~)')
        xlim(shorttimelim)
        ylabel('Power (W)')
        FontWidthandPos
        if 1
            save_figs(end+1) = struct('fig',gcf,'savepath',fullfile(savedir,[plotname '.fig']));
        end

        %%%%%%%% Subsystem Power vs. Time %%%%%%%%%%%%%
        oneperlim = [0 4];
        subfields = fieldnames(SubPowStruct);
        % subfields = subfields(2:end);
        plotname = 'SubsystemPower_Time';
        figure(10)
        for ids = 1:length(subfields)
            plot(plot_time,SubPowStruct.(subfields{ids}),'-.')
            hold on
        end
        hold off
        grid on
        legend(subfields)
        title('Instant Power by Subsystem')
        xlabel('Orbits (~)')
        ylabel('Power (W)')
        xlim(oneperlim)
        FontWidthandPos
        if 1
            save_figs(end+1) = struct('fig',gcf,'savepath',fullfile(savedir,[plotname '.fig']));
        end
        

        %%%%%%%%%%% Battery Storage vs. Time %%%%%%%%%%
        plotname = 'BatteryStore_Time';
        figure(4)
        plot(plot_time,batt_pow(1:end-1)./3600) %Whr
        grid on
        title('Battery Power')
        xlabel('Orbits (~)')
        xlim(longtimelim)
        ylabel('Power (Whr)')
        FontWidthandPos
        if 1
            save_figs(end+1) = struct('fig',gcf,'savepath',fullfile(savedir,[plotname '.fig']));
        end

        %%%%%%%%%% Data Storage vs. Time %%%%%%%%%%%
        plotname = 'DataStore_Time';
        figure(2)
        plot(plot_time,data_stored(1:end-1)./8/1000/1000)
        hold on
        plot(plot_time,telem_data./8/1000/1000)
        plot(plot_time,collect_data./8/1000/1000)
        hold off
        grid on
        title('Data Collection')
        legend({'On Satellite','Telemetered','Collected'})
        ylabel('Data (GB)')
        xlabel('Orbits (~)')
        xlim(longtimelim)
        FontWidthandPos
        if 1
            save_figs(end+1) = struct('fig',gcf,'savepath',fullfile(savedir,[plotname '.fig']));
        end

        %%%%%%%%%% Momentum Buildup vs. Time %%%%%%%%%
        plotname = 'MomentumBuildup_Time';
        figure(5)
        plot(plot_time,torque_build(1:end-1)./1e3)
        grid on
        title('Reaction Wheel Momentum Buildup')
        xlabel('Orbits (~)')
        xlim(longtimelim)
        ylabel('Momentum (N-m-sec)')
        FontWidthandPos
        if 1
            save_figs(end+1) = struct('fig',gcf,'savepath',fullfile(savedir,[plotname '.fig']));
        end
        
        %%%%%%%%%% Operational Mode vs. Time %%%%%%%%%%
        plotname = 'OperateMode_Time';
        %Get unique modes and indexing for plotting
        [unq_modes, ~, plot_modes] = unique(sat_mode);
        figure(8)
        plot(plot_time,plot_modes)
        grid on
        title('Operational Modes vs. Time')
        xlabel('Orbits (~)')
        xlim(longtimelim)
        ylabel('Mode')
        ylim([0 length(unq_modes)+1])
        set(gca,'YTick',[0, 1:length(unq_modes)+1] ,'YTickLabel', [{''} unq_modes {''}])
        FontWidthandPos
        if 1
            save_figs(end+1) = struct('fig',gcf,'savepath',fullfile(savedir,[plotname '.fig']));
        end

        %%%%%%%%%%%%% Ground Track %%%%%%%%%%%%%%%%
        statecell = struct2cell(state); statefields = fieldnames(state);
        latname = 'lat';
        latid = find(strcmp(statefields,latname));
        lonid = latid+1;
        lats = reshape(cell2mat(statecell(latid,:,:)),[1 n]);
        lons = reshape(cell2mat(statecell(lonid,:,:)),[1 n]);
        %Plot the ground track
        N_off = 1; %Orbit offsets in the groundtrack plot
        N_now = 1;
        minichunk = 1:timesteps;
        chunk = [];
        while N_now <=N
            if N_now == 1 || rem(N_now,N_off)==0
                chunk = [chunk minichunk+(N_now-1)*timesteps];
            end
            N_now = N_now+1;
        end
        figure(1)
        itvec = 1:10:length(chunk);
        if itvec(end) ~= length(chunk)
            itvec(end+1) = length(chunk);
        end
        for i =  length(chunk)%1:30:length(chunk)
%             keyboard
            clf
            subchunk = chunk(1:i);
            geoaxes();
            geoscatter(lats(subchunk),lons(subchunk),'.b')
            hold on
            %plot the GS windows
            for wid = 1:length(gs_windows)
                capture_radius = perigee_altitude*tand(gs_windows(wid).BW/2);
                [telemcirc_lats, telemcirc_lons] = geocircle(gs_windows(wid).Lat,gs_windows(wid).Lon,capture_radius);
                geoscatter(telemcirc_lats,telemcirc_lons,'.r')
            end
            %Plot the radar windows
            for wid = 1:length(radar_windows)
                [radcirc_lats, radcirc_lons] = geocircle(radar_windows(wid).Lat,radar_windows(wid).Lon,radar_windows(wid).Rad);
                geoscatter(radcirc_lats,radcirc_lons,'.m')
            end
            hold off
            geolimits([-90 90], [-180 180])
            title(sprintf('Ground Track: t = %.1f min.',time(subchunk(end))/60))
            % title(sprintf('Ground Track for %.0f Orbits',N))
            FontWidthandPos
            drawnow
        end
    end
    %Save figures
    if 0
        for id = 2:length(save_figs)
            saveas(save_figs(id).fig,save_figs(id).savepath);
        end
    end
end

%% Constellation code from AE 725
if 0
    %% Cubesat Orbit Validation Code
    numsats = 50;
    sat_alongsepkm = 10/1e3; %km Total along-track separation
    sat_crosssepkm = 1/1e3; %km Total cross-track separation
    sat_sepkm = sqrt(sat_alongsepkm^2+sat_crosssepkm^2); %km Total separation
    num_periods = 10; %Periods for orbit integration
    RE = 6378.145; %km Radius of the Earth

    %% Orbital Settings
    r_obs = 400; %km Radius above the Earth during observations (=periapsis)
    e_orbit = 0.01; %~ Eccentricity
    i_orbit = deg2rad(86); %radians Inclination
    nu_orbit = 0; %radians True Anomaly
    RA_orbit = 0; %radians Right Ascension of the Ascending Node
    w_orbit = 0; %radians Argument of Perigee
    mu_orbit = 3.986004e5; %km^3/s^2

    %Do the calculations and integrations
        %Changes in true anomaly likely due to waiting to push the cubesats out
            %of the rocket
        kickspeed = 2; %m/s  %Based on CubeSat PPOD
        waittime = sat_alongsepkm*1e3/kickspeed; %seconds Time to wait between launching sats
        nuang = atan((kickspeed*waittime)*1e-3/(r_obs+RE));
        %Changes in inclination for cross-track separation
        iang = atan(sat_crosssepkm/(r_obs+RE)); %Radians Change in inclination
    %Apply changes and iterate for each sat
    for j = 1:numsats
        ichng = iang*(j-1);
        nuchng = -nuang*(j-1);
        [R1data{j}, R2data{j}, R3data{j}, tvec, V1data{j}, V2data{j}, V3data{j}]...
            = OrbitRun(r_obs, e_orbit, i_orbit+ichng, nu_orbit+nuchng, ...
            RA_orbit, w_orbit, mu_orbit, num_periods);
    end
end

%% Ground Track Plot
%Mercator

%Polar Stereographic

%% Post processing