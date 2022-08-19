clc
% clearvars
tic
addpath(genpath('Subfunctions'))
%% Settings
%Simulation Settings
N =20; %Number of orbits
timesteps = 206; %Time steps per orbit

Settings_General; %Loads paths
Settings_Radar; %Gives us wavelength and valid array spacing
Settings_Constellation; %Gives us const_struct
Settings_Targets; %Gives geocircles for targets
% Settings_STK

%% Propagate Orbit
if 1
    %Matlab Simulation
    [sats] = Sim_Constellation_matlab_J2Drag(sats,init_utcvec,N,timesteps);
elseif 0
    %STK Simulation
    [sats] = Sim_Constellation_STK_J2Drag(sats,init_utcvec,N,timesteps);
end

%% Do Post Calculations
driving_sat_id = ceil(length(const_vec)/2);%find(strcmp('TX',const_vec)==1,1,'first');
num_states = length(sats(driving_sat_id).states);

%Do calculations
[out_struct,fail_mat,mindist,mindistt] = constellation_calcs(sats,driving_sat_id,avg_array_rangewl*wavelength,fail_min_off,fail_prob);
%Parse outputs
state_pos_cells = out_struct.state_pos;
alt_mat = out_struct.alt;
along_pos_mat = out_struct.along_pos;
cross_pos_mat = out_struct.cross_pos;
elev_pos_mat = out_struct.elev_pos;
inplanemag_pos_mat = out_struct.inplane_mag;
mag_pos_mat = out_struct.total_mag;
avg_mag_m = out_struct.avg_total_mag;
avg_inplanemag_m = out_struct.avg_inplane_mag;
Grecheck = out_struct.Greencheck;
Antcheck = out_struct.Antcheck;
datacheck = out_struct.datacheck;
datametric_totallinekm = out_struct.data_total; %line-km of collected data total
datametric_Grelinekm = out_struct.data_Green; %line-km of data from Greenland
datametric_Antlinekm = out_struct.data_Ant; %line-km of data from Antarctica

%Get close approach points
closeapproach_dist_m = 2;
closecheck = zeros(size(datacheck));
for statid = 1:length(avg_mag_m)
    if avg_mag_m(statid)<closeapproach_dist_m
        closecheck(statid) = 1;
    end
end

%% Print results
if any(any(fail_mat))
    fprintf('\nThere was a failure for this constellation configuration.\n')
end
[~,tparse,tunitstr] = parsetoc(mindistt);
fprintf('The minimum distance of %.2f meters occured at %.f %s\n',mindist,tparse,tunitstr)
%% Plot results
if 1
    %Find average offsets over greenland
    avg_off_Gre = zeros(3,num_sats);
    avg_Gre_check = find(and(Grecheck,datacheck)==1,1,'first');
    avg_off_Ant = zeros(3,num_sats);
    avg_Ant_check = and(Antcheck,datacheck);
    for satid = 1:num_sats
        avg_off_Gre(:,satid)= [mean(along_pos_mat(avg_Gre_check,satid)) ...
            mean(cross_pos_mat(avg_Gre_check,satid)), mean(elev_pos_mat(avg_Gre_check,satid))].*1000;
        avg_off_Ant(:,satid)= [mean(along_pos_mat(avg_Ant_check,satid)),...
            mean(cross_pos_mat(avg_Ant_check,satid)), mean(elev_pos_mat(avg_Ant_check,satid))].*1000;
    end
    
    %Set config plot limits
    crosslims = [-15 15];
    alonglims = [-200 200];
    elevlims = [-1 1];
    %Plot the average offsets for Greenland
    figure(22); clf
        %Overhead
    subplot(2,1,1)
    leg = {};
    for satid = 1:num_sats
        apnow = plot(avg_off_Gre(2,satid),avg_off_Gre(1,satid),'o');
        set(apnow,'MarkerFaceColor',apnow.Color)
        hold on
        leg{end+1} = sats(satid).name;
    end
    hold off
    grid on
    ylabel('Along Track Offset (m)')
    xlabel('Cross Track Offset (m)')
    title(sprintf('Along-track and Cross-track Offset over Greenland\nFor %.0f Orbits',N))
    xlim(crosslims)
    ylim(alonglims)
    legend(leg,'Location','best')
    
        %In plane
    subplot(2,1,2)
    leg = {};
    for satid = 1:num_sats
        apnow = plot(avg_off_Gre(2,satid),avg_off_Gre(3,satid),'o');
        set(apnow,'MarkerFaceColor',apnow.Color)
        hold on
        leg{end+1} = sats(satid).name;
    end
    hold off
    grid on
    ylabel('Elevation Offset (m)')
    xlabel('Cross Track Offset (m)')
    title(sprintf('Elevation and Cross-track Offset over Greenland\nFor %.0f Orbits',N))
    xlim(crosslims)
    ylim(elevlims)
    legend(leg,'Location','best')
    
    %Plot the average offsets for Antarctica
    figure(23); clf
        %Overhead
    subplot(2,1,1)
    leg = {};
    for satid = 1:num_sats
        apnow = plot(avg_off_Ant(2,satid),avg_off_Ant(1,satid),'o');
        set(apnow,'MarkerFaceColor',apnow.Color)
        hold on
        leg{end+1} = sats(satid).name;
    end
    hold off
    grid on
    ylabel('Along Track Offset (m)')
    xlabel('Cross Track Offset (m)')
    title(sprintf('Along-track and Cross-track Offset over Antarctica\nFor %.0f Orbits',N))
    xlim(crosslims)
    ylim(alonglims)
    legend(leg,'Location','best')
        %In plane
    subplot(2,1,2)
    leg = {};
    for satid = 1:num_sats
        apnow = plot(avg_off_Ant(2,satid),avg_off_Ant(3,satid),'o');
        set(apnow,'MarkerFaceColor',apnow.Color)
        hold on
        leg{end+1} = sats(satid).name;
    end
    hold off
    grid on
    ylabel('Elevation Offset (m)')
    xlabel('Cross Track Offset (m)')
    title(sprintf('Elevation and Cross-track Offset over Antarctica\nFor %.0f Orbits',N))
    xlim(crosslims)
    ylim(elevlims)
    legend(leg,'Location','best')
    
    %Plot the closest approach
    [close_avg_m, closeid] = min(avg_mag_m);
    figure(24); clf
        %Overhead
    subplot(2,1,1)
    leg = {};
    for satid = 1:num_sats
        apnow = plot(cross_pos_mat(closeid,satid)*1000,along_pos_mat(closeid,satid)*1000,'o');
        set(apnow,'MarkerFaceColor',apnow.Color)
        hold on
        leg{end+1} = sats(satid).name;
    end
    hold off
    grid on
    ylabel('Along Track Offset (m)')
    xlabel('Cross Track Offset (m)')
    title(sprintf('Along-track and Cross-track Offset\nFor Closest Approach'))
    xlim(crosslims)
    ylim(alonglims)
    legend(leg,'Location','best')
        %In plane
    subplot(2,1,2)
    leg = {};
    for satid = 1:num_sats
        apnow = plot(cross_pos_mat(closeid,satid)*1000,elev_pos_mat(closeid,satid)*1000,'o');
        set(apnow,'MarkerFaceColor',apnow.Color)
        hold on
        leg{end+1} = sats(satid).name;
    end
    hold off
    grid on
    ylabel('Elevation Offset (m)')
    xlabel('Cross Track Offset (m)')
    title(sprintf('Elevation and Cross-track Offset\nFor Closest Approach'))
    xlim(crosslims)
    ylim(elevlims)
    legend(leg,'Location','best')
    
end
    


%% Plot the outputs
%Initialize relative motion figures
if 1
    fig_ground = figure(12); clf;
    %Start Ground track plot
    statecell = struct2cell(sats(driving_sat_id).states); 
    statefields = fieldnames(sats(driving_sat_id).states);
    latname = 'lat';
    latid = find(strcmp(statefields,latname));
    lonid = latid+1;
    lats = reshape(cell2mat(statecell(latid,:,:)),[1 num_states]);
    lons = reshape(cell2mat(statecell(lonid,:,:)),[1 num_states]);

    set(0,'currentfigure',fig_ground);
    geoaxes();
    %Plot all points
%     geoscatter(lats,lons,'.b')
    hold on
    %Plot Antarctica Circle
    geoscatter(Antcirc_lats,Antcirc_lons,'.m')
    %Plot Greenland Circle
    geoscatter(Grecirc_lats,Grecirc_lons,'.m')
    %Plot Close Approach Points 
%     geoscatter(lats(logical(closecheck)),lons(logical(closecheck)),'.r');
    %Plot Data Collection Points over Greenland
    geoscatter(lats(and(Grecheck,datacheck)),lons(and(Grecheck,datacheck)),'.g');
    %Plot Data Collection Points over Antarctica
    geoscatter(lats(and(Antcheck,datacheck)),lons(and(Antcheck,datacheck)),'.c');
    
    hold off
    geolimits([-90 90], [-180 180])
    title('Ground Track')
    FontWidthandPos
end
if 0
    
    fig_offset = figure(13); clf;
    %Start distance magnitude plot
    timename = 't';
    timeid = find(strcmp(statefields,timename));
    timevec = reshape(cell2mat(statecell(timeid,:,:)),[1 num_states]);
    posname = 'R';
    drivingsatcell = struct2cell(sats(driving_sat_id).states);
    posid = find(strcmp(statefields,posname));
    pos_driving_sat=reshape(cell2mat(drivingsatcell(posid,:,:)),[3 num_states]);
    set(0,'currentfigure',fig_offset);
    %Do latitude plot
    subplot(2,1,1)
    plot(timevec,lats)
    grid on
    ylabel('Latitude (deg)')
    ylim([-90 90])
    for sat_id = 1:num_sats
        statecellnow = struct2cell(sats(sat_id).states);
        %Do magnitude of offsets
        subplot(2,1,2)
        pos_now_sat=reshape(cell2mat(statecellnow(posid,:,:)),[3 num_states]);
        pos_now_off = pos_now_sat-pos_driving_sat;
        mag_pos_off = vecnorm(pos_now_off);
        hold on
        plot(timevec,mag_pos_off*1e3)
        hold on
    end
    hold off
    grid on
    xlabel('Time (sec)')
    ylabel('Offset Magnitude (m)')
    ylim([0 50])
end

if 0
    fig_Gre = figure(10); clf;
    fig_Ant = figure(11); clf;
    %iterate
    plotnow = true; plotmaxticks = 5; plottick = 0;
    green_linekm = 0; arctic_linekm = 0;
    for sid = 1:num_states
        %tick
        plottick = plottick+1;
        %Get driving sat lat lon
        latnow = sats(driving_sat_id).states(sid).lat;
        lonnow = sats(driving_sat_id).states(sid).lon;
        %Determine if over Greenland
        Greencheck = abs(geodist(radar_Gre.Lat,radar_Gre.Lon,latnow,lonnow))<radar_Gre.Rad;

        %Determine if over Antarctica
        if ~Greencheck
            Antcheck = abs(geodist(radar_Ant.Lat,radar_Ant.Lon,latnow,lonnow))<radar_Ant.Rad;
        else
            Antcheck = false;
        end

        %Check if plotting
        if Greencheck || Antcheck || (plottick-plotmaxticks)>=0
            plottick = 0;

            %Plot on ground track
            set(0,'currentfigure',fig_ground);
            delete(findall(gcf,'Tag','sat'))
            hold on
            for sat_id = 1:num_sats
                %Plot new points
                grndpoint = geoscatter(sats(driving_sat_id).states(sid).lat,sats(driving_sat_id).states(sid).lon,'.g');
                set(grndpoint,'Tag','sat','SizeData',1000)
            end
            drawnow
            %Plot on offsets
            if Antcheck
                tickcol = 'b';
            elseif Greencheck
                tickcol = 'g';
            else
                tickcol = 'k';
            end
            set(0,'currentfigure',fig_offset);
            delete(findall(gcf,'Tag','tick'))
            %Plot on the latitude
            subplot(2,1,1)
            hold on
            plot([timevec(sid) timevec(sid)],[-200 200],tickcol,'tag','tick')
            hold off
            %Plot on the offset magnitude
            subplot(2,1,2)
            hold on
            plot([timevec(sid) timevec(sid)],[-1000 1000],tickcol,'tag','tick')
        %     ylim([-10 10])
            hold off
            drawnow
        end
        %REMOVE 
        Greencheck = true;
        %Plot checks
        if Greencheck || Antcheck
            %Plot relative positions based on check
            if Greencheck
                set(0,'currentfigure',fig_Gre); %Set to current figure
            elseif Antcheck
                set(0,'currentfigure',fig_Ant); %Set to current figure
            end
            subplot(1,4,1)
            plot(cross_pos_mat(sid,:)*1e3,elev_pos_mat(sid,:)*1e3,'--o')
            title('Cross-track & Elevation')
            xlabel('Cross-track (m)')
            ylabel('Elevation (m)')
            ylim([-10 10])
            xlim([-10 10])
            grid on
            subplot(1,4,2)
            plot(inplanemag_pos_mat(sid,:)*1e3/wavelength,'--o')
            title('In-Plane Magnitude')        
            xlabel('Satellite #')
            ylabel('Offset Magnitude (wavelength)')
            ylim([0 10])
            grid on
            subplot(1,4,3)
            plot(cross_pos_mat(sid,:)*1e3,along_pos_mat(sid,:)*1e3,'--o')
            title('Cross-track & Along-track')        
            xlabel('Cross-track (m)')
            ylabel('Along-track (m)')
            ylim([-250 250])
            xlim([-10 10])
            grid on
            subplot(1,4,4)
            plot(mag_pos_mat(sid,:)*1e3,'--o')
            title('Total Magnitude')        
            xlabel('Satellite #')
            ylabel('Offset Magnitude (m)')
            ylim([0 250])
            grid on
            drawnow            
        end
    end
end
toc