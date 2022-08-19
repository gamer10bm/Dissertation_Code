%% Export propagation data
%Initialize outputs
clearvars sat_states
sat_states(num_states) = struct('t',0,'OE',COE0,'R',[0;0;0],'V',[0;0;0],'lat',0,'lon',0,'alt',0);
if 1
    %% Grab satellite position data for each propagation time step (Try 2: Create Report)
    timeendsimstr = datestr(dt_init+seconds(time(end)),'dd mmmm yyyy HH:MM:SS.FFF');
    reportoptsstr = sprintf('TimeStep %.3f TimePeriod "%s" "%s"',stepsize,timestartstr,timeendsimstr);
    %Generate and send call for position
    poscmdstr = sprintf('Report_RM */Satellite/%s Style "J2000 Position Velocity" %s',...
        satname,reportoptsstr);
    satposresult = root.ExecuteCommand(poscmdstr);
    %Generate and send call for position
    llacmdstr = sprintf('Report_RM */Satellite/%s Style "LLA Position" %s',...
        satname,reportoptsstr);
    satllaresult = root.ExecuteCommand(llacmdstr);
    %Parse results for each time step
    for stateid = 1:num_states
        sat_states(stateid).t = time(stateid);
        %Parse position data
        splitposresult = strsplit(satposresult.Item(stateid),',');
        %1: Time String 2:X 3:Y 4:Z 5:vX 6:vY 7:vZ
        sat_states(stateid).timestr = splitposresult{1};
        posvec = [str2double(splitposresult{2}); str2double(splitposresult{3}); str2double(splitposresult{4})];
        sat_states(stateid).R = posvec; %km
        velvec = [str2double(splitposresult{5}); str2double(splitposresult{6}); str2double(splitposresult{7})];
        sat_states(stateid).V = velvec; %km/s
        %Parse lla data
        splitllaresult = strsplit(satllaresult.Item(stateid),',');
        %1: Time String 2:Lat 3:Lon 4:Alt 5:Latrate 6:Lonrate 7:Altrate
        sat_states(stateid).lat = str2double(splitllaresult{2});
        sat_states(stateid).lon = str2double(splitllaresult{3});
        sat_states(stateid).alt = str2double(splitllaresult{4})/1000; %km
        %Get current orbital elements
        sat_states(stateid).OE = RV2coe(posvec,velvec,SCnow.mu);
    end
else
    %% Grab satellite position data for each propagation time step (Try 1: Position Command)
    for stateid = 1:num_states
        %Get current time
        sat_states(stateid).t = time(stateid);
        secfromstart = time(stateid); %sec
        %Generate new time string
        dt_now = dt_init+seconds(secfromstart);
        timestr = datestr(dt_now,'dd mmmm yyyy HH:MM:SS.FFF');
        %Do Position Calls for this satellite
        resultnow = root.ExecuteCommand(sprintf('Position */Satellite/%s "%s"',satname,timestr));
        %Parse result
        resultstr = resultnow.Item(0);
        sepstrs = strsplit(resultstr,' ');
        %1-6: Are LLA %7-12 Fixed %13-18 Inertial Position
        %Load result to data
        sat_states(stateid).lat = str2double(sepstrs{1}); %deg
        sat_states(stateid).lon = str2double(sepstrs{2}); %deg
        sat_states(stateid).alt = str2double(sepstrs{3})/1000; %km
        posvec = [str2double(sepstrs{13}); str2double(sepstrs{14}); str2double(sepstrs{15})];
        sat_states(stateid).R = posvec/1000; %km
        velvec = [str2double(sepstrs{16}); str2double(sepstrs{17}); str2double(sepstrs{18})];
        sat_states(stateid).V = velvec/1000; %km
        %Get current orbital elements
        sat_states(stateid).OE = RV2coe(posvec,velvec,SCnow.mu);
    end
end