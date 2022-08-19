%% Settings
%Ground Station details
gs_windows = struct('Lat',38.971669,'Lon',-95.23525,'Elev',0,'BW',131,'Bands',{{'UHF','S','X'}},'Name','HawksNest'); %Lawrence, KS
gs_windows(end+1) = struct('Lat',64+48/60,'Lon',-147-39/60,'Elev',0,'BW',131,'Bands',{{'S','X'}},'Name','NorthPole'); %North Pole

%Sounding Target details
radar_Ant = struct('Lat',-87,'Lon',90,'Rad',2500); %Antartica
radar_Gre = struct('Lat',71,'Lon',-43,'Rad',1300); %Greenland

%% Automated Calculations
%Find geocircles for radar targets
[Antcirc_lats, Antcirc_lons] = geocircle(radar_Ant.Lat,radar_Ant.Lon,radar_Ant.Rad);
[Grecirc_lats, Grecirc_lons] = geocircle(radar_Gre.Lat,radar_Gre.Lon,radar_Gre.Rad);
