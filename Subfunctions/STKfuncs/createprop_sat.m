%% Create Satellite Object
if ~scenario.Children.Contains('eSatellite',satname)
    %Add satellite object
    scenario.Children.New('eSatellite',satname);
    %Grab satname
    satnames{satid} = satname;
else
    %Message
    root.ExecuteCommand(sprintf('Message / 3 "Satellite %s Already Exists"',satname));
    %Grab satellite object
    scenario.Children.Item(satname);
    satnames{satid} = satname;            
end

%% Define Model file
if contains(satname,'RX')
    modelfnnow = RXmodelfn;
elseif contains(satname,'TX')
    modelfnnow = TXmodelfn;
end
modelcmdstr = sprintf('VO */Satellite/%s Model File "%s"',satname,fullfile(modelbasedir,modelfnnow));
root.ExecuteCommand(modelcmdstr);


%% Propagate Mission STK
switch statetype
    case 'Classical'
        %Grab orbital elements
        a_str = sprintf('%.6f',COE0.a_SemiMajorAxis*1000); %m
        e_str = sprintf('%.6f',COE0.e_Eccentricity); %~
        i_str = sprintf('%.6f',rad2deg(COE0.i_Inclination)); %deg
        w_str = sprintf('%.6f',rad2deg(COE0.w_ArgumentofPerigee)); %deg
        RA_str = sprintf('%.6f',rad2deg(COE0.RA_RightAscension)); %deg
        M_str = sprintf('%.6f',rad2deg(COE0.M_MeanAnomaly)); %deg

        %Set satellite orbital elements and propagator type
        setstatestr = sprintf('SetState */Satellite/%s %s %s "%s" "%s" %.0f %s "%s" %s %s %s %s %s %s',...
            satname,statetype,propagatortype,timestartstr,timeendstr,stepsize,coordsyststr,timestartstr,...
            a_str,e_str,i_str,w_str,RA_str,M_str);

    case 'Cartesian'
        %Grab state details
        x_str = sprintf('%.6f',Rnow(1)*1000); %m
        y_str = sprintf('%.6f',Rnow(2)*1000); %m
        z_str = sprintf('%.6f',Rnow(3)*1000); %m
        xv_str = sprintf('%.6f',Vnow(1)*1000); %m/s
        yv_str = sprintf('%.6f',Vnow(2)*1000); %m/s
        zv_str = sprintf('%.6f',Vnow(3)*1000); %m/s

        %Set satellite orbital elements and propagator type
        setstatestr = sprintf('SetState */Satellite/%s %s %s "%s" "%s" %.3f %s "%s" %s %s %s %s %s %s',...
            satname,statetype,propagatortype,timestartstr,timeendstr,stepsize,coordsyststr,timestartstr,...
            x_str,y_str,z_str,xv_str,yv_str,zv_str);
end

%Send propagation command
root.ExecuteCommand(setstatestr);
%Send message about command
root.ExecuteCommand(['Message / 1 "Satellite ' satname ' Loaded"']);

%% Turn on Drag
switch propagatortype
    case 'HPOP'
        areamassratio = SCnow.A*1000*1000/SCnow.m; %m^2/kg
        dragcmdstr = sprintf('HPOP */Satellite/%s Drag On %.2f %.4f "NRLMSISE 2000" Manual %.0f %.0f %.01f',...
            satname,SCnow.CD,areamassratio,F10abs,F10avg,max(Apvec));
        root.ExecuteCommand(dragcmdstr);
end