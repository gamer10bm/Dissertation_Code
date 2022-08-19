function [total_torque_inst] = momentumbuildup(SC,R,V,eclipse_chk)
% [total_torque_inst] = momentumbuildup(SC,R,V,eclipse_chk)
%
% Inputs: SC - Spacecraft structure needed for constants
%             R - Position vector in ECEF (km)
%             V - Velocity vector in ECEF (km/s)
%             eclipse_chk - Boolean value to consider solar torque
%                                  (true = in eclipse, false = in sunlight)
%
% Outputs: total_torque_inst - Instantneous torque (N-m)
%
% Created by Bailey Miller 2/12/2021
% See also: initiate_SC_model.m, ECEF2latlon.m,

%% Settings
%Torques to consider
magt_bool = true;
dragt_bool = true;
gravt_bool = true;
vertorflat = 'vert';

if nargin<4
    solt_bool = true;
else
    solt_bool = ~eclipse_chk;
end

%Earth settings
B0 = 3.1e-5; %tesla Magnetic field at Earth surface
Is = 1376; %W/m^2 incident solar radiation @ Earth
c = physconst('Lightspeed'); %m/s
Re = SC.Re*1e3;%m Earth Radius
mue = SC.mu*1e3^3;%m^3/s^2 gravitational parameter at Earth
p0 = 101325; %Pa sea level pressure
T0 = 288.15; %K sea level temperature
g0 = 9.80665; %m/s^2 sea level acceleration due to gravity
Rconst = 8.31447; %J/(mol*K) universal gas constant
Mdryair = 0.0289654; %kg/mol molar mass of dry air
Ltemp = 0.0065; %K/m temperature lapse rate

%Satellite Settings
long_dim = .3; %m
short_dim = .1; %m
mass = SC.m;%kg

%Magnetic settings
Mmters = .1; %A-m^2 estimate for magnetorquers residual dipole

%Solar settings
alpha = 0; %deg Solar incidence angle (conservative to be 0)
q_craft = 0.3; %reflectance of spacecraft (between 0 and 1)
L_solarm = 0.15; %Since the full surface is perpendicularly torqued

%Drag settings
Cd = SC.CD; 
f107a = 109.8;
f107d = 102.6;
Ap = 50;

%Gravity Gradient settings
theta = pi; %as this is modulated below

%% Calculations
%Do quick calculations based on orientation
if strcmp('vert',vertorflat)
    z_dim = long_dim;
    y_dim = short_dim;
    x_dim = short_dim;
    thetmters = [0 90 90]; %deg %Angle between the magnetic field lines and perpendicular to the coil
else
    z_dim = short_dim;
    y_dim = long_dim;    
    x_dim = short_dim;
    thetmters = [90 0 0]; %deg 
end
Area_cross = x_dim*y_dim; %m^2 
L_drag = y_dim/2; %m estimate


%Use state to find orbital parameters
% [~,h,i,raan,e,arg_per,nu,a,T] = RV2coe(R,V,mue);
Rmag = norm(R)*1e3; %m
Vmag = norm(V)*1e3; %m
alt =Rmag-Re; %m height above the surface
[lat, lon] = ECEF2latlon(R); %degrees latitude


% Momentum buildup code based on Brown Design Book (560) ch. 5
    %Total torque = Drag + Solar + Gravity + magnetic
if magt_bool 
    %% Magnetic Torque
    %Find the magnetic field strength
    B = B0*Re^3/Rmag^3*sqrt(3*sind(lat)^2+1); %tesla
    %Find the magnetic torque
    magT = 0;
    for m = 1:length(thetmters)
        %Do for each magnetorquer
        magT = magT + Mmters*B*sind(thetmters(m));
    end
else
    magT = 0;
end

if solt_bool
    %% Solar Torque
    %Calculate solar pressure
    Ps = Is/c; %N/m^2 solar pressure
    %Calculate the projected area A*cos(alpha)
    Aproj = Area_cross*cosd(alpha);
    %Calculate the solar torque
    solT = Ps*Aproj*L_solarm*(1+q_craft);
else
    solT = 0;
end

if dragt_bool
    %% Drag Torque
    %Determine atmospheric density based on exponential model
    [T, Rho] = atmosnrlmsise00(Rmag-Re,lat,lon,2022,01,1,f107a,f107d,Ap);
    density = Rho(6);

%     p = p0*(1-L_drag*alt/T0)^(g0*Mdryair/(Rconst*Ltemp));
%     T = T0-L_drag*alt;
%     rho = p*Mdryair/(Rconst*T); %kg/m^3
    %Calculate drag force
    Drag = .5*density*Vmag^2*Cd*Area_cross; %N
    %Calculate the torque
    dragT = Drag*L_drag;
else
    dragT = 0;
end

if gravt_bool
    %% Gravity Gradient Torque
    %Find the moment of inertias http://dynref.engr.illinois.edu/rem.html 
    Ix = mass/12*(y_dim^2+z_dim^2); %.0217;
    Iy = mass/12*(z_dim^2+x_dim^2); %.0010;
    Iz = mass/12*(y_dim^2+x_dim^2); %.0217;
    if Ix<Iy
        Idiff = abs(Iz-Ix);
    else
        Idiff = abs(Iz-Iy);
    end
    %Calculate the torque
    gravT = 3*mue*Idiff*theta/Rmag^3;
else
    gravT = 0;
end

%Total torque
total_torque_inst = gravT+dragT+magT+solT;
end