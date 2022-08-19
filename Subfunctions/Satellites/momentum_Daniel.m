function [total_torque_inst] = momentum_Daniel(SC,R,V,epoch,t)
% [total_torque_inst] = momentum_Daniel(SC,R,V,eclipse_chk)
%
% Inputs: SC - Spacecraft structure needed for constants
%             R - Position vector in ECEF (km)
%             V - Velocity vector in ECEF (km/s)
%
% Outputs: total_torque_inst - Instantneous torque (N-m)
%
% Created by Daniel and adapted by Bailey Miller 7/1/2021
% See also: initiate_SC_model.m, ECEF2latlon.m,

%% Settings
%Torques to consider
magt_bool = false; %Not created by Daniel (not in original code)
dragt_bool = true;
gravt_bool = true;
solt_bool = true;
vertorflat = 'vert';

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
% [lat, lon] = ECEF2latlon(R); %degrees latitude
%% Daniel Code (Adapted)
xArea = 0.1*0.3405;
yArea = 0.1*0.1;
zArea = 0.1*0.3405;
% Attitude
C = eye(3);
% Rotate R vector into S/C body frame
Rb = C*R;
Vb = C*V;
Rsun = C*sun_vector(epoch, t);
% Not sure if I trust these values...
J = diag([1082793.489*20; 4354.008742*223; 1085536.436*20])/1000^3;
%J = diag([1082793.489; 4354.008742; 1085536.436])/1000^3;
CoM = [1.14122; 19.92782; 0.54375]*0.001; %CoM offset

% Momentum buildup code developed by Daniel Owens
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
    p0 = 4.58e-6;
    
    eclipse = 1;
    eclipseAng = atan(Re/Rmag);
    if dot(Rb, Rsun) > 0
        ang = acos(dot(Rb, Rsun)/Rmag);
        if ang < eclipseAng
            eclipse = 0;
        end
    else
        ang = 3;
    end

    Fsrp = p0*[Rsun(1)*xArea; Rsun(2)*yArea; Rsun(3)*zArea];
    solT = eclipse*Fsrp.*CoM;
else
    solT = 0;
end

if dragt_bool
    %% Drag Torque
    
    % Drag
    rho = HighAltDrag(Rmag-Re);
    q = 0.5*rho*Vmag^2;
%     cd = 1.28;               % Flat plate
    Fd = q*Cd*[Vb(1)/Vmag*xArea; Vb(2)/Vmag*yArea; Vb(3)/Vmag*zArea];
    
    dragT = Fd.*CoM;
else
    dragT = 0;
end

if gravt_bool
    %% Gravity Gradient Torque
    gravT = 3*mue/Rmag^5*cross(Rb, J*Rb);
else
    gravT = 0;
end
%Total torque
total_torque_inst = norm(gravT+dragT+magT+solT);
end

function rho = HighAltDrag(h)
    hKm = h/1000;
    
    if(hKm > 600)
        hKm = 600;
    end
    
    h0 = -1;
    rho0 = -1;
    H = -1;
    
    if (hKm < 25) 
        h0 = 0;
        rho0 = 1.225;
        H = 7.249;
    elseif (hKm < 30)
        h0 = 25;
        rho0 = 3.899e-2;
        H = 6.349;
    elseif (hKm < 40)
        h0 = 30;
        rho0 = 1.774e-2;
        H = 6.682;
    elseif (hKm < 50)
        h0 = 40;
        rho0 = 3.972e-3;
        H = 7.554;
    elseif (hKm < 60)
        h0 = 50;
        rho0 = 1.057e-3;
        H = 8.382;
    elseif (hKm < 70)
        h0 = 60;
        rho0 = 3.206e-4;
        H = 7.714;
    elseif (hKm < 80)
        h0 = 70;
        rho0 = 8.770e-5;
        H = 6.549;
    elseif (hKm < 90)
        h0 = 80;
        rho0 = 1.905e-5;
        H = 5.799;
    elseif (hKm < 100)
        h0 = 90;
        rho0 = 3.396e-6;
        H = 5.382;
    elseif (hKm < 110)
        h0 = 100;
        rho0 = 5.297e-7;
        H = 5.877;
    elseif (hKm < 120)
        h0 = 100;
        rho0 = 9.661e-7;
        H = 7.263;
    elseif (hKm < 130)
        h0 = 120;
        rho0 = 2.438e-8;
        H = 9.473;
    elseif (hKm < 140)
        h0 = 130;
        rho0 = 8.484e-9;
        H = 12.636;
    elseif (hKm < 150)
        h0 = 140;
        rho0 = 3.845e-9;
        H = 16.149;
    elseif (hKm < 180)
        h0 = 150;
        rho0 = 2.070e-9;
        H = 22.523;
    elseif (hKm < 200)
        h0 = 180;
        rho0 = 5.464e-10;
        H = 29.740;
    elseif (hKm < 250)
        h0 = 200;
        rho0 = 2.789e-10;
        H = 37.105;
    elseif (hKm < 300)
        h0 = 250;
        rho0 = 7.248e-11;
        H = 45.546;
    elseif (hKm < 350)
        h0 = 300;
        rho0 = 2.418e-11;
        H = 53.628;
    elseif (hKm < 400)
        h0 = 350;
        rho0 = 9.518e-12;
        H = 53.298;
    elseif (hKm < 450)
        h0 = 400;
        rho0 = 3.725e-12;
        H = 58.515;
    elseif (hKm < 500)
        h0 = 450;
        rho0 = 1.585e-12;
        H = 60.828;
    elseif (hKm <= 600)
        h0 = 500;
        rho0 = 6.967e-13;
        H = 63.822;
    end

    rho = rho0*exp(-(hKm - h0)/H);
    
end