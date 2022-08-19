function [Xdot] = OrbitDerivFunc_J2Drag(X,RE,muE,J2,CD,A,m,rho0,r0,H,thetadot,yr,dayofyear,dt,F10avg,F10day,Ap)
%Derivative function to be passed to the numerical solver
%X should have the form [r1;r2;r3;v1;v2;v3] in km and km/s 

%Constants for Earth as the central body
if nargin<2
    RE = 6378.145; %km
    muE =  3.986004e5; %km^3/s^2
    J2 = 0.00108248;
    CD = 2; %Drag Coefficient
    A = .36/1000/1000;%km^2
    m = 10; %kg
    rho0 = 4e-4; %kg/km^3 %check
    r0 = 7298145/1000; %km
    H = 200; %km
    thetadot = 7.29211585530066e-5;%rad/s
end

%Recalculate time
datenow = datetime([yr 1 1 1 0 0])+days(dayofyear)+seconds(dt);
dayofyear = day(datenow,'dayofyear');
yr = year(datenow);
dt = second(datenow,'secondofday');

%Position terms
x = X(1); y = X(2); z = X(3); dx = X(4); dy = X(5); dz = X(6);
r = sqrt(x^2+y^2+z^2);

%Velocity terms
Xdot(1) = dx;
Xdot(2) = dy;
Xdot(3) = dz;

%Acceleration terms
%2 Body
ax = -muE*x/r^3;
ay = -muE*y/r^3;
az = -muE*z/r^3;

%With J2
ax = ax + 15*muE*RE^2*J2*z^2*x/(2*r^7)-3*muE*RE^2*J2*x/(2*r^5);
ay = ay + 15*muE*RE^2*J2*z^2*y/(2*r^7)-3*muE*RE^2*J2*y/(2*r^5);
az = az + 15*muE*RE^2*J2*z^2*z/(2*r^7)-9*muE*RE^2*J2*z/(2*r^5);

%With Drag and J2
%Calculate VA
VAx = dx+thetadot*y;
VAy = dy-thetadot*x;
VAz = dz;
VA = sqrt((VAx)^2+(VAy)^2+VAz^2);
if 0
    %Exponential density
    %Calculate rhoA
    rhoA = rho0*exp(-(r-r0)/H);
else
    %NRLMSISE density
    %Calculate Lat, Lon, and altitude above ellipse
    [lat,lon,rellps] = ECEF2latlon(X(1:3),dt);
    rellps = rellps*1000; %m from km
    %Calculate rhoA
    [~, rhovec] = atmosnrlmsise00(rellps,lat,lon,yr,dayofyear,dt,F10avg,F10day,Ap);
    rhoA = rhovec(6); %kg/m^3
    rhoA = rhoA*1e9; %kg/km^3
end

%Solve for drag components
adrag = -1/2*CD*(A/m)*rhoA*VA;

ax = ax + adrag*VAx;
ay = ay + adrag*VAy;
az = az + adrag*VAz;


Xdot(4:6) = [ax;ay;az];
Xdot = Xdot';
end