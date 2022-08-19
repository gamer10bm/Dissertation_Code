function [lat, lon, h_ellps] = ECEF2latlon(R,dt)
% [lat, lon] = ECEF2latlon(R,dt)
%
% This function converts a position vector (R) from earth centered earth
% fixed frame to latitude and longitude
%
% Inputs:
%       R:      Position Vector in Earth Centered Earth Fixed (km)
%       dt:     Time since simulation start (sec)
%
% Outputs:
%       lat:    Latitude (deg)
%       lon:    Longitude (deg)
%       H_ellps: Height above the ellipse (km)
%
% Adapted from: 
% Sources: 
%   Fundamentals of Astrodynamics and Applications by Vallado (Sec 3.4.4)
%       -Base algorithm for transformation
%   https://smallsats.org/2013/01/20/satellite-ground-track-iss/
%       -For earth rotation
% Created by: Bailey Miller 10/4/2021
%

%% Check inputs
if nargin < 2
    dt = 0; %Also considered as no earth rotation
end
Rnow = [0;0;0];
if size(R,2) >1 && size(R,1) == 1 %transpose vector
    Rnow = transpose(R);
elseif size(R,2) >1 && size(R,1) == 3 %Not supported
    error('Vectors not yet supported')
else
    Rnow = R;
end
if 1
    %% Rotate R to geocentric equatorial frame (smallsat source)
    %Consider Earth's rotation
    we = 360*(1 + 1/365.25)/(3600*24);      % Earth's rotation [deg/s]
    fi_earth = we*dt; %degrees of earth spin since sim start
    %Calculate rotation matrix
    Rotmat = [cosd(fi_earth), sind(fi_earth),0;...
            -sind(fi_earth),cosd(fi_earth),0;0,0,1];
    %Calculate R in geocentric equatorial frame
    Rgeo = Rotmat*Rnow;


    %% Convert Rgeo to lat and lon (Vallado source)
    %Parse R vector
    rx = Rgeo(1); ry = Rgeo(2); rz = Rgeo(3);
    %Calculate rdsat
    rdsat = sqrt(rx^2+ry^2);
    %Find alpha
    alpha = acos(rx/rdsat); %radians
    if ry < 0
        alpha = -alpha;
    end
    %Set longitude to alpha
    lon = rad2deg(alpha);
    %Find d
    d = atan2(rz,rdsat);
    %Iterate to find latitude
    epsilon = 1e-12;
    phi_gd = d;
    phi_gd_old = phi_gd+epsilon*8;
    Re = 6378.1363; %km from Vallado
    e_earth=0.0033528131; %~ from Vallado
    i = 0;
    while abs(phi_gd-phi_gd_old)>epsilon
        i = i+1;
        phi_gd_old = phi_gd;
        Cbody = Re/(sqrt(1-e_earth^2*sin(phi_gd)^2));
        phi_gd = atan2((rz+Cbody*e_earth^2*sin(phi_gd)),rdsat);
    end
    lat = rad2deg(phi_gd);
    %Find the height above the ellipse
    h_ellps = rdsat./(cosd(lat))-Cbody;
else
    %Code developed in AE 765 using Vallado Algorithm 12
    %Earth Constants
    RE = 6378.137; %km
    eE = 0.081819221456;
    rb = 6356.7516005629; %km
    r_body = RE;
    e_body = eE;
    r_polar = rb;
    if nargin <2
        X = R;
        R = X(1:3,:);
    else
        R = R;
    end

    %Rename for ease of use
    rB = r_body;
    eb = e_body;
    rb = r_polar;


    %Get vector components
    ri = R(1,:);
    rj = R(2,:);
    rk = R(3,:);

    %Find rd_sat
    rd_sat = sqrt(ri.^2+rj.^2);
    %% Algorithm 12
    %Find longitude
    lon = asind(rj./rd_sat); 
    %Find del
    del = atan2(rk,rd_sat);
    %Iterate until a solution is found
    lat = del;
    error = 1;
    while error > 1e-6
        C = rB./(sqrt(1-eb.^2*sin(lat).^2));
        Latnew = atan2((rk+C.*eb.^2.*sin(lat)),rd_sat);
        error = abs(Latnew-lat);
        lat = Latnew;
    end
    lat = rad2deg(lat);
    %Find the height above the ellips
    h_ellps = rd_sat./(cosd(lat))-C;
end