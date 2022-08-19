function [lats, lons] = geocircle(lat0_deg,lon0_deg,radius_km,n)
% [lats, lons] = geocircle(lat0_deg,lon0_deg,radius_km,n)
%
% Inputs: lat0_deg - Latitude in degrees
%             lon0_deg - Longitude in degrees
%             radius_km - Radius of circle in km
%             n - Length of output vectors
%
% Outputs: lats - Latitudes of circle in degrees
%                lons - Longitudes of circle in degrees
%
% Created by Bailey Miller
if nargin<4
    n = 100;
end
Re = 6378.137; %Earth Equatorial Radius (km)
max_ang = atand(radius_km/Re);

%Initialize
lats = zeros(1,n);
lons = zeros(1,n);
for i = 1:n
    ang = (i-1)/n*360; %bearing clockwise from north
    %Determine circle coordinates
    %http://www.movable-type.co.uk/scripts/latlong.html 
    lats(i) = asind(sind(lat0_deg)*cosd(max_ang)+cosd(lat0_deg)*sind(max_ang)*cosd(ang));
    lons(i) = lon0_deg+rad2deg(atan2(sind(ang)*sind(max_ang)*cosd(lat0_deg),cosd(max_ang)-sind(lat0_deg)*sind(lats(i))));
end
%Modulate to keep between -180 and 180
lons(lons>180) = lons(lons>180)-360;
end