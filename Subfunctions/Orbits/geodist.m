function [d_km] = geodist(lat1,lon1,lat2,lon2)

Re = 6371; %km
d_km=2*asin(sqrt((sind((lat1-lat2)/2))^2 + cosd(lat1)*cosd(lat2)*(sind((lon1-lon2)/2))^2))*Re;