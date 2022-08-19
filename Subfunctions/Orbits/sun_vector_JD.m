function [RS_au, RS_km] = sun_vector_JD(init_utcvec,t_sec)
%
% Function based on Algorithm 29 from Vallado
%
% Created by: Bailey Miller

%Determine the Julian Date for the given initial time and seconds after
yr_init = init_utcvec(1);
mnth_init = init_utcvec(2);
day_init = init_utcvec(3);
hr_init = init_utcvec(4);
min_init = init_utcvec(5);
sec_init = init_utcvec(6)+t_sec;
JD_UTI = UTC2JD(yr_init,mnth_init,day_init,hr_init,min_init,sec_init);

%Apply Algorithm 29
TUTI = (JD_UTI-2451545)/36525; %Julian centuries
lamb_Msun = 280.460+36000.771*TUTI; %degrees
Msun = 357.5291092+35999.05034*TUTI; %degrees
lamb_eclip = lamb_Msun+1.914666471*sind(Msun)+0.019994643*sind(2*Msun); %degrees
rsun = 1.000140612-0.016708617*cosd(Msun)-0.000139589*cosd(2*Msun); %AU
eps = 23.439291 - 0.0130042*TUTI; %degrees

RS_au = [rsun*cosd(lamb_eclip); rsun*cosd(eps)*sind(lamb_eclip); rsun*sind(eps)*sind(lamb_eclip)]; %AU

RS_km = RS_au.*149597871; %km
