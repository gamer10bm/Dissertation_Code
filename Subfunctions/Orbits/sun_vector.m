function rS = sun_vector(epoch, t)
% This function determines the unit vector to the Sun from Earth. The
% output coordinates are in Earth-Centered Equatorial Frame.

jd = epoch.JD + t/3600/24;
[~, rE, ~] = planet_elements_and_sv_JD(3, jd);
rS = -rE'./norm(rE);

% vector is currently in ecliptic coordinates, need to rotate by obliquity
% of the ecliptic
ee = 23.43669*pi/180;

C = [1, 0, 0; 0, cos(ee), -sin(ee); 0, sin(ee), cos(ee)];

rS = C*rS;

end