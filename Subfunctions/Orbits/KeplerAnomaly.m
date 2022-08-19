function [Mean_Anom, Mean_Motion, Ecc_Anom, TimeSincePeriapsis_s] = KeplerAnomaly(h_AngularMomentumorStruct,e_Eccentricity,...
    nu_TrueAnomaly,a_SemiMajorAxis,mu)
% [Mean_Anom, Mean_Motion, Ecc_Anom, TimeSincePeriapsis_s] = KeplerAnomaly(h_AngularMomentumorStruct,e_Eccentricity,...
%     nu_TrueAnomaly,a_SemiMajorAxis,mu)
%
%Solve kepler's equation with the given COE's 
if nargin<2
    S = h_AngularMomentumorStruct;
    h = S.h_AngularMomentum;
    e = S.e_Eccentricity;
    a = S.a_SemiMajorAxis;
    nu = S.nu_TrueAnomaly; 
    mu = S.mu; %km^3/s^2
else
    h = h_AngularMomentumorStruct; %km^2/s
    e = e_Eccentricity;
    a = a_SemiMajorAxis; %km
    nu = nu_TrueAnomaly; %Radians
end
%Be able do this for hyperbolic (e>1), parabolic (e=0), or elliptical
Mean_Anom = 0;
Ecc_Anom = 0;
    if e>0&&e<1 %Elliptical orbit
%         fprintf('\n<<<<Elliptical orbit>>>>>\n\n')
        %Calculate initial Anomalies
        Ecc_Anom = acos((e+cos(nu))/(1+e*cos(nu))); %radians Equation 2-9 Vallado
        Mean_Anom = Ecc_Anom-e*sin(Ecc_Anom); %radians Equations 2-6 Vallado
        
    elseif abs(1-e)<1e-3 %Parabolic orbit
%         error('\n<<<<Parabolic orbit>>>>>\n\n')
        
    elseif e>1 %Hyperbolic orbit
%         warning('\n<<<<Hyperbolic orbit>>>>>\n\tEccenctric Anomaly == Hyperbolic Anomaly\n')
        %Find the initial hyperbolic anomaly
        Ecc_Anom = asinh((sin(nu)*sqrt(e^2-1))/(1+e*cos(nu))); %Equation 2-31 Vallado
        
        %Find initial mean anomaly
        Mean_Anom = e*sinh(Ecc_Anom)-Ecc_Anom;
    else %e == 0
        Ecc_Anom = nu;
        Mean_Anom = Ecc_Anom-e*sin(Ecc_Anom);
    end
%Calculate mean motion
Mean_Motion = sqrt(mu/(a)^3);
%Find the time since periapsis
TimeSincePeriapsis_s = Mean_Anom/Mean_Motion;
end