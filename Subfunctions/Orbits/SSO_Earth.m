function OE = SSO_Earth(h,SC)
%  This function returns the classical orbital parameters (including an
%  empty spot for true anomaly) of an orbit corresponding to a
%  sun-synchronous orbit around Earth.
%
%  The input "h" is an altitude given in km
%  The output is an array including:
%    OE(1) = specific angular momentum (km^2/s)
%    OE(2) = eccentricity 
%    OE(3) = inclination (rad)
%    OE(4) = Right Ascension of Ascending Node (RAAN) (rad)
%    OE(5) = Argument of Perigee (rad)
%    OE(6) = True Anomaly (rad)
%
%  A 6-element vector output is used to maintain compatibility with other
%  routines. This routine, in its original form, outputs values for a
%  circular orbit, therefore elements 3 and 5 are returned as zeros.
%  Additionally, elements 4 and 6 are returned as zero, as these are
%  undetermined by the given input.
%  version A.0; saved B. Kaplinger, 1/3/2020
%  Edited by Bailey Miller 2/10/2021

if nargin <2 || (~isfield(SC,'mu') && ~isfield(SC,'Re') && ~isfield(SC,'J2'))
    mu = 3.98600435e5; % km^3/s^2
    R = 6378.137; % km
    J2 = 1.08262668e-3;
else
    mu = SC.mu; R = SC.Re; J2 = SC.J2;
end
    P = 365.25*24*3600; %sec Earth Period about Sun
    
    % sun axis rotation rate (rad/s)
    W = 2*pi/(P);

    % convert given altitude to semimajor axis and momentum
    a = R + h;
    H = sqrt(mu*a);
    
    % compute inclination to match sun rate
    cosi = -2*W*a^(7/2)/3/J2/R^2/sqrt(mu);
    i = acos(cosi);
    
    % set output vector
    OE = [H; 0; i; 0; 0; 0];
end