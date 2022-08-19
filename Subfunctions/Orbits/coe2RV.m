function [X, R, V] = coe2RV(h_AngularMomentumorStruct,e_Eccentricity,i_Inclination,...
    nu_TrueAnomaly,RA_RightAscension,w_ArgumentofPerigee,mu)
% [X, R, V] = coe2RV(h_AngularMomentumorStruct,e_Eccentricity,i_Inclination,...
%     nu_TrueAnomaly,RA_RightAscension,w_ArgumentofPerigee,mu)
%
%
%
if nargin<2
    S = h_AngularMomentumorStruct;
    h = S.h_AngularMomentum;
    a = S.a_SemiMajorAxis;
    e = S.e_Eccentricity;
    i = S.i_Inclination;
    RA = S.RA_RightAscension;
    w = S.w_ArgumentofPerigee;
    nu = S.nu_TrueAnomaly;
    try
        mu = S.mu; %km^3/s^2
    end
elseif nargin<5
    h = h_AngularMomentumorStruct; %km^2/s
    e = e_Eccentricity; 
    i = i_Inclination; %Radians
    nu = nu_TrueAnomaly; %Radians
    %Circular (e=0) equitorial (i=0) orbit
        %Argument of perigee and RA are undefined
    w =0; 
    RA = 0;
else
    h = h_AngularMomentumorStruct; %km^2/s
    e = e_Eccentricity; 
    i = i_Inclination; %Radians
    RA = RA_RightAscension; %Radians
    w = w_ArgumentofPerigee; %Radians
    nu = nu_TrueAnomaly; %Radians
end
%Check if mu exists
if nargin < 7
    mu = 3.986004e5; %km^3/s^2
end

%Check if Right ascension is undefined
if isnan(RA)
    RA = 0;
end

%Check if argument of perigee is undefined
if isnan(w)
    w = 0;
end



%% Build the base state vectors in the perifocal frame

r_coef = h^2/(mu*(1+e*cos(nu)));
r_vec_peri = r_coef*[cos(nu);sin(nu);0];

v_coef = mu/h;
v_vec_peri = v_coef*[-sin(nu);(e+cos(nu)); 0];

%% Build the transformation matrix
cRA = cos(RA);
sRA = sin(RA);
cw = cos(w);
sw = sin(w);
ci = cos(i);
si = sin(i);

%Column 1
A_1 = [(cRA*cw-sRA*sw*ci); (sRA*cw+cRA*ci*sw); (si*sw)];

%Column 2
A_2 = [(-cRA*sw-sRA*ci*cw); (-sRA*sw+cRA*ci*cw); (si*cw)];

%Column 3
A_3 = [(sRA*si); (-cRA*si); ci];

%Whole thing
A = [A_1, A_2, A_3];

%% Multiply this to your perifocal coordinates to get centric coordinates
R = A*r_vec_peri;
V = A*v_vec_peri;
X = [R; V];
end