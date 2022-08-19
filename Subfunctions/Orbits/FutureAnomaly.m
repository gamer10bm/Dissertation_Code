function [nu_new,COEstruct] = FutureAnomaly(dt_sec,h_AngularMomentumorStruct,e_Eccentricity,...
    nu_TrueAnomaly,a_SemiMajorAxis,mu)
%Solve kepler's equation with the given COE's to get the true anomaly at a
%future time
if nargin<3
    COEstruct = h_AngularMomentumorStruct;
    h = COEstruct.h_AngularMomentum;
    e = COEstruct.e_Eccentricity;
    a = COEstruct.a_SemiMajorAxis;
    nu = COEstruct.nu_TrueAnomaly; 
    mu = COEstruct.mu; %km^3/s^2
else
    h = h_AngularMomentumorStruct; %km^2/s
    e = e_Eccentricity;
    a = a_SemiMajorAxis; %km
    nu = nu_TrueAnomaly; %Radians
end

%Be able do this for hyperbolic (e>1), parabolic (e=1), or elliptical
    if e>=0&&e<1 %Elliptical orbit
%         fprintf('\n<<<<Elliptical orbit>>>>>\n\n')
        %Calculate initial Anomalies
        E_init = acos((e+cos(nu))/(1+e*cos(nu))); %radians Equation 2-9 Vallado
        M_init = E_init-e*sin(E_init); %radians Equations 2-6 Vallado
        
        %Calculate mean motion
        n = sqrt(mu/a^3); %radians/s
        
        %Calculate new Mean anomaly and modulate
        M_end = dt_sec*n+M_init;
        M_end =mod(M_end,2*pi);
        %Use the new Mean anomaly to find the new eccentric anomaly
        E_end = Elliptical_E_Anomaly(M_end,e);
        
        %Calculate the new true anomaly
        nu_new = acos((cos(E_end)-e)/(1-e*cos(E_end))); %radians
        %Do quadrant check
        if E_end >pi
            nu_new = 2*pi-nu_new;
        end
        
    elseif abs(1-e)<1e-3 %Parabolic orbit
        fprintf('\n<<<<Parabolic orbit>>>>>\n\n')
        %Find mean motion
            %Get semilatus rectum
            p = h^2/mu;
            
            %Evaluate
            n_p = 2*sqrt(mu/p^3); %rads/sec
        %Solve Barkers Equation
            %Find s
            s = acot(3/2*n_p*dt_sec)/2;
            
            %Find w
            w = atan((tan(s))^(1/3));
            
            %Find B at the end
            B_end = 2*cot(2*w); %radians
            
        %Find the new true anomaly
        nu_new = atan(B_end)*2;
        
    elseif e>1 %Hyperbolic orbit
        fprintf('\n<<<<Hyperbolic orbit>>>>>\n\n')
        %Find the initial hyperbolic anomaly
        H_init = asinh((sin(nu)*sqrt(e^2-1))/(1+e*cos(nu))); %Equation 2-31 Vallado
        
        %Find initial mean anomaly
        M_init = e*sinh(H_init)-H_init;
        
        %Calculate mean motion
        n = sqrt(mu/(-a)^3);
        
        %Find next mean anomaly
        M_end = dt_sec*n + M_init;
        
        %Solve for end hyperbolic anomaly (Algorithm 4 Vallado)
        H_end = Hyperbolic_H_Anomaly(M_end,e);
        
        %Find the new true anomaly
        nu_new = 2*atan(sqrt((e+1)/(e-1))*tanh(H_end/2));
    else %Uh-oh
        error('Uh-oh, the eccentricity of %.2f you have calculated doesn''t fit our criteria.\n',e)
    end
    %Update COE structure if given
    if exist('S','var')
        COEstruct.nu_TrueAnomaly = nu_new;
    end
end

%% Subfunction section
function E_end = Elliptical_E_Anomaly(M_end,e,tolerance)
%Find the eccentric anomaly for an elliptical orbit
    if nargin<3
        tolerance = 1e-8;
    end
    %Inital guess of E_end
    if (M_end<0 && M_end>-pi) || M_end > pi
        E_end = M_end-e;
    else
        E_end = M_end+e;
    end

    %Iterate until dE is less than the tolerance
    dE = 100; %Initialize dE
    while abs(dE) > tolerance
        f_E = (M_end-E_end+e*sin(E_end)); %Original Equation
        f_E_prime = (1-e*cos(E_end)); %Equation derivative
        dE = f_E/f_E_prime; %Newton-Rhapson Method
        E_end = E_end + dE;
    end
end

%Hyperbolic anomaly
function H_end = Hyperbolic_H_Anomaly(M_end,e,tolerance)
%Find the Hyperbolic anomaly for a hyperbolic orbit
    if nargin<3
        tolerance = 1e-8;
    end
%Solve for end hyperbolic anomaly (Algorithm 4 Vallado)
    %Make an initial guess
    if e < 1.6
        if (M_end<0&&M_end>-pi)||M_end>pi
            H_end = M_end -e;
        else
            H_end = M_end +e;
        end
    else
        if e < 3.6 && abs(M_end)>pi
            H_end = M_end-sign(M_end)*e;
        else
            H_end = M_end/(e-1);
        end
    end
    %Initialize the dH
    dH = 100;
    %Apply Newton Rhapson method
    while abs(dH)>tolerance
        f_H = M_end - e*sinh(H_end)+H_end;
        f_H_prime = e*cosh(H_end)-1;
        dH = f_H/f_H_prime;
        H_end = H_end + dH;
    end
end