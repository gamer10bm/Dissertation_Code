function [Rmat] = ECI2ConstFrame(COEstruct)

%This function transforms the coordinates from the ECI frame (orbit
%propagation) to the RSW frame (relative analysis)
%Algorithm based on Pg. 169 of Vallado
%See also RV2coe.m and coe2RV.m


%Grab inputs
Om = COEstruct.RA_RightAscension; %Right Ascension of Ascending Node (RAAN) (rad)
w = COEstruct.w_ArgumentofPerigee; %Argument of Perigee (rad)
i = COEstruct.i_Inclination; %inclination (rad)

%Do precalculations
cO = cos(Om); sO = sin(Om);
cw = cos(w); sw = sin(w);
ci = cos(i); si = sin(i);

%Form rotation matrix

Rmat = [(cO*cw-sO*ci*sw), (sO*cw+cO*ci*sw), (si*sw);
        (-cO*sw-sO*ci*cw), (-sO*cw+cO*ci*cw), (si*cw);
        (sO*si), (-cO*si), ci];
end