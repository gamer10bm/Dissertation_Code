function [SP, illum_logical] = solar_illum(R_sat,R_sun,Re, Q, SC)
%
% Function based on SIGHT Algorithm (35) by Vallado
%
%

%Initialize
illum_logical = false;
%Take dot product between position vectors
RdR = -dot(R_sat,R_sun);

tau_min = (norm(R_sat)^2+RdR)/(norm(R_sat)^2+norm(R_sun)^2+2*RdR);

if tau_min >1 || tau_min <0
    illum_logical = true;
else
    ctaumin = (1-tau_min)*norm(R_sat)^2-RdR*tau_min;
    if ctaumin >= Re^2
        illum_logical = true;
    end
end

%Determine solar power generated
SP = 0;
if illum_logical
    % Determine transformation from body axes to inertial frame
    C = [1-2*(Q(2)^2+Q(3)^2), 2*(Q(1)*Q(2)+Q(3)*Q(4)), 2*(Q(1)*Q(3)-Q(2)*Q(4));
    2*(Q(1)*Q(2)-Q(3)*Q(4)), 1-2*(Q(1)^2+Q(3)^2), 2*(Q(2)*Q(3)+Q(1)*Q(4));
    2*(Q(1)*Q(3)+Q(2)*Q(4)), 2*(Q(2)*Q(3)-Q(1)*Q(4)), 1-2*(Q(1)^2+Q(2)^2)];
    C = C';

    % extract number of s/c sides
    n = length(SC.side);
    
    R_sun_au = R_sun/149597871;
    for i = 1:n
        s = C*SC.side{i}; % convert side normal to inertial coords
        del = dot(s,R_sun_au); % cosine of side normal and sun vector
        if (del > 0) % side is in sunlight
            SP = SP + SC.SP(i)*del; % add power to total with cosine loss
        end
    end % panel loop k
    
end