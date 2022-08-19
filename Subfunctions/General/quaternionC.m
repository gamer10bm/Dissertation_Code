function Q = quaternionC(roll, pitch, yaw, R, V, Q_near)
%  This function converts desired roll(X), pitch(Y), and yaw(Z) angles to a
%  desired quaternion vector. The Euler angles are relative to a local
%  radial, tangential, normal (RTN) coordinate system. In this case, the X
%  axis is treated as along-track(T), Y is cross-track(-N), and Z is
%  nadir(-R)
%  The transformation assumes a rotation order (from RTN) of 1-2-3
% Version A1. Saved 1/23/2020 by B. Kaplinger
    
    % if nearby quaternion is not given for reference
    if (nargin < 6)
        Q_near = [0; 0; 0; 1];
    end

    % find transformation from RTN to inertial
    RTN = zeros(3,3);
    rr = norm(R); RTN(:,3) = -R./rr;
    H = cross(R,V); RTN(:,2) = -H./norm(H);
    RTN(:,1) = cross(RTN(:,2),RTN(:,3));
    
    % Find transformation components from RTN to s/c body frame
    C1 = [1, 0, 0; 0, cos(roll), sin(roll); 0, -sin(roll), cos(roll)];
    C2 = [cos(pitch), 0, -sin(pitch); 0, 1, 0; sin(pitch), 0, cos(pitch)];
    C3 = [cos(yaw), sin(yaw), 0; -sin(yaw), cos(yaw), 0; 0, 0, 1];
    
    % complete transformation from inertial to body frame
    C = C3*C2*C1*RTN';
    
    % convert transformation matrix to quaternion representation
    Q = zeros(4,1);
    Q(4) = 0.5*sqrt(1+C(1,1)+C(2,2)+C(3,3));
    Q(1) = (C(2,3)-C(3,2))/4/Q(4);
    Q(2) = (C(3,1)-C(1,3))/4/Q(4);
    Q(3) = (C(1,2)-C(2,1))/4/Q(4);
    
    % check to see if equivalent quaternion is closer
    n1 = norm(Q_near - Q);
    n2 = norm(Q_near + Q);
    if (n2 < n1)
        Q = -Q;
    end
end
    
    