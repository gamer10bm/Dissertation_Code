function state = initiate_SC_state(SC)
%  This function initiates the SC state object utilized by other
%  routines as a Matlab structure. This is the place to set new state
%  values, but these could also be added by other routines. Setting 
%  initial values here may result in more contiguous memory usage later.
%  version B.3; saved B. Kaplinger, 3/5/2020
%  Edited by Bailey Miller 2/10/2021

% set state time
state.t = 0;

% set state position and velocity vectors (actual)
state.R = zeros(3,1); state.V = zeros(3,1);

% set state position and velocity vector estimates
state.R_est = zeros(3,1); state.V_est = zeros(3,1);

% set state osculating orbit elements (actual)
state.OE = zeros(6,1);

% set state orbital elements estimate
state.OE_est = zeros(6,1);

% set state attitude (currently represented by inertial to body unit
% quaternion aka Euler parameters) (actual)
state.Q = [0; 0; 0; 1];

% set state attitude estimate
state.Q_est = [0; 0; 0; 1];


% set state current attitude commanded
state.Q_com = [0; 0; 0; 1];

% set state angular velocity
state.W = zeros(3,1);

% set RW spin states in angular velocity units (rad/s) and RPM
state.RW = zeros(SC.N_RW,1); state.RW_RPM = zeros(SC.N_RW,1);

% set misc ADCS parameters
state.RCS_gainP = 4e-4;
state.RCS_gainD = 1.5e-3;

% set battery charge state in Watt Seconds
state.batt = SC.Batt_store_max*0.80; %Assuming some unintentional discharge during transition and launch

end