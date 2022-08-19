function SC = initiate_TX_model()
%  This function initiates the SC vehicle object utilized by other
%  routines as a Matlab structure. This is the place to reset SC global 
%  parameters (or to read them from a config file, for example).
%  Edited by Bailey Miller 2/10/2021


% Physical Problem Parameters
SC.mu = 3.98600435e5;   % standard gravitational parameter 
                        % of problem (km^3/s^2)
SC.Re = 6378.137;       % Equatorial radius of body for problem (km)
SC.J2 = 1.08262668e-3;  % Oblateness parameter for long-term perturbations
SC.SRP = 4.59e-6;       % Solar radiation pressure in Pa

%Drag parameters
SC.CD = 1.28; %Drag Coefficient
SC.A = .36/1000/1000; %km^2 
SC.m = 3.33; %kg
SC.rho0 = 4e-4; %kg/km^3
SC.r0 = 7298145/1000;%km
SC.H = 200; %km
SC.thetadot = 7.29211585530066e-5;%rad/s

% Vehicle Coordinate Parameters
SC.CoM = [50.495; 161.02; 47.97].*1e-3; % Mean CoM (m) from corner
SC.MoI = [37500; 7500; 37500].*1e-6;    % Mean MoI (kg-m^2)

%Side Description
    % 1) 1Ux3U Face (Side)
    % 2) 1Ux3U Face (Side)
    % 3) 1Ux2U Face (Top)
    % 4) 1Ux2U Face (Bottom)
    % 5) 2Ux3U Face (Front)
    % 6) 2Ux3U Face (Back)
% Vehicle Side Unit Vectors
SC.Nsides = 6;
SC.side{1} = [-1; 0; 0];
SC.side{2} = [1; 0; 0];
SC.side{3} = [0; -1; 0];
SC.side{4} = [0; 1; 0];
SC.side{5} = [0; 0; -1];
SC.side{6} = [0; 0; 1];

% Vehicle Side Areas (m^2)
SC.area = zeros(SC.Nsides,1);
SC.area(1) = 300*100*1e-6;
SC.area(2) = 300*100*1e-6;
SC.area(3) = 100*100*1e-6;
SC.area(4) = 100*100*1e-6;
SC.area(5) = 300*100*1e-6;
SC.area(6) = 300*100*1e-6;

% Vehicle Side Centroids (m)
SC.centroid{1} = [0; 150; 50].*1e-3;
SC.centroid{2} = [100; 150; 50].*1e-3;
SC.centroid{3} = [50; 0; 50].*1e-3;
SC.centroid{4} = [50; 300; 50].*1e-3;
SC.centroid{5} = [50; 150; 0].*1e-3;
SC.centroid{6} = [50; 150; 100].*1e-3;

% Vehicle Solar Panel Collection (in W)
SC.SP = zeros(SC.Nsides,1);
SC.SP(1) = 8.4; %Endurosat 3U
SC.SP(2) = 8.4; % %Endurosat 3U
SC.SP(3) = 0;
SC.SP(4) = 0;
SC.SP(5) = 19.2; %Endurosat 6U
SC.SP(6) = 0;

% Vehicle Side Reflectivity
SC.reflect = zeros(SC.Nsides,1);
SC.reflect(1) = 0.5;
SC.reflect(2) = 0.333;
SC.reflect(3) = 0.5;
SC.reflect(4) = 0.5;
SC.reflect(5) = 0.3;
SC.reflect(6) = 0.333;

% SC Antenna Pointing (primary axis)
SC.antenna = SC.side{3};

% SC Antenna Attitude (using primary axis as X) - quaternion
th = pi/2; e = [0; 0; -1];
SC.antennaQ = [e(1)*sin(th/2); e(2)*sin(th/2); e(3)*sin(th/2); cos(th/2)];
SC.antennaQ = SC.antennaQ./norm(SC.antennaQ);

% SC Antenna Parameters
SC.antenna_P_max = 3.5; % maximum transmission power in W

% SC Star Tracker Pointing (primary axis)
SC.star{1} = SC.side{2};
SC.star{2} = SC.side{6};

% SC Star Tracker Attitude (using primary axis as X) - quaternion
th = 0; e = [1; 0; 0];
SC.starQ{1} = [e(1)*sin(th/2); e(2)*sin(th/2); e(3)*sin(th/2); cos(th/2)];
SC.starQ{1} = SC.starQ{1}./norm(SC.starQ{1});
th = pi/2; e = [0; -1; 0];
SC.starQ{2} = [e(1)*sin(th/2); e(2)*sin(th/2); e(3)*sin(th/2); cos(th/2)];
SC.starQ{2} = SC.starQ{2}./norm(SC.starQ{2});

% SC Star Tracker Parameters (FoV assumed Y,Z axes respectively) in deg
SC.Nstar = 2;
SC.star_FoV{1} = [16.8; 12.6];
SC.star_FoV{2} = [16.8; 12.6];
SC.star_P_max = zeros(SC.Nstar,1); % maximum consumed power in W
SC.star_P_max(1) = 0.5;
SC.star_P_max(2) = 0.5;

% SC Reaction Wheel Orientations
SC.N_RW = 3;
SC.RW{1} = [1; 0; 0];
SC.RW{2} = [0; 1; 0];
SC.RW{3} = [0; 0; 1];

% SC Reaction Wheel Parameters
SC.RW_mom_max = zeros(SC.N_RW,1); % maximum momentum storage in N-m-s
SC.RW_mom_max(1) = 10.0e-3; 
SC.RW_mom_max(2) = 10.0e-3; 
SC.RW_mom_max(3) = 10.0e-3; 
SC.RW_T_max = zeros(SC.N_RW,1); % maximum applied Torque N-m
SC.RW_T_max(1) = 3.0e-3; 
SC.RW_T_max(2) = 3.0e-3; 
SC.RW_T_max(3) = 3.0e-3; 
SC.RW_rpm_max = zeros(SC.N_RW,1); % maximum wheel speed in RPM
SC.RW_rpm_max(1) = 10e3; 
SC.RW_rpm_max(2) = 10e3; 
SC.RW_rpm_max(3) = 10e3; 
SC.RW_MoI = zeros(SC.N_RW,1); % wheel moments of inertia in kg-m^2(derived)
SC.RW_MoI(1) = SC.RW_mom_max(1)/(SC.RW_rpm_max(1)*pi/30); 
SC.RW_MoI(2) = SC.RW_mom_max(2)/(SC.RW_rpm_max(2)*pi/30);  
SC.RW_MoI(3) = SC.RW_mom_max(3)/(SC.RW_rpm_max(3)*pi/30);  
SC.RW_P_mean = zeros(SC.N_RW,1); % mean consumed power in W (at 5k RPM)
SC.RW_P_mean(1) = 0.2;
SC.RW_P_mean(2) = 0.2;
SC.RW_P_mean(3) = 0.2;

% SC Torque Rod Orientations
SC.N_TR = 3;
SC.TR{1} = [1; 0; 0];
SC.TR{2} = [0; 1; 0];
SC.TR{3} = [0; 0; 1];

% SC Torque Rod Parameters
SC.TR_D_max = zeros(SC.N_TR,1); % maximum dipole moments in A-m^2
SC.TR_D_max(1) = 0.1;
SC.TR_D_max(2) = 0.1;
SC.TR_D_max(3) = 0.1;
SC.TR_P_max = zeros(SC.N_TR,1); % max consumed power in W 
SC.TR_P_max(1) = 0.75;
SC.TR_P_max(2) = 0.75;
SC.TR_P_max(3) = 0.75;

% SC gyro parameters
SC.gyro_bias_stab = 8.333e-4;  % bias instability in deg/s
SC.gyro_rand_walk = 1.6667e-3;  % random walk in deg/s^(1/2)
SC.gyro_P_max = 0.1;  % max power consumption in W

% SC Battery Parameters
SC.Batt_store_max = 40*3600;   % W-s
SC.Batt_input_max = 40;   % W
SC.Batt_output_max = 80;   % W
end